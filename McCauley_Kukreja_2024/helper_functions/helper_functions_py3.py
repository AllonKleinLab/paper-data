import numpy as np
import scipy
import scipy.stats
import scipy.sparse 
from sklearn.decomposition import PCA,TruncatedSVD
from sklearn.neighbors import NearestNeighbors
import time
import os
import json
from datetime import datetime
import matplotlib.pyplot as plt
import pandas as pd
import sklearn.metrics as skm
import seaborn as sn
import scanpy as sc
import scanpy.external as sce
import  statsmodels.sandbox.stats




########## LOADING DATA
def load_genes(filename, delimiter='\t', column=0, skip_rows=0):
    gene_list = []

    with open(filename) as f:
        for iL,l in enumerate(f):
            if iL >= skip_rows:
                gene_list.append(l.strip('\n').split(delimiter)[column])
    return gene_list

def make_genes_unique(orig_gene_list):
    gene_list = []
    gene_dict = {}

    for gene in orig_gene_list:
        if gene in gene_dict:
            gene_dict[gene] += 1
            gene_list.append(gene + '__' + str(gene_dict[gene]))
            if gene_dict[gene] == 2:
                i = gene_list.index(gene)
                gene_list[i] = gene + '__1'
        else:
           gene_dict[gene] = 1
           gene_list.append(gene)
    return gene_list

def load_pickle(fname):
    '''
    Load .pickle data
    '''
    try:
        import pickle
        dat = pickle.load(open(fname, 'rb'))
    except:
        import pandas as pd
        dat = pd.read_pickle(fname)
    return dat

### loading counts

def file_opener(filename):
    if filename.endswith('.gz'):
        fileData = open(filename, 'rb')
        import gzip
        outData = gzip.GzipFile(fileobj = fileData, mode = 'rb')
    elif filename.endswith('.zip'):
        fileData = open(filename, 'rb')
        import zipfile
        zipData = zipfile.ZipFile(fileData, 'r')
        fnClean = filename.strip('/').split('/')[-1][:-4]
        outData = zipData.open(fnClean)
    else:
        fileData = open(filename, 'r')
        outData = fileData
    return outData

def load_mtx(file_data):
    ''' Reads mtx file or, supposedly, an open file object
        Returns scipy.sparse.coo_matrix (if sparse)'''
    return scipy.io.mmread(file_data).tocsc()

def load_npz(file_data):
    return scipy.sparse.load_npz(file_data).tocsc()

def load_npy(file_data):
    return scipy.sparse.csc_matrix(np.load(file_data))

def load_text(file_data, delim='\t', start_column=None, start_row=None, print_row_interval=None):
    X_data = []
    X_row = []
    X_col = []
    
    ncol = None

    for row_ix, dat in enumerate(file_data):
        if print_row_interval is not None:
            if (row_ix+1) % print_row_interval == 0:
                print('Row {}'.format(row_ix+1))
        if type(dat) == bytes:
            dat = dat.decode('utf-8')
        dat = dat.strip('\n').split(delim)

        if start_row is None:
            current_col = 0
            found_float = False
            while not found_float and current_col < len(dat):
                try: 
                    tmp = float(dat[current_col])
                    
                    try:
                        rowdat = np.array(list(map(float, dat[current_col:])))
                        ncol = len(rowdat)
                        col_ix = np.nonzero(rowdat)[0]

                        found_float = True
                        start_row = row_ix
                        start_column = current_col

                        X_col.extend(col_ix)
                        X_row.extend([row_ix - start_row] * len(col_ix))
                        X_data.extend(rowdat[col_ix])

                    except:
                        current_col += 1

                except:
                    current_col += 1
        elif row_ix >= start_row:
            rowdat = np.array(list(map(float, dat[start_column:])))
            if ncol is None:
                ncol = len(rowdat)
            else:
                if len(rowdat) != ncol:
                    return 'ERROR: Rows have different numbers of numeric columns.'
            col_ix = np.nonzero(rowdat)[0]
            X_col.extend(col_ix)
            X_row.extend([row_ix - start_row] * len(col_ix))
            X_data.extend(rowdat[col_ix])

    if start_row is None:
        return 'ERROR: no numeric values found'

    nrow = row_ix - start_row + 1
    E = scipy.sparse.coo_matrix((X_data, (X_row, X_col)), dtype=float, shape=(nrow, ncol)).tocsc()
    
    return E

########## USEFUL SPARSE FUNCTIONS

def sparse_var(E, axis=0):
    ''' variance across the specified axis '''

    mean_gene = E.mean(axis=axis).A.squeeze()
    tmp = E.copy()
    tmp.data **= 2
    return tmp.mean(axis=axis).A.squeeze() - mean_gene ** 2

def sparse_multiply(E, a):
    ''' multiply each row of E by a scalar '''

    nrow = E.shape[0]
    w = scipy.sparse.lil_matrix((nrow, nrow))
    w.setdiag(a)
    return w * E

def sparse_zscore(E, gene_mean=None, gene_stdev=None):
    ''' z-score normalize each column of E '''

    if gene_mean is None:
        gene_mean = E.mean(0)
    if gene_stdev is None:
        gene_stdev = np.sqrt(sparse_var(E))
    return sparse_multiply((E - gene_mean).T, 1/gene_stdev).T


########## CELL FILTERING
def filter_dict(d, filt):
    for k,v in d.items():
        if k != 'meta':
            if len(v.shape) == 1:
                d[k] = v[filt]
            else:
                d[k] = v[filt,:]
    return d

########## GENE FILTERING

def runningquantile(x, y, p, nBins):

    ind = np.argsort(x)
    x = x[ind]
    y = y[ind]

    dx = (x[-1] - x[0]) / nBins
    xOut = np.linspace(x[0]+dx/2, x[-1]-dx/2, nBins)

    yOut = np.zeros(xOut.shape)

    for i in range(len(xOut)):
        ind = np.nonzero((x >= xOut[i]-dx/2) & (x < xOut[i]+dx/2))[0]
        if len(ind) > 0:
            yOut[i] = np.percentile(y[ind], p)
        else:
            if i > 0:
                yOut[i] = yOut[i-1]
            else:
                yOut[i] = np.nan

    return xOut, yOut


def get_vscores(E, min_mean=0, nBins=50, fit_percentile=0.1, error_wt=1):
    '''
    Calculate v-score (above-Poisson noise statistic) for genes in the input counts matrix
    Return v-scores and other stats
    '''

    ncell = E.shape[0]

    mu_gene = E.mean(axis=0).A.squeeze()
    gene_ix = np.nonzero(mu_gene > min_mean)[0]
    mu_gene = mu_gene[gene_ix]

    tmp = E[:,gene_ix]
    tmp.data **= 2
    var_gene = tmp.mean(axis=0).A.squeeze() - mu_gene ** 2
    del tmp
    FF_gene = var_gene / mu_gene

    data_x = np.log(mu_gene)
    data_y = np.log(FF_gene / mu_gene)

    x, y = runningquantile(data_x, data_y, fit_percentile, nBins)
    x = x[~np.isnan(y)]
    y = y[~np.isnan(y)]

    gLog = lambda input: np.log(input[1] * np.exp(-input[0]) + input[2])
    h,b = np.histogram(np.log(FF_gene[mu_gene>0]), bins=200)
    b = b[:-1] + np.diff(b)/2
    max_ix = np.argmax(h)
    c = np.max((np.exp(b[max_ix]), 1))
    errFun = lambda b2: np.sum(abs(gLog([x,c,b2])-y) ** error_wt)
    b0 = 0.1
    b = scipy.optimize.fmin(func = errFun, x0=[b0], disp=False)
    a = c / (1+b) - 1


    v_scores = FF_gene / ((1+a)*(1+b) + b * mu_gene);
    CV_eff = np.sqrt((1+a)*(1+b) - 1);
    CV_input = np.sqrt(b);

    return v_scores, CV_eff, CV_input, gene_ix, mu_gene, FF_gene, a, b

def filter_genes(E, base_ix = [], min_vscore_pctl = 85, min_counts = 3, min_cells = 3, show_vscore_plot = False, sample_name = ''):
    ''' 
    Filter genes by expression level and variability
    Return list of filtered gene indices
    '''

    if len(base_ix) == 0:
        base_ix = np.arange(E.shape[0])

    Vscores, CV_eff, CV_input, gene_ix, mu_gene, FF_gene, a, b = get_vscores(E[base_ix, :])
    ix2 = Vscores>0
    Vscores = Vscores[ix2]
    gene_ix = gene_ix[ix2]
    mu_gene = mu_gene[ix2]
    FF_gene = FF_gene[ix2]
    min_vscore = np.percentile(Vscores, min_vscore_pctl)
    ix = (((E[:,gene_ix] >= min_counts).sum(0).A.squeeze() >= min_cells) & (Vscores >= min_vscore))
    
    if show_vscore_plot:
        import matplotlib.pyplot as plt
        x_min = 0.5*np.min(mu_gene)
        x_max = 2*np.max(mu_gene)
        xTh = x_min * np.exp(np.log(x_max/x_min)*np.linspace(0,1,100))
        yTh = (1 + a)*(1+b) + b * xTh
        plt.figure(figsize=(8, 6));
        plt.scatter(np.log10(mu_gene), np.log10(FF_gene), c = [.8,.8,.8], alpha = 0.3, edgecolors='');
        plt.scatter(np.log10(mu_gene)[ix], np.log10(FF_gene)[ix], c = [0,0,0], alpha = 0.3, edgecolors='');
        plt.plot(np.log10(xTh),np.log10(yTh));
        plt.title(sample_name)
        plt.xlabel('log10(mean)');
        plt.ylabel('log10(Fano factor)');
        plt.show()

    return gene_ix[ix]

def remove_corr_genes(E, gene_list, exclude_corr_genes_list, test_gene_idx, min_corr = 0.1):
    seed_ix_list = []
    for l in exclude_corr_genes_list:
        seed_ix_list.append(np.array([i for i in range(len(gene_list)) if gene_list[i] in l], dtype=int))

    exclude_ix = []
    for iSet in range(len(seed_ix_list)):
        seed_ix = seed_ix_list[iSet][E[:,seed_ix_list[iSet]].sum(axis=0).A.squeeze() > 0]

        tmp = sparse_zscore(E[:,seed_ix])
        tmp = tmp.sum(1).A.squeeze()

        c = np.zeros(len(test_gene_idx))
        for iG in range(len(c)):
            c[iG],_ = scipy.stats.pearsonr(tmp, E[:,test_gene_idx[iG]].A.squeeze())

        exclude_ix.extend([test_gene_idx[i] for i in range(len(test_gene_idx)) if (c[i]) >= min_corr])
    exclude_ix = np.array(exclude_ix)

    return np.array([g for g in test_gene_idx if g not in exclude_ix], dtype=int)


########## CELL NORMALIZATION

def tot_counts_norm(E, exclude_dominant_frac = 1, included = [], target_mean = 0):
    ''' 
    Cell-level total counts normalization of input counts matrix, excluding overly abundant genes if desired.
    Return normalized counts, average total counts, and (if exclude_dominant_frac < 1) list of genes used to calculate total counts 
    '''

    E = E.tocsc()
    ncell = E.shape[0]
    if len(included) == 0:
        if exclude_dominant_frac == 1:
            tots_use = E.sum(axis=1)
        else:
            tots = E.sum(axis=1)
            wtmp = scipy.sparse.lil_matrix((ncell, ncell))
            wtmp.setdiag(1. / tots)
            included = np.asarray(~(((wtmp * E) > exclude_dominant_frac).sum(axis=0) > 0))[0,:]
            tots_use = E[:,included].sum(axis = 1)
            print('Excluded %i genes from normalization' %(np.sum(~included)))
    else:
        tots_use = E[:,included].sum(axis = 1)

    if target_mean == 0:
        target_mean = np.mean(tots_use)

    w = scipy.sparse.lil_matrix((ncell, ncell))
    w.setdiag(float(target_mean) / tots_use)
    Enorm = w * E

    return Enorm.tocsc(), target_mean, included

########## DIMENSIONALITY REDUCTION


def get_pca(E, base_ix=[], numpc=50, keep_sparse=False, normalize=True):
    '''
    Run PCA on the counts matrix E, gene-level normalizing if desired
    Return PCA coordinates
    '''
    # If keep_sparse is True, gene-level normalization maintains sparsity
    #     (no centering) and TruncatedSVD is used instead of normal PCA.

    if len(base_ix) == 0:
        base_ix = np.arange(E.shape[0])

    if keep_sparse:
        if normalize:
            zstd = np.sqrt(sparse_var(E[base_ix,:]))
            Z = sparse_multiply(E.T, 1 / zstd).T
        else:
            Z = E
        pca = TruncatedSVD(n_components=numpc)

    else:
        if normalize:
            zmean = E[base_ix,:].mean(0)
            zstd = np.sqrt(sparse_var(E[base_ix,:]))
            Z = sparse_multiply((E - zmean).T, 1/zstd).T
        else:
            Z = E
        pca = PCA(n_components=numpc)

    pca.fit(Z[base_ix,:])
    return pca.transform(Z)


def preprocess_and_pca(E, total_counts_normalize=True, norm_exclude_abundant_gene_frac=1, min_counts=3, min_cells=5, min_vscore_pctl=85, gene_filter=None, num_pc=50, sparse_pca=False, zscore_normalize=True, show_vscore_plot=False):
    '''
    Total counts normalize, filter genes, run PCA
    Return PCA coordinates and filtered gene indices
    '''

    if total_counts_normalize:
        print('Total count normalizing')
        E = tot_counts_norm(E, exclude_dominant_frac = norm_exclude_abundant_gene_frac)[0]

    if gene_filter is None:
        print('Finding highly variable genes')
        gene_filter = filter_genes(E, min_vscore_pctl=min_vscore_pctl, min_counts=min_counts, min_cells=min_cells, show_vscore_plot=show_vscore_plot)

    print('Using %i genes for PCA' %len(gene_filter))
    PCdat = get_pca(E[:,gene_filter], numpc=num_pc, keep_sparse=sparse_pca, normalize=zscore_normalize)

    return PCdat, gene_filter

########## GRAPH CONSTRUCTION

def get_knn_graph(X, k=5, dist_metric='euclidean', approx=False, return_edges=True):
    '''
    Build k-nearest-neighbor graph
    Return edge list and nearest neighbor matrix
    '''

    t0 = time.time()
    if approx:
        try:
            from annoy import AnnoyIndex
        except:
            approx = False
            print('Could not find library "annoy" for approx. nearest neighbor search')
    if approx:
        #print('Using approximate nearest neighbor search')

        if dist_metric == 'cosine':
            dist_metric = 'angular'
        npc = X.shape[1]
        ncell = X.shape[0]
        annoy_index = AnnoyIndex(npc, metric=dist_metric)

        for i in range(ncell):
            annoy_index.add_item(i, list(X[i,:]))
        annoy_index.build(10) # 10 trees

        knn = []
        for iCell in range(ncell):
            knn.append(annoy_index.get_nns_by_item(iCell, k + 1)[1:])
        knn = np.array(knn, dtype=int)

    else:
        #print('Using sklearn NearestNeighbors')

        if dist_metric == 'cosine':
            nbrs = NearestNeighbors(n_neighbors=k, metric=dist_metric, algorithm='brute').fit(X)
        else:
            nbrs = NearestNeighbors(n_neighbors=k, metric=dist_metric).fit(X)
        knn = nbrs.kneighbors(return_distance=False)

    if return_edges:
        links = set([])
        for i in range(knn.shape[0]):
            for j in knn[i,:]:
                links.add(tuple(sorted((i,j))))

        t_elapse = time.time() - t0
        #print('kNN graph built in %.3f sec' %(t_elapse))

        return links, knn
    return knn

def build_adj_mat(edges, n_nodes):
    A = scipy.sparse.lil_matrix((n_nodes, n_nodes))
    for e in edges:
        i, j = e
        A[i,j] = 1
        A[j,i] = 1
    return A.tocsc()

########## CLUSTERING

def get_spectral_clusters(A, k):
    from sklearn.cluster import SpectralClustering
    spec = SpectralClustering(n_clusters=k, random_state = 0, affinity = 'precomputed', assign_labels = 'discretize')
    return spec.fit_predict(A)


def get_louvain_clusters(nodes, edges):
    import networkx as nx
    import community
    
    G = nx.Graph()
    G.add_nodes_from(nodes)
    G.add_edges_from(edges)
    
    return np.array(list(community.best_partition(G).values()))

########## 2-D EMBEDDING

def get_force_layout(links, n_cells, n_iter=100, edgeWeightInfluence=1, barnesHutTheta=2, scalingRatio=1, gravity=0.05, jitterTolerance=1, verbose=False):
    from fa2 import ForceAtlas2
    import networkx as nx

    G = nx.Graph()
    G.add_nodes_from(range(n_cells))
    G.add_edges_from(list(links))

    forceatlas2 = ForceAtlas2(
                  # Behavior alternatives
                  outboundAttractionDistribution=False,  # Dissuade hubs
                  linLogMode=False,  # NOT IMPLEMENTED
                  adjustSizes=False,  # Prevent overlap (NOT IMPLEMENTED)
                  edgeWeightInfluence=edgeWeightInfluence,

                  # Performance
                  jitterTolerance=jitterTolerance,  # Tolerance
                  barnesHutOptimize=True,
                  barnesHutTheta=barnesHutTheta,
                  multiThreaded=False,  # NOT IMPLEMENTED

                  # Tuning
                  scalingRatio=scalingRatio,
                  strongGravityMode=False,
                  gravity=gravity,
                  # Log
                  verbose=verbose)

    positions = forceatlas2.forceatlas2_networkx_layout(G, pos=None, iterations=n_iter)
    positions = np.array([positions[i] for i in sorted(positions.keys())])
    return positions

def get_umap(X, n_neighbors):
    import umap
    # UMAP: https://github.com/lmcinnes/umap
    embedding = umap.UMAP(n_neighbors=n_neighbors).fit_transform(X)
    return embedding

########## GENE ENRICHMENT

def rank_enriched_genes(E, gene_list, cell_mask, min_counts=3, min_cells=3):
    gix = (E[cell_mask,:]>=min_counts).sum(0).A.squeeze() >= min_cells
    print('%i cells in group' %(sum(cell_mask)))
    print('Considering %i genes' %(sum(gix)))
    
    gene_list = gene_list[gix]
    
    z = sparse_zscore(E[:,gix])
    scores = z[cell_mask,:].mean(0).A.squeeze()
    o = np.argsort(-scores)
    
    return gene_list[o], scores[o]

########## SPRING PREP

def save_hdf5_genes(E, gene_list, filename):
    '''SPRING standard: filename = main_spring_dir + "counts_norm_sparse_genes.hdf5"'''
    
    import h5py
    
    E = E.tocsc()
    
    hf = h5py.File(filename, 'w')
    counts_group = hf.create_group('counts')
    cix_group = hf.create_group('cell_ix')

    hf.attrs['ncells'] = E.shape[0]
    hf.attrs['ngenes'] = E.shape[1]

    for iG, g in enumerate(gene_list):
        counts = E[:,iG].A.squeeze()
        cell_ix = np.nonzero(counts)[0]
        counts = counts[cell_ix]
        counts_group.create_dataset(g, data = counts)
        cix_group.create_dataset(g, data = cell_ix)

    hf.close()
    
def save_hdf5_cells(E, filename):
    '''SPRING standard: filename = main_spring_dir + "counts_norm_sparse_cells.hdf5" '''
    import h5py
    
    E = E.tocsr()
    
    hf = h5py.File(filename, 'w')
    counts_group = hf.create_group('counts')
    gix_group = hf.create_group('gene_ix')

    hf.attrs['ncells'] = E.shape[0]
    hf.attrs['ngenes'] = E.shape[1]

    for iC in range(E.shape[0]):
        counts = E[iC,:].A.squeeze()
        gene_ix = np.nonzero(counts)[0]
        counts = counts[gene_ix]
        counts_group.create_dataset(str(iC), data = counts)
        gix_group.create_dataset(str(iC), data = gene_ix)

    hf.close()
    
def save_sparse_npz(E, filename, compressed = False):
    ''' SPRING standard: filename = main_spring_dir + "/counts_norm.npz"'''
    E = E.tocsc()
    scipy.sparse.save_npz(filename, E, compressed = compressed)

def write_graph(filename, n_nodes, edges):
    nodes = [{'name':int(i),'number':int(i)} for i in range(n_nodes)]
    edges = [{'source':int(i), 'target':int(j), 'distance':0} for i,j in edges]
    out = {'nodes':nodes,'links':edges}
    open(filename,'w').write(json.dumps(out,indent=4, separators=(',', ': ')))

def write_edges(filename, edges):
    with open(filename, 'w') as f:
        for e in edges:
            f.write('%i;%i\n' %(e[0], e[1]))

def write_color_tracks(ctracks, fname):
    out = []
    for name,score in ctracks.items():
        line = name + ',' + ','.join(['%.3f' %x for x in score])
        out += [line]
    out = sorted(out,key=lambda x: x.split(',')[0])
    open(fname,'w').write('\n'.join(out))

def frac_to_hex(frac):
    rgb = tuple(np.array(np.array(plt.cm.jet(frac)[:3])*255,dtype=int))
    return '#%02x%02x%02x' % rgb

def get_color_stats_genes(color_stats, E, gene_list):
    means = E.mean(0).A.squeeze()
    stdevs = np.sqrt(sparse_var(E, 0))
    mins = E.min(0).todense().A1
    maxes = E.max(0).todense().A1
    
    pctl = 99.6
    pctl_n = (100-pctl) / 100. * E.shape[0]    
    pctls = np.zeros(E.shape[1], dtype=float)
    for iG in range(E.shape[1]):
        n_nonzero = E.indptr[iG+1] - E.indptr[iG]
        if n_nonzero > pctl_n:
            pctls[iG] = np.percentile(E.data[E.indptr[iG]:E.indptr[iG+1]], 100 - 100 * pctl_n / n_nonzero)
        else:
            pctls[iG] = 0
        color_stats[gene_list[iG]] = tuple(map(float, (means[iG], stdevs[iG], mins[iG], maxes[iG], pctls[iG])))
    return color_stats

def get_color_stats_custom(color_stats, custom_colors):
    for k,v in custom_colors.items():
        color_stats[k] = tuple(map(float, (np.mean(v),np.std(v),np.min(v),np.max(v),np.percentile(v,99))))
    return color_stats

def save_color_stats(filename, color_stats):
    with open(filename,'w') as f:
        f.write(json.dumps(color_stats, indent=4, sort_keys=True))

def build_categ_colors(categorical_coloring_data, cell_groupings):
    for k,labels in cell_groupings.items():
        label_colors = {l:frac_to_hex(float(i)/len(set(labels))) for i,l in enumerate(list(set(labels)))}
        categorical_coloring_data[k] = {'label_colors':label_colors, 'label_list':labels}
    return categorical_coloring_data

def save_cell_groupings(filename, categorical_coloring_data):
    with open(filename,'w') as f:
        #f.write(json.dumps(categorical_coloring_data,indent=4, sort_keys=True).decode('utf-8'))
        f.write(json.dumps(categorical_coloring_data, indent=4, sort_keys=True))

def save_spring_dir_sparse_hdf5(E,gene_list,project_directory, edges, custom_colors={}, cell_groupings={}):

    if not os.path.exists(project_directory):
        os.makedirs(project_directory)

    if not project_directory[-1] == '/': 
        project_directory += '/'

    # save custom colors
    custom_colors['Uniform'] = np.zeros(E.shape[0])
    write_color_tracks(custom_colors, project_directory+'color_data_gene_sets.csv')

    # create and save a dictionary of color profiles to be used by the visualizer
    color_stats = {}
    color_stats = get_color_stats_genes(color_stats, E, gene_list)
    color_stats = get_color_stats_custom(color_stats, custom_colors)
    save_color_stats(project_directory + 'color_stats.json', color_stats)

    # save cell labels
    categorical_coloring_data = {}
    categorical_coloring_data = build_categ_colors(categorical_coloring_data, cell_groupings)
    save_cell_groupings(project_directory+'categorical_coloring_data.json', categorical_coloring_data)

    # write graph
    write_graph(project_directory + 'graph_data.json', E.shape[0], edges)
    write_edges(project_directory + 'edges.csv', edges)


#========================================================================================#

def make_spring_subplot(E, gene_list, save_path, base_ix = None, normalize = True, 
    exclude_dominant_frac = 1.0, min_counts = 3, min_cells = 5, min_vscore_pctl = 75,show_vscore_plot = False, exclude_gene_names = None, num_pc = 30, sparse_pca = False, pca_norm = True, k_neigh = 4, cell_groupings = {}, num_force_iter = 100, output_spring = True, precomputed_pca = None, gene_filter = None, custom_colors = {}, exclude_corr_genes_list = None, exclude_corr_genes_minCorr = 0.2, dist_metric = 'euclidean', use_approxnn=False, run_doub_detector = False, dd_k=50, dd_frac=5, dd_approx=True, tot_counts_final = None):
    
    out = {}

    E = E.tocsc()
    if base_ix is None:
        base_ix = np.arange(E.shape[0])

    # total counts normalize
    if tot_counts_final is None:
        tot_counts_final = E.sum(1).A.squeeze()
    out['tot_counts_final'] = tot_counts_final

    if normalize:
        #print 'Normalizing'
        E = tot_counts_norm(E, exclude_dominant_frac = exclude_dominant_frac)[0]

    if precomputed_pca is None:
        if gene_filter is None:
            # Get gene stats (above Poisson noise, i.e. V-scores)
            #print 'Filtering genes'
            if (min_counts > 0) or (min_cells > 0) or (min_vscore_pctl > 0): 
                gene_filter = filter_genes(E, base_ix, min_vscore_pctl=min_vscore_pctl, min_counts=min_counts,min_cells=min_cells,show_vscore_plot = show_vscore_plot)
            else:
                gene_filter = np.arange(E.shape[1])

            if len(gene_filter) == 0:
                print('Error: No genes passed filter')
                sys.exit(2)
                #print 'Error: All genes have mean expression < '+repr(min_exp) + ' or CV < '+repr(min_cv)
            #print 'Using %i genes' %(len(gene_filter))

            if not exclude_corr_genes_list is None:
                gene_filter = remove_corr_genes(E, gene_list, exclude_corr_genes_list, gene_filter, min_corr = exclude_corr_genes_minCorr)
                if len(gene_filter) == 0:
                    print('Error: No genes passed filter')
                    sys.exit(2)

            # Remove user-excluded genes from consideration
            if not exclude_gene_names is None:
                keep_ix = np.array([ii for ii,gix in enumerate(gene_filter) if gene_list[gix] not in exclude_gene_names])
                #print 'Excluded %i user-provided genes' %(len(gene_filter)-len(keep_ix))
                gene_filter = gene_filter[keep_ix]
                if len(gene_filter) == 0:
                    print('Error: No genes passed filter')
                    sys.exit(2)

        out['gene_filter'] = gene_filter
        # RUN PCA
        # if method == 'sparse': normalize by stdev
        # if method == anything else: z-score normalize
        #print 'Running PCA'
        num_pc = min(len(gene_filter), num_pc)
        out['num_pc'] = num_pc
        Epca = get_pca(E[:,gene_filter], base_ix=base_ix, numpc=num_pc, keep_sparse=sparse_pca, normalize = pca_norm)
        out['Epca'] = Epca
    # else:
    #     print 'Using user-supplied PCA coordinates'
    #     Epca = precomputed_pca

    #print 'Building kNN graph'

    links, knn_graph = get_knn_graph(Epca, k=k_neigh, dist_metric = dist_metric, approx=use_approxnn)
    out['knn_graph'] = knn_graph

    if run_doub_detector:
        import doublet_detector as woublet
        #print 'Running woublet'
        doub_score, doub_score_full, doub_labels = woublet.detect_doublets([], counts=tot_counts_final, doub_frac=dd_frac, k=dd_k, use_approxnn=dd_approx, precomputed_pca=Epca)
        out['doub_score'] = doub_score
        out['doub_score_sim'] = doub_score_sim

    if output_spring:

        if not os.path.exists(save_path):
            os.makedirs(save_path)

        #print 'Saving SPRING files to %s' %save_path
        custom_colors['Total Counts'] = tot_counts_final
        np.savez_compressed(save_path + '/intermediates.npz', Epca = Epca, gene_filter = gene_filter, total_counts = tot_counts_final)

        if run_doub_detector:
            custom_colors['Doublet Score'] = doub_score

        if len(cell_groupings) > 0:
            save_spring_dir_sparse_hdf5(E, gene_list, save_path, list(links),
                            custom_colors = custom_colors,
                            cell_groupings = cell_groupings)
        else:
            save_spring_dir_sparse_hdf5(E, gene_list, save_path, list(links),
                            custom_colors = custom_colors)


    if num_force_iter > 0:
        positions = get_force_layout(links, Epca.shape[0], n_iter=num_force_iter, 
            edgeWeightInfluence=1, barnesHutTheta=2, scalingRatio=1, gravity=0.05, 
            jitterTolerance=1, verbose=False)
        positions = positions / 5.0
        positions = positions - np.min(positions, axis = 0) - np.ptp(positions, axis = 0) / 2.0
        positions[:,0] = positions[:,0]  + 750
        positions[:,1] = positions[:,1]  + 250
        out['coordinates'] = positions

    if output_spring:
        if num_force_iter > 0:
            np.savetxt(save_path + '/coordinates.txt',
                       np.hstack((np.arange(positions.shape[0])[:,None], positions)), fmt='%i,%.5f,%.5f')
         
        info_dict = {}
        info_dict['Date'] = '%s' %datetime.now()
        info_dict['Nodes'] = Epca.shape[0]
        info_dict['Filtered_Genes'] = len(gene_filter)
        info_dict['Gene_Var_Pctl'] = min_vscore_pctl
        info_dict['Min_Cells'] = min_cells
        info_dict['Min_Counts'] = min_counts
        info_dict['Num_Neighbors'] = k_neigh
        info_dict['Num_PCs'] = num_pc
        info_dict['Num_Force_Iter'] = num_force_iter
        with open(save_path+'/run_info.json','w') as f:
            #f.write(json.dumps(info_dict,indent=4, sort_keys=True).decode('utf-8'))
            f.write(json.dumps(info_dict, indent=4, sort_keys=True))
         
    return out

#========================================================================================#


############# PLOTTING

def gene_plot(x, y, E, gene_list, gene_name, col_range=(0,100), order_points=False, x_buffer=0, y_buffer=0,
        fig_size=(5,5), point_size=15, colormap='Reds', bg_color=[1,1,1], ax='', smooth_operator = []):
    '''
    Plot gene expression values on a scatter plot.

    Input
        x : x coordinates for scatter plot
        y : y coordinates for scatter plot
        E : gene expression counts matrix (cells x genes)
        gene_list (list of strings, length=n_cells): full list of gene names
        gene_name (string): name of gene to visualize
        col_range (float tuple, length=2): (color_floor, color_ceiling) percentiles
        order_points (boolean): if True, plot points with higher color values on top of points with lower values
        x_buffer (float): white space to add to x limits
        y_buffer (float): white space to add to y limits
        fig_size (float tuple, length=2): size of figure
        point_size (float): size of scatter plot points
        colormap: color scheme for coloring the scatter plot
        bg_color (RGB/HEX/color name): background color

    Output
        fig: figure handle
        ax: axis handle
        pl: scatter plot handle
    '''
    # get gene index and color data
    import matplotlib.pyplot as plt

    gene_ix = gene_list.index(gene_name)
    colordat = E[:,gene_ix].toarray()[:,0]

    if len(smooth_operator) > 0:
        colordat = np.dot(smooth_operator, colordat)

    # get min and max color values
    cmin = np.percentile(colordat, col_range[0])
    cmax = np.percentile(colordat, col_range[1])
    if cmax == 0:
        cmax = max(colordat)

    # order points by intensity, if desired
    if order_points:
        plot_ord = np.argsort(colordat)
    else:
        plot_ord = np.arange(len(colordat))

    # make the plot
    return_all = False
    if ax == '':
        return_all = True
        fig, ax = plt.subplots(1, 1, figsize = fig_size)

    pl = ax.scatter(x[plot_ord], y[plot_ord], c=colordat[plot_ord], s=point_size, edgecolor='none',
                    cmap=colormap, vmin=cmin, vmax=cmax)

    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlim((min(x) - x_buffer, max(x) + x_buffer))
    ax.set_ylim((min(y) - y_buffer, max(y) + y_buffer))
    ax.patch.set_color(bg_color)

    if return_all:
        return fig, ax, pl
    else:
        return pl

########## PLOTTING STUFF

def darken_cmap(cmap, scale_factor):
    cdat = np.zeros((cmap.N, 4))
    for ii in range(cdat.shape[0]):
        curcol = cmap(ii)
        cdat[ii,0] = curcol[0] * scale_factor
        cdat[ii,1] = curcol[1] * scale_factor
        cdat[ii,2] = curcol[2] * scale_factor
        cdat[ii,3] = 1
    cmap = cmap.from_list(cmap.N, cdat)
    return cmap

def custom_cmap(rgb_list):
    import matplotlib.pyplot as plt
    rgb_list = np.array(rgb_list)
    cmap = plt.cm.Reds
    cmap = cmap.from_list(rgb_list.shape[0],rgb_list)
    return cmap

def plot_groups(x, y, groups, lim_buffer = 50, saving = False, fig_dir = './', fig_name = 'fig', res = 300, close_after = False, title_size = 12, point_size = 3, ncol = 5):
    import matplotlib.pyplot as plt

    n_col = int(ncol)
    ngroup = len(np.unique(groups))
    nrow = int(np.ceil(ngroup / float(ncol)))
    fig = plt.figure(figsize = (14, 3 * nrow))
    for ii, c in enumerate(np.unique(groups)):
        ax = plt.subplot(nrow, ncol, ii+1)
        ix = groups == c

        ax.scatter(x[~ix], y[~ix], s = point_size, c = [.8,.8,.8], edgecolors = '')
        ax.scatter(x[ix], y[ix], s = point_size, c = [0,0,0], edgecolors = '')
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlim([min(x) - lim_buffer, max(x) + lim_buffer])
        ax.set_ylim([min(y) - lim_buffer, max(y) + lim_buffer])

        ax.set_title(str(c), fontsize = title_size)

    fig.tight_layout()

    if saving:
        if not os.path.exists(fig_dir):
            os.makedirs(fig_dir)
        plt.savefig(fig_dir + '/' + fig_name + '.png', dpi=res)

    if close_after:
        plt.close()


############ Kalki's functions ###################

# genes correlated a list of genes

def corr_coeff(adata, program_genes, cutoff):
    
    '''
    -------------------------------------------------------------------
    arguments
    
    adata         : adata with only highly variable genes
    program_genes : list of genes in a program (like cell cycle genes)
    cutoff        : cutoff for correlation coefficient
    -------------------------------------------------------------------
    returns
    
    corr_coeff    : correlation coefficient dataframe
    corr_genes    : list of genes correlated to the program genes with 
                    a corr_coeff of greater than cutoff
    -------------------------------------------------------------------- 
    '''
    
    # genes X cell matrix for variable genes
    gene_all = adata.X.T # dimensions - Ga x N
    
    #genes X cell matrix for cell cycle gene list
    gene_subset = adata.raw.X[:,np.in1d(adata.raw.var_names,program_genes)].T # Gb x N
    
    # find correlation coeff:
    centered_subset = gene_subset - gene_subset.mean(axis = 1) # x - mean_x
    centered_subset = np.array(centered_subset[:,np.newaxis,:]) # numpy broadcasting
                                                                # dimensions: Ga x 1 x N

    centered_all = gene_all - gene_all.mean(axis = 1) # y - mean_y
    centered_all = np.array(centered_all[np.newaxis,:,:]) # numpy broadcasting 
                                                          # dimensions: 1 x Gb x N
    
    # get numerator of the pearson_coorelation: 
    # Ga x 1 x N * 1 x Gb x N = Ga x Gb x N
    # sum along z axis: Ga x Gb
    corr_num = (centered_subset*centered_all).sum(axis = 2) 

    # get denominator of the pearson correlation term:
    corr_denom = np.sqrt(np.square(centered_subset).sum(axis=2)*np.square(centered_all).sum(axis=2))

    # find correlation coefficient:
    corr_coeff = corr_num/corr_denom
    
    corr_coeffpd = pd.DataFrame(corr_coeff,
                            index = np.array(program_genes)[np.in1d(program_genes,adata.raw.var_names)],
                           columns= adata.var_names)
    
    # get list of genes highly correlated to the given list:
    corr_genes = np.array(adata.var_names)[abs(corr_coeffpd).max() > cutoff]
 
    return corr_coeffpd, corr_genes
    
    
# baysian classifier

def bay_classifier(obs_counts, state_centroids):
    # obs_counts -> observed counts. it is the cellxgene matrix of your counts
    # state_centroids -> it is a state x gene np matrix
    
    state_prob = state_centroids/state_centroids.sum(axis=1) 
    #state_probabilities which is centroids normalized to 1 so they act as probabilities.
    
    bay = obs_counts.dot(np.log10(state_prob).T)
    
    return bay
    

# differential gene expression analysis between 2 sets of samples

def dge(dataset1, dataset2, alpha, pseudocount):

    '''
    Input:
    --------
    dataset1 : anndata for control data
    dataset2 : anndata for perturbed data
    alpha    : p-value threshold for significance
    
    Returns:
    --------
    log_fold_change      : log2 fold difference between perturbed and control
    pvals                : pvals after multiple hypothesis correction
    mask_differential    : get mask for differentially expressed genes
    diff_expressed_genes : list of differentially expressed genes
    
    '''
    import pandas as pd
    import scipy
    import numpy as np
    import statsmodels.sandbox.stats.multicomp
    
    
    # find mean and log fold change
    means1 = dataset1.X.mean(0).A.squeeze() + pseudocount # control mean
    means2 = dataset2.X.mean(0).A.squeeze() + pseudocount # perturbed mean
    
    log_fold_change = np.log2(means2) - np.log2(means1) # perturbed - control
    
    # find pvalue:                                           
    pvals = np.zeros(dataset1.shape[1])
    for ii in range(dataset1.shape[1]):
        pvals[ii] = scipy.stats.ranksums(dataset1.X[:,ii].A.squeeze(), dataset2.X[:,ii].A.squeeze())[1]

        
    pvals = statsmodels.sandbox.stats.multicomp.multipletests(pvals, alpha=0.05, method='fdr_bh',)[1]
  
    
    # find differentially expressed genes:
    mask1 = log_fold_change  > 1 # greater than 2 fold difference
    mask2 = log_fold_change  < -1 # less than 2 fold difference
    mask3 = pvals < alpha 

    mask_downregulated = mask1*mask3
    mask_upregulated = mask2*mask3

    mask_differential = mask_downregulated + mask_upregulated
    
    diff_expressed_genes = dataset1.var_names[mask_differential]
    
    return {'FC':log_fold_change, 'p_value':pvals, 
    'mask': mask_differential, 'diff_genes':diff_expressed_genes}


# differential gene expression analysis between different sets of samples when the treatments are not coupled
def dge_replicate_unlabeled(dataset1, dataset2, no_of_replicates, alpha):
    '''
    Input:
    --------
    dataset1         : anndata for control data
    dataset2         : anndata for perturbed data
    alpha            : p-value threshold for significance
    no_of_replicates : total number of replicates
    
    Returns:
    --------
    log_fold_change      : log2 fold difference between perturbed and control
    pvals                : pvals after multiple hypothesis correction
    mask_differential    : get mask for differentially expressed genes
    diff_expressed_genes : list of differentially expressed genes
    
    '''
    
    rep_dict_control = {}
    rep_dict_perturbed = {}
    dge_dict = {}
    rep_array = (np.array(range(no_of_replicates)) + 1).astype(str)
    for replicate in rep_array:
        rep_dict_control[replicate] = dataset1[dataset1.obs['replicate'] == replicate]
        rep_dict_perturbed[replicate] = dataset2[dataset2.obs['replicate'] == replicate]


    for i in range(1,no_of_replicates+1):
        j = i
        while j<=no_of_replicates:
            dge_dict['{}_{}'.format(i,j)] = dge(rep_dict_control['{}'.format(i)], rep_dict_perturbed['{}'.format(i)], alpha)
            print ('{}_{} dge done'.format(i,j))
            j = j +1
    dge_dict['all']=dge(dataset1,dataset2,alpha)
    print ('all done')
    
    return dge_dict
    
# differential gene expression analysis between paired conditions

def dge_replicate(dataset1, dataset2, no_of_replicates, alpha):

    '''
    Input:
    --------
    dataset1         : anndata for control data, should have replicate annotations in the obs
    dataset2         : anndata for perturbed data, should have replicate annotations in obs
    alpha            : p-value threshold for significance
    no_of_replicates : total number of replicates
    
    Returns:
    --------
    log_fold_change      : log2 fold difference between perturbed and control
    pvals                : pvals after multiple hypothesis correction
    mask_differential    : get mask for differentially expressed genes
    diff_expressed_genes : list of differentially expressed genes
    
    '''

    rep_dict_control = {}
    rep_dict_perturbed = {}
    dge_dict = {}
    rep_array = (np.array(range(no_of_replicates)) + 1).astype(str)
    for replicate in rep_array:
        rep_dict_control[replicate] = dataset1[dataset1.obs['replicate'] == replicate]
        rep_dict_perturbed[replicate] = dataset2[dataset2.obs['replicate'] == replicate]


    for i in range(1,no_of_replicates+1):
        dge_dict['{}_{}'.format(i,i)] = dge(rep_dict_control['{}'.format(i)], rep_dict_perturbed['{}'.format(i)], alpha)
        print ('{}_{} dge done'.format(i,i))
        

    return dge_dict


# function to find marker genes:

def find_markers(norm_counts, gene_list, groups, min_frac_expr=0.01, min_fold_change=2, pseudocount=0.1, max_p=0.05):
    '''
    norm_counts: normalized counts matrix (scipy.sparse with shape (n_cells, n_genes))
    gene_list: numpy array of gene names (length = n_genes)
    groups: 1-D numpy array with group labels for each cell (shape (n_cells,))
    min_frac_expr: only test genes expressed in this fraction of cells in the cluster of interest
    min_fold_change: min fold-change between highest avg expression and 2nd-highest avg expression
    pseudocount: pseudocount to use when calculating fold-change
    max_p: maximum multiple hypothesis-corrected p-value
    '''
    import pandas as pd, scipy
    import numpy as np
    import statsmodels.sandbox.stats.multicomp
    
    ##########
    min_fold_change = np.log2(min_fold_change)

    # cluster labels for each cell
    clusts_use = groups.copy() 

    # initialize results
    results = {'group': [], 'gene': [], 'p-value': [], 'log2_fold_change': []}

    for c1 in np.unique(clusts_use):

        #############    SETUP 

        # Find cells in test cluster which is the cluster whose marker genes we need to find
        i1 = clusts_use == c1 # give indices of cells with a particular cluster
        i2 = ~i1 # indices of all other cells


        # Find genes expressed in some percent of cells in test cluster
        # basically you get the genes boolean mask with fractional expression value greater than min_frac_expression
        gix = (norm_counts.A[i1,:]>0).sum(0)/i1.sum() > min_frac_expr
        print('Testing {} genes for cluster {}'.format(gix.sum(), c1))


        #############    P-VALUES


        E1 = norm_counts.A[i1,:][:,gix] # expression profile of cells in test cluster
        E2 = norm_counts.A[i2,:][:,gix] # expression profile of cells in all other clusters
        pv = np.zeros(gix.sum())  # initialize p-value matrix
        for ii,iG in enumerate(np.nonzero(gix)[0]): 
            pv[ii] = scipy.stats.ranksums(E1[:,ii], E2[:,ii])[1] # compare genes expression in two clusters
        pv = statsmodels.sandbox.stats.multicomp.multipletests(pv, alpha=0.05, method='fdr_bh',)[1] # to get fdr from p-value

        del E1,E2 # to that the computer memory is not eaten up as this code runs


        #############    HIGHEST VS. 2ND-HIGHEST MEANS


        m1 = norm_counts.A[i1,:][:,gix].mean(0) # mean expression of genes for the test cluster
        m2 = np.zeros(gix.sum()) # initializing mean gene expression matrix for maximum expression of the gene from other cluster
        for c in np.unique(clusts_use):
            if c != c1:
                i_test = clusts_use == c # get cell boolean mask for the cluster to which comparison is being made
                m_test = norm_counts.A[i_test,:][:,gix].mean(0) #  mean gene expression of cells in the cluster
                gix_change = m_test-m2>0   # for updating m2 (if new cluster's genes are higher in expression than the previous cluster) 
                m2[gix_change] = m_test[gix_change] # update m2

        logFc = np.log2((m1+pseudocount)/(m2+pseudocount)) # log fold change of your cluster vs maximum expression from other clusters

        #############    IDENTIFY SIGNIFICANT MARKER GENES
        hi = (pv < max_p) & (logFc > min_fold_change) # for genes significantly over 2 fold from their current cluster
        if hi.sum() > 0:
            results['group'].extend([c1 for i in range(hi.sum())]) 
            results['gene'].extend(list(gene_list[gix][hi]))
            results['p-value'].extend(list(pv[hi]))
            results['log2_fold_change'].extend(list(logFc[hi]))

        print('cluster={}, n_cells={}, n_diff={}'.format(c1, i1.sum(), hi.sum()))

    results = pd.DataFrame(results)[['group', 'gene', 'p-value', 'log2_fold_change']]
    return results

    
    ########### classifiers #########

def make_training(adata, fold, rand_seed=0):

    training_set_size = int(adata.shape[0]-(adata.shape[0]/fold))
    np.random.seed(rand_seed)
    training_index = np.random.choice( np.array(range(adata.shape[0]), dtype = int) , size =  training_set_size, replace = False)
    adata_training = adata[training_index]    
    adata_test = adata_test = adata[list(set(np.array(range(adata.shape[0]), dtype = int))-set(training_index))]
        
    
    return  adata_training, adata_test

# SVM parameters : 'https://medium.com/all-things-ai/in-depth-parameter-tuning-for-svc-758215394769'

def multi_class_svm(train_x, train_y, test_x, C_param=1,n_iterations = 1000):
    """
    Trains a linear SVM for multiclass classifciation using a one-vs-rest strategy

    Args:
        train_x - (n, d) NumPy array (n datapoints each with d features)
        train_y - (n, ) NumPy array containing the labels (int) for each training data point
        test_x - (m, d) NumPy array (m datapoints each with d features)
    Returns:
        pred_test_y - (m,) NumPy array containing the labels (int) for each test data point
    """
    
    
    from sklearn.svm import LinearSVC

    clf = LinearSVC(C= C_param, max_iter = n_iterations)
    clf.fit(train_x, train_y)
    pred_test_y = clf.predict(test_x)
    return pred_test_y

    raise NotImplementedError

def softmax_regression(train_x, train_y, test_x,  c_param=1, n_iterations=1000, seed=0):

    from sklearn.linear_model import LogisticRegression

    lr = LogisticRegression(C=c_param, max_iter = n_iterations, random_state=seed)
    lr.fit(train_x, train_y)
    return lr.predict(test_x)


    
def get_confusion(true_states, predicted_states, state_labels, figuresize = (15,10), linewidth=.2 ):
    confusion = skm.confusion_matrix(true_states, predicted_states,
                                     labels = state_labels)
    f = plt.figure(figsize = figuresize)
    sn.heatmap(confusion/confusion.sum(axis = 1)[:,np.newaxis], 
               xticklabels=state_labels, yticklabels=state_labels, linewidths = linewidth)
    
    return confusion


    
def model_accuracy(confusion):
    return np.diagonal(confusion).sum()/confusion.sum().sum()

    
def model_precision(confusion): 
    return np.diag(confusion)/np.sum(confusion, axis = 0)

def model_recall(confusion):
    return np.diag(confusion)/np.sum(confusion, axis = 1)

def state_sizes(adata, obs_label, figuresize = (5,7), plot = True):
    groups = pd.DataFrame(adata.obs[obs_label]).groupby(obs_label).groups

    statesize = {}
    for keys, values in groups.items():
        statesize[keys] = len(values)

    if plot == True:


        f = plt.figure(figsize = figuresize)

        plt.scatter(statesize.values(),statesize.keys())
    return statesize

def make_training_f1(adata, fold, rand_seed = 0 ):
    
    adata_training, adata_test = make_training(adata, fold, rand_seed)
    
    # use highly variable genes:
    filter_result = sc.pp.filter_genes_dispersion(

    adata_training.X, min_mean=0.0125, max_mean=3, min_disp=0.5)
    print('{} genes passing filter'.format(filter_result['gene_subset'].sum()))

    sc.pl.filter_genes_dispersion(filter_result)
    adataf1_training = adata_training[:,np.array(pd.DataFrame(filter_result)['gene_subset'])]
    adataf1_test = adata_test[:,np.array(pd.DataFrame(filter_result)['gene_subset'])]
    
    return adataf1_training, adataf1_test

def get_PCA_coordinates(adata_training, adata_test, ncomponents = 20):
    #### Total count normalization

    sc.pp.normalize_per_cell(adata_training, counts_per_cell_after=adata_training.obs['n_counts'].mean())

    #### Store unfiltered data as `adata_training_dan.raw` - used in later steps
    adata_training.raw = adata_training

    #### Remove lowly expressed genes - filter genes mask 1
    filter_low_express = sc.pp.filter_genes(adata_training.X, min_cells=3)
    filter_genes_mask1 = filter_low_express[0]
    adata_training = adata_training[:,filter_genes_mask1]

    #### Find variable genes - filter genes mask 2

    filter_result = sc.pp.filter_genes_dispersion(
        adata_training.X, min_mean=0.0125, max_mean=3, min_disp=0.5)

    print('{} genes passing filter'.format(filter_result['gene_subset'].sum()))

    sc.pl.filter_genes_dispersion(filter_result)
    filter_genes_mask2 = np.array(pd.DataFrame(filter_result)['gene_subset'])
    adata_training = adata_training[:,filter_genes_mask2]

    ### Remove cell cycle and house keeping genes - filter genes mask3
    # cell cycle genes:
    cc_genes_dan = ['cdk1','mcm2','mcm7','rrm2','cenpa', 'cdc6', 'ccnf', 'cdca4','ccnd1', 'kif4']

    cc_iter1  = corr_coeff(adata_training , cc_genes_dan, 0.3)[1]
    print ('first iteration of cell cycle genes search done')
    cc_iter2  = corr_coeff(adata_training , cc_iter1, 0.3)[1]
    print ('second iteration of cell cycle genes search done')
    cc_to_remove  = list(set(list(cc_iter1) + list(cc_iter2)))
    house_keeping = ['hmgb1b', 'hmgb3a', 'hspd1', 'hspa9', 'rplp0', 'hnrnpaba', 'rps2', 'rps12', 'rpl12', 'rps13', 'rps14', 'rps15a', 'rpl18', 'rps3a', 'rpl31', 'rpl37', 'rps6', 'rpl9', 'rpl11', 'rpl34', 'rpl13', 'rpl36a', 'rpl26', 'rps8a', 'rpl21', 'rps27.1', 'rpl27a', 'cirbpb'] 

    hk_iter1  = corr_coeff(adata_training , house_keeping, 0.3)[1]
    print ('first iteration of house keeping genes done')
    hk_iter2  = corr_coeff(adata_training , hk_iter1 , 0.3)[1]
    print ('first iteration of house keeping genes done')
    hk_to_remove  = list(set(list(hk_iter1) + list(hk_iter2)))

    to_remove  = np.array(cc_to_remove + hk_to_remove )
    gene_removal  =np.in1d(adata_training.var_names,to_remove ) # remove ribo, mito or cell cycle genes

    filter_genes_mask3 = (1-gene_removal).astype(bool)
    adata_training  = adata_training[:,filter_genes_mask3] # update adata_training 


    ## Get back to your data to filtering and further processing

    #### Total count normalization:
    adata_test.obs['n_counts'].mean()

    sc.pp.normalize_per_cell(adata_test, counts_per_cell_after=adata_training.obs['n_counts'].mean())


    #### Store unfiltered data as `adata_test.raw` - used in later steps
    adata_test.raw = adata_test

    #### Filter data using ref data masks
    adata_test = adata_test[:,filter_genes_mask1]
    adata_test = adata_test[:,filter_genes_mask2]
    adata_test = adata_test[:,filter_genes_mask3]


    ### scale data, find PCs and tranform data in the PC space
    m_ref = adata_training.X.mean(0) # mean
    s_ref = np.sqrt(sparse_var(adata_training.X)) # standard deviation
    Z_ref = sparse_multiply((adata_training.X - m_ref).T, 1/s_ref).T
    Z_test = sparse_multiply((adata_test.X - m_ref).T, 1/s_ref).T

    from sklearn.decomposition import PCA

    pca = PCA(n_components = ncomponents)

    pca.fit(Z_ref)


    return  adata_training, adata_test, pca.transform(Z_ref), pca.transform(Z_test)


def classifier_pca_2022(adata_ref, ref_states, adata_query, classifier, counts_per_cell = 10000, 
                        ncomponents=50, n_iterations=1000, seed=0):
    
    # process reference data:
    sc.pp.filter_genes(adata_ref, min_cells=3)
    sc.pp.filter_genes(adata_ref, min_counts=6)
    adata_ref.layers['raw_counts'] = adata_ref.X.copy()
    sc.pp.normalize_per_cell(adata_ref, counts_per_cell_after=1e4)
    sc.pp.log1p(adata_ref)
    sc.pp.highly_variable_genes(adata_ref, flavor='cell_ranger', n_top_genes=4000)
    sc.pl.highly_variable_genes(adata_ref)

    # scale reference data:
    e_ref = adata_ref[:,np.array(adata_ref.var['highly_variable'])].X
    m_ref = e_ref.mean(0) # mean
    s_ref = np.sqrt(sparse_var(e_ref)) # standard deviation
    Z_ref = sparse_multiply( (e_ref - m_ref).T, 1/s_ref).T  # zscore reference 

    # process and scale query data:
    sc.pp.normalize_per_cell(adata_query, counts_per_cell_after=1e4)
    var_genes = adata_ref.var_names[adata_ref.var['highly_variable']]
    e_query = adata_query[:,var_genes].X
    Z_query = sparse_multiply( (e_query - m_ref).T, 1/s_ref).T  # zscore query 

    # find PCs and project data on pcs:
    from sklearn.decomposition import PCA
    pca = PCA(n_components = 50) # find PCs
    pca.fit(Z_ref) 

    pca_ref = pca.transform(Z_ref) # tranform reference data on the new coordinates
    pca_query = pca.transform(Z_query) # transform query data on the new coordinates

    state_query = classifier(pca_ref, ref_states, pca_query, n_iterations, seed)

    return state_query

def classifier_pca(adata_ref, ref_states, adata_query, classifier, unimportant_genes, counts_per_cell = 10000, 
    remove_unimportant_genes=True, ncomponents=50, n_iterations=1000, seed=0, **args):
    '''
    input arguments:
    ----------------
    adata_ref         : reference anndata
    ref_states        : states labels of the reference data
    adata_query       : the data that needs to be classified
    unimportant_genes : if you have a list of genes that you want to remove by finding highly correlated genes
    n_component       : number of principle components

    returns:
    -------
    state_query    : states of the query dataset
    
    '''
    ### process reference data and make masks to be used for query data
    
    # normalized data:
    sc.pp.normalize_per_cell(adata_ref, counts_per_cell_after=counts_per_cell)

    # save raw data:
    adata_ref.raw = adata_ref
    
    # remove lowly expressed genes - filter genes mask 1
    filter_low_express = sc.pp.filter_genes(adata_ref.X, min_cells=3)
    filter_genes_mask1 = filter_low_express[0]
    adata_ref = adata_ref.copy()[:,filter_genes_mask1]

    # find variable genes - filter genes mask 2
    filter_result = sc.pp.filter_genes_dispersion(
        adata_ref.X, min_mean=0.0125, max_mean=3, min_disp=0.5)
    print('{} genes passing filter'.format(filter_result['gene_subset'].sum()))
    #sc.pl.filter_genes_dispersion(filter_result)
    filter_genes_mask2 = np.array(pd.DataFrame(filter_result)['gene_subset'])
    adata_ref = adata_ref.copy()[:,filter_genes_mask2]
    
    # Remove cell cycle and house keeping genes - filter genes mask3 (optional)
    if remove_unimportant_genes == True:
        unimportant_genes_iter1  = corr_coeff(adata_ref , unimportant_genes, 0.3)[1]
        print ('first iteration of searching unimportant genes done')
        unimportant_genes_iter2  = corr_coeff(adata_ref , unimportant_genes_iter1, 0.3)[1]
        print ('second iteration of searching unimportant genes done')
        genes_to_remove  = list(set(list(unimportant_genes_iter1) + list(unimportant_genes_iter1)))
        
        unimportant_genes_mask  = np.in1d(adata_ref.var_names,genes_to_remove )
        filter_genes_mask3 = (1-unimportant_genes_mask).astype(bool)
        adata_ref  = adata_ref.copy()[:,filter_genes_mask3] # update adata_ref 

    ### Get back to your data and process it using reference data
    
    # Total count normalization:
    sc.pp.normalize_per_cell(adata_query, counts_per_cell_after=10000)

    # Filter data using ref data masks
    adata_query = adata_query.copy()[:,filter_genes_mask1]
    adata_query = adata_query.copy()[:,filter_genes_mask2]
    if remove_unimportant_genes == True:
        adata_query = adata_query.copy()[:,filter_genes_mask3]


    ### scale data, find PCs and tranform data in the PC space (essentially wolbatch)
    m_ref = adata_ref.X.mean(0) # mean
    s_ref = np.sqrt(sparse_var(adata_ref.X)) # standard deviation
    Z_ref = sparse_multiply((adata_ref.X - m_ref).T, 1/s_ref).T  # zscore reference 
    Z_ref = np.array(Z_ref)
    Z_test = sparse_multiply((adata_query.X - m_ref).T, 1/s_ref).T # zscore test 
    Z_test = np.array(Z_test)


    from sklearn.decomposition import PCA
    pca = PCA(n_components = ncomponents) # find PCs
    pca.fit(Z_ref) 

    pca_ref = pca.transform(Z_ref) # tranform reference data on the new coordinates
    pca_query = pca.transform(Z_test) # transform query data on the new coordinates

    state_query = classifier(pca_ref, ref_states, pca_query, n_iterations, seed)
    
    return state_query



def process_data(adata, unimportant_genes, remove_unimportant_genes = True, counts_per_cell=10000, neigbhor_number =5,
                 batch_correction = True, batch_key='replicate',  ncomponents = 30, seed = 0, **args):

    '''
    description: normalizes data, filters genes expressed in less than 3 cells, 
    removes unimportant genes, find k nearest neighbors, PCs, louvain clusters, umap coordinates

    input arguments:
    ----------------
    adata             : reference anndata
    unimportant_genes : if you have a list of genes that you want to remove by finding highly correlated genes
    n_component       : number of principle components

    returns:
    -------
    adata_raw    : raw adata
    adata        : processed adata
    '''

    # normalize data:
    sc.pp.normalize_per_cell(adata, counts_per_cell_after=10000)
    adata_raw = adata.copy()
    adata.raw = adata
    adata.layers['raw_counts'] = adata.X.copy()

    sc.pp.filter_genes(adata, min_cells=3)

    # highly variable genes:
    filter_result = sc.pp.filter_genes_dispersion(adata.X, min_mean=0.01, max_mean=3, min_disp=0.5)
    print('{} genes passing filter'.format(filter_result['gene_subset'].sum()))
    #sc.pl.filter_genes_dispersion(filter_result)

    adata = adata[:,np.array(pd.DataFrame(filter_result)['gene_subset'])] # adata with only variable genes


    # find unimportant genes:
    if remove_unimportant_genes == True:

        unimportant_genes_iter1  = corr_coeff(adata, unimportant_genes, 0.3)[1]
        print ('first iteration of finding unimportant genes done')
        unimportant_genes_iter2 = corr_coeff(adata, unimportant_genes_iter1, 0.3)[1]
        print ('second iteration of finding unimportant genes done')
        genes_to_remove  = list(set(unimportant_genes + list(unimportant_genes_iter1) + list(unimportant_genes_iter2)))
        unimportant_genes_mask  = np.in1d(adata.var_names,genes_to_remove )
        adata  = adata[:,(1-unimportant_genes_mask).astype(bool)] # update adata 

    # scale data and find PCs
    sc.pp.scale(adata)
    sc.tl.pca(adata, n_comps= ncomponents, random_state = seed)

    # knn, louvain, umap
    if batch_correction == True:
        import bbknn
        bbknn.bbknn(adata,batch_key)
    
    else:
        sc.pp.neighbors(adata, n_neighbors=neigbhor_number, use_rep='X_pca')
    
    sc.tl.umap(adata, random_state = seed)
    sc.tl.leiden(adata)
    #sc.tl.louvain(adata)

    return adata


def start_subplot_figure(n_subplots, n_columns=5, fig_width=14, row_height=3, dpi=75):
    n_rows = int(np.ceil(n_subplots / float(n_columns)))
    fig = plt.figure(figsize = (fig_width, n_rows * row_height), dpi=dpi)
    return fig, n_rows, n_columns








