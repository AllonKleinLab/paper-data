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
import scanpy as sc
import pandas as pd
from distutils.version import LooseVersion, StrictVersion


########## LOADING DATA
def load_genes(filename, delimiter='\t', column=0, skip_rows=0):
    ''' Load gene list from a file

    Arguments
    - filename : str
        Name of file containing gene names
    - delimiter : str, optional (default: "\t")
        Column delimiter
    - column : int, optional (default: 0)
        Column containing gene names
    - skip_rows : int, optional (default: 0)
        Number of rows to skip at beginning of file

    Returns
    - gene_list : list, length n_genes
        List of gene names

    '''
    gene_list = []

    with open(filename) as f:
        for iL,l in enumerate(f):
            if iL >= skip_rows:
                gene_list.append(l.strip('\n').split(delimiter)[column])
    return gene_list

def make_genes_unique(orig_gene_list):
    ''' Make gene names unique by adding "__1", "__2", etc. to end of duplicate gene names
    
    Arguments
    - orig_gene_list : list, length n_genes
        List of gene names possibly containing duplicates

    Returns
    - gene_list : list, length n_genes
        List of unique gene names
    '''
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

def load_pickle(filename):
    '''Load data from pickle file
    Attempts to load pickle file using pickle library, falling back on pandas read_pickle if there 
    are version issues.

    Arguments
    - filename : str
        Pickle file name

    Returns
    - dat : object
        Object loaded from filename
    '''
    try:
        import pickle
        dat = pickle.load(open(filename, 'rb'))
    except:
        import pandas as pd
        dat = pd.read_pickle(filename)
    return dat

### loading counts

def file_opener(filename):
    '''Open file and return a file object, automatically decompressing zip and gzip 

    Arguments
    - filename : str
        Name of input file

    Returns
    - outData : file object
        (Decompressed) file data

    '''
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
        outData = open(filename, 'r')
    return outData

def load_mtx(file_data):
    ''' Load mtx file as scipy.sparse.csc_matrix

    Arguments
    - file_data : str or file object
        Name of input file or a file object

    Returns
    - scipy.sparse.csc_matrix

    '''
    return scipy.io.mmread(file_data).tocsc()

def load_npz(file_data):
    ''' Load scipy.sparse npz file as scipy.sparse.csc_matrix

    Arguments
    - file_data : str or file object
        Name of input file or a file object

    Returns
    - scipy.sparse.csc_matrix
    '''
    return scipy.sparse.load_npz(file_data).tocsc()

def load_npy(file_data):
    ''' Load npy file, converting to scipy.sparse.csc_matrix

    Arguments
    - file_data : str or file object
        Name of input file or a file object

    Returns
    - scipy.sparse.csc_matrix
    '''
    return scipy.sparse.csc_matrix(np.load(file_data))

def load_text(file_data, delim='\t', start_column=None, start_row=None, print_row_interval=None):
    '''Load text file as scipy.sparse.csc_matrix
    If start_column is not specificied, attempts to automatically identify the counts matrix
    (all numeric columns). 

    '''
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


def load_text2(file_data, delim='\t', start_column=0, start_row=0, row_indices=None, column_indices=None, print_row_interval=None):
    '''Load text file as scipy.sparse.csc_matrix
    Can load a user-specified subset of rows and/or columns
    '''

    X_data = []
    X_row = []
    X_col = []
    
    ncol = None
    nrow = 0

    for row_ix, dat in enumerate(file_data):
        if print_row_interval is not None:
            if (row_ix+1) % print_row_interval == 0:
                print('Row {}'.format(row_ix+1))
        if type(dat) == bytes:
            dat = dat.decode('utf-8')
        dat = dat.strip('\n').split(delim)

        if row_ix >= start_row:
            read_row = True
            if row_indices is not None:
                if (row_ix - start_row) not in row_indices:
                    read_row = False

            if read_row:
                rowdat = dat[start_column:]
                if ncol is None:
                    ncol = len(rowdat)
                else:
                    if len(rowdat) != ncol:
                        return 'ERROR: Rows have different numbers of numeric columns.'

                if column_indices is None:
                    column_indices = np.arange(ncol)

                rowdat = np.array(list(map(float, rowdat)))[column_indices]
                col_ix = np.nonzero(rowdat)[0]
                X_col.extend(col_ix)
                X_row.extend([nrow] * len(col_ix))
                X_data.extend(rowdat[col_ix])

                nrow += 1

    ncol = len(column_indices)
    print(nrow,ncol)
    E = scipy.sparse.coo_matrix((X_data, (X_row, X_col)), dtype=float, shape=(nrow, ncol)).tocsc()
    
    return E


def load_annotated_text(file_data, delim='\t', read_row_labels=False, read_column_labels=False, transpose=False, chunk_size=2000):
    '''Load text file as scipy.sparse.csc_matrix, returning column and/or row labels if desired.
    Loads rows in chunks to ease memory demands.
    '''

    X_data = []
    X_row = []
    X_col = []
    
    ncol = None
    nrow = 0

    row_labels = []
    column_labels = []

    E_chunks = []

    for row_ix, dat in enumerate(file_data):
        if type(dat) == bytes:
            dat = dat.decode('utf-8')
            
        dat = dat.strip('\n').split(delim)

        if read_column_labels and row_ix == 0:
            if read_column_labels:
                column_labels = dat[1:]
            else:
                column_labels = dat
        else:
            if read_row_labels:
                row_labels.append(dat[0])
                rowdat = dat[1:]
            else:
                rowdat = dat[0:]
            if ncol is None:
                ncol = len(rowdat)
            else:
                if len(rowdat) != ncol:
                    return 'ERROR: Line {} has {} columns. Previous line(s) had {}'.format(row_ix, len(rowdat), ncol)
            rowdat = np.array(list(map(float, rowdat)))
            col_ix = np.nonzero(rowdat)[0]
            X_col.extend(col_ix)
            X_row.extend([nrow] * len(col_ix))
            X_data.extend(rowdat[col_ix])
            nrow += 1

            if chunk_size is not None:
                if nrow % chunk_size == 0:
                    E_chunks.append(scipy.sparse.coo_matrix((X_data, (X_row, X_col)), dtype=float, shape=(nrow, ncol)))
                    X_data = []
                    X_row = []
                    X_col = []
                    nrow = 0

    if nrow > 0:
        E_chunks.append(scipy.sparse.coo_matrix((X_data, (X_row, X_col)), dtype=float, shape=(nrow, ncol)))

    E = scipy.sparse.vstack(E_chunks)
    if transpose: 
        E = E.T
    
    return E.tocsc(), np.array(row_labels), np.array(column_labels)

def load_cellranger_h5_v2(filename, genome):
    import h5py
    import scipy.sparse as ssp
    
    f = h5py.File(filename, 'r')
    barcodes = np.array(f.get(genome).get('barcodes')).astype(str)
    gene_names = np.array(f.get(genome).get('gene_names')).astype(str)

    
    data = np.array(f.get(genome).get('data'))
    indices = np.array(f.get(genome).get('indices'))
    indptr = np.array(f.get(genome).get('indptr'))
    shape = np.array(f.get(genome).get('shape'))

    # Make sparse expression matrix
    E = ssp.csc_matrix((data, indices, indptr), shape=shape).T.tocsc()

    f.close()

    return E, barcodes, gene_names

def load_cellranger_h5_v3(filename):
    import h5py
    import scipy.sparse as ssp
    
    f = h5py.File(filename, 'r')
    barcodes = np.array(f.get('matrix').get('barcodes')).astype(str)
    gene_names = np.array(f.get('matrix').get('features').get('name')).astype(str)

    
    data = np.array(f.get('matrix').get('data'))
    indices = np.array(f.get('matrix').get('indices'))
    indptr = np.array(f.get('matrix').get('indptr'))
    shape = np.array(f.get('matrix').get('shape'))

    # Make sparse expression matrix
    E = ssp.csc_matrix((data, indices, indptr), shape=shape).T.tocsc()

    f.close()

    return E, barcodes, gene_names


########## USEFUL SPARSE FUNCTIONS

def sparse_var(E, axis=0):
    ''' calculate variance across the specified axis of a sparse matrix'''

    mean_gene = E.mean(axis=axis).A.squeeze()
    tmp = E.copy()
    tmp.data **= 2
    return tmp.mean(axis=axis).A.squeeze() - mean_gene ** 2

def sparse_rowwise_multiply(E, a):
    ''' multiply each row of sparse matrix by a scalar '''

    nrow = E.shape[0]
    w = scipy.sparse.lil_matrix((nrow, nrow))
    w.setdiag(a)
    return w * E

def mean_center(E, column_means=None):
    ''' mean-center columns of a sparse matrix '''

    if column_means is None:
        column_means = E.mean(axis=0)
    return E - column_means

def normalize_variance(E, column_stdevs=None):
    ''' variance-normalize columns of a sparse matrix '''

    if column_stdevs is None:
        column_stdevs = np.sqrt(sparse_var(E, axis=0))
    return sparse_rowwise_multiply(E.T, 1 / column_stdevs).T

def sparse_zscore(E, gene_mean=None, gene_stdev=None):
    ''' z-score normalize each column of a sparse matrix '''
    if gene_mean is None:
        gene_mean = E.mean(0)
    if gene_stdev is None:
        gene_stdev = np.sqrt(sparse_var(E))
    return sparse_rowwise_multiply((E - gene_mean).T, 1/gene_stdev).T


########## CELL FILTERING
def filter_dict(d, filt):
    ''' filter 1-D and 2-D entries in a dictionary '''
    for k,v in d.items():
        if k != 'meta':
            if len(v.shape) == 1:
                d[k] = v[filt]
            else:
                d[k] = v[filt,:]
    return d

########## GENE FILTERING

def runningquantile(x, y, p, nBins):
    ''' calculate the quantile of y in bins of x '''

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
    Calculate v-score (above-Poisson noise statistic) for genes in the input sparse counts matrix
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
        plt.figure(figsize=(4, 3));
        plt.scatter(np.log10(mu_gene), np.log10(FF_gene), c = [[.8,.8,.8]], alpha = 0.3, s = 3);
        plt.scatter(np.log10(mu_gene)[ix], np.log10(FF_gene)[ix], c = [[0,0,0]], alpha = 0.3,  s= 3);
        plt.plot(np.log10(xTh),np.log10(yTh));
        plt.title(sample_name)
        plt.xlabel('log10(mean)');
        plt.ylabel('log10(Fano factor)');
        plt.show()

    return gene_ix[ix]

def remove_corr_genes(E, gene_list, exclude_corr_genes_list, test_gene_idx, min_corr = 0.1):
    ''' remove signature-correlated genes from a list of test genes 
    Arguments:
    E: scipy.sparse.csc_matrix, shape (n_cells, n_genes)
        - full counts matrix
    gene_list: numpy array, shape (n_genes,)
        - full gene list
    exclude_corr_genes_list: list of list(s)
        - Each sublist is used to build a signature. Test genes correlated
          with this signature will be removed
    test_gene_idx: 1-D numpy array
        - indices of genes to test for correlation with the 
          gene signatures from exclude_corr_genes_list
    min_corr: float (default=0.1)
        - Test genes with a Pearson correlation of min_corr or higher 
          with any of the gene sets from exclude_corr_genes_list will
          be excluded

    Returns:
        numpy array of gene indices (subset of test_gene_idx) that 
        are not correlated with any of the gene signatures
    '''
    seed_ix_list = []
    for l in exclude_corr_genes_list:
        seed_ix_list.append(np.array([i for i in range(len(gene_list)) if gene_list[i] in l], dtype=int))

    exclude_ix = []
    for iSet in range(len(seed_ix_list)):
        seed_ix = seed_ix_list[iSet][E[:,seed_ix_list[iSet]].sum(axis=0).A.squeeze() > 0]
        if type(seed_ix) is int:
            seed_ix = np.array([seed_ix], dtype=int)
        elif type(seed_ix[0]) is not int:
            seed_ix = seed_ix[0]
        indat = E[:, seed_ix]
        tmp = sparse_zscore(indat)
        tmp = tmp.sum(1).A.squeeze()

        c = np.zeros(len(test_gene_idx))
        for iG in range(len(c)):
            c[iG],_ = scipy.stats.pearsonr(tmp, E[:,test_gene_idx[iG]].A.squeeze())

        exclude_ix.extend([test_gene_idx[i] for i in range(len(test_gene_idx)) if (c[i]) >= min_corr])
    exclude_ix = np.array(exclude_ix)

    return np.array([g for g in test_gene_idx if g not in exclude_ix], dtype=int)

########## SAMPLE PROJECTION
def classifier_find_shared_genes(adata_ref, adata_query):
    ''' Find genes present in both AnnData objects '''
    genes1 = adata_ref.var_names.values.astype(str)
    genes2 = adata_query.var_names.values.astype(str)
    shared_genes = np.intersect1d(genes1, genes2)

    adata_ref.var['valid_for_classification'] = False
    adata_ref.var.loc[shared_genes, 'valid_for_classification'] = True
    adata_query.var['valid_for_classification'] = False
    adata_query.var.loc[shared_genes, 'valid_for_classification'] = True
    
    return


def classifier_learn_pca_projection(adata_ref, highvar_genes_column='highly_variable_projection', cell_normalize=True, find_highvar_genes=True, mean_center=True, n_pcs=30, n_neighbors=10, knn_metric='euclidean', random_state=1, **kwargs):
    ''' Preprocess reference data to "learn" PCA/TruncatedSVD loadings and initialize KNN
    Requires running `projection_find_shared_genes()` first. 

    Total counts normalization, gene filtering, and decomposition. Automatically stores
    intermediates needed to apply same transformations to a "query" dataset.
    - kwargs: for filter_genes function (min_counts, min_cells, min_vscore_pctl)

    Most important results stored in:
    - adata_ref.uns['projector']
    - adata_ref.var
    '''
    valid_gene_names = adata_ref.var_names[adata_ref.var['valid_for_classification']]
    adata_tmp = adata_ref.copy()
    adata_tmp = adata_tmp[:, valid_gene_names].copy()

    if find_highvar_genes is False:
        if highvar_genes_column in adata_tmp.var.columns:
            highvar_gene_names = adata_tmp.var_names[adata_tmp.var[highvar_genes_column]]
        else:
            print(f'ERROR: could not find "{highvar_genes_column}" in `adata_ref.var` columns, but `find_highvar_genes` is `False`.')
            print('Please specify a valid column name for `highvar_genes_column` or set `find_highvar_genes` to `True`.')
            return
    else:
        if cell_normalize:
            sc.pp.normalize_per_cell(adata_tmp, counts_per_cell_after=adata_tmp.X.sum(1).A.squeeze().mean())
        highvar_gene_names = valid_gene_names[filter_genes(adata_tmp.X, **kwargs)]
        adata_ref.var[highvar_genes_column] = False
        adata_ref.var.loc[highvar_gene_names, highvar_genes_column] = True

    if cell_normalize:
        sc.pp.normalize_per_cell(adata_tmp, counts_per_cell_after=1e4)

    X = adata_tmp[:, highvar_gene_names].X.copy()
    del adata_tmp
    gene_means = X.mean(0).A.squeeze()
    gene_stdevs = np.sqrt(sparse_var(X))
    
    adata_ref.var['mean_projection'] = np.nan
    adata_ref.var['stdev_projection'] = np.nan
    adata_ref.var.loc[highvar_gene_names, 'mean_projection'] = gene_means
    adata_ref.var.loc[highvar_gene_names, 'stdev_projection'] = gene_stdevs
    
    adata_ref.uns['projector'] = {'params': {'cell_normalize': cell_normalize, 
                                             'mean_center': mean_center, 
                                             'highvar_genes_column': highvar_genes_column,
                                             'n_neighbors': n_neighbors
                                            }}
    
    if mean_center:
        adata_ref.uns['projector']['decomp'] = PCA(n_components=n_pcs, random_state=random_state)
        Z = sparse_zscore(X, gene_mean=gene_means, gene_stdev=gene_stdevs)
    else:
        adata_ref.uns['projector']['decomp'] = TruncatedSVD(n_components=n_pcs, random_state=random_state)
        Z = normalize_variance(X, column_stdevs=gene_stdevs)

    adata_ref.uns['projector']['decomp'].fit(Z)
    adata_ref.obsm['projection'] = adata_ref.uns['projector']['decomp'].transform(Z)
    
    adata_ref.uns['projector']['knn'] = NearestNeighbors(n_neighbors=n_neighbors, metric=knn_metric).fit(adata_ref.obsm['projection'])

    return

def classifier_apply_pca_projection(adata_ref, adata_query):
    ''' Projection query into ref stored in `adata_ref` and find nearest neighbors 
    Requires running `projection_learn()` first.

    Results stored in: 
    - adata_query.obsm['projection']
    - adata_query.obsm['projection_neighbors']
    - adata_query.obsm['projection_distances']
    '''

    proj_par = adata_ref.uns['projector']['params']
    hvg_column = proj_par['highvar_genes_column']
    
    adata_tmp = adata_query.copy()
    adata_tmp = adata_tmp[:, adata_tmp.var['valid_for_classification']].copy()
    
    if proj_par['cell_normalize']:
        sc.pp.normalize_per_cell(adata_tmp, counts_per_cell_after=1e4)
    
    highvar_gene_names = adata_ref.var_names[adata_ref.var[hvg_column]].values.astype(str)
    gene_means = adata_ref.var.loc[highvar_gene_names, 'mean_projection'].values
    gene_stdevs = adata_ref.var.loc[highvar_gene_names, 'stdev_projection'].values
    
    adata_tmp = adata_tmp[:, highvar_gene_names]
    X = adata_tmp.X.copy()
    del adata_tmp
    
    if proj_par['mean_center']:
        Z = sparse_zscore(X, gene_mean=gene_means, gene_stdev=gene_stdevs)
    else:
        Z = normalize_variance(X, column_stdevs=gene_stdevs)
    
    adata_query.obsm['projection'] = adata_ref.uns['projector']['decomp'].transform(Z)
    
    knn_dist, knn_neigh = adata_ref.uns['projector']['knn'].kneighbors(
        adata_query.obsm['projection'], 
        return_distance=True
    )
    adata_query.obsm['projection_neighbors'] = knn_neigh
    adata_query.obsm['projection_distances'] = knn_dist
    
    return


def classifier_learn_max_likelihood(adata_ref, column_name, use_raw=None):
    if use_raw is None:
        use_raw = True if adata_ref.raw is not None else False

    adata_ref.uns['max_likelihood'] = {
        'params': {
            'use_raw': use_raw,
            'column_name': column_name
        }
    }

    if use_raw:
        adata_use = sc.AnnData(adata_ref.raw.X.copy())
        adata_use.var_names = adata_ref.raw.var_names.values.astype(str)
    else:
        adata_use = sc.AnnData(adata_ref.X.copy())
        adata_use.var_names = adata_ref.var_names.values.astype(str)

    groups = adata_ref.obs[column_name].values.astype(str)
    group_centroids = scipy.sparse.lil_matrix((len(set(groups)), adata_use.shape[1]))
    group_numbers = np.zeros(len(groups), dtype=int)
    group_labels = []
    for iG,g in enumerate(set(groups)):
        ix = groups == g
        group_centroids[iG,:] = adata_use.X[ix,:].sum(0)
        group_numbers[groups == g] = iG
        group_labels.append(g)
    group_centroids = group_centroids.toarray()
    adata_ref.uns['max_likelihood']['profiles'] = pd.DataFrame(
        data=group_centroids, 
        index=group_labels, 
        columns=adata_use.var_names.values.astype(str))
    return

def classifier_apply_max_likelihood(adata_ref, adata_query, query_use_raw=None, pseudocount=1, postnorm_total=1e4, min_shared_genes=100):
    if query_use_raw is None:
        query_use_raw = True if adata_query.raw is not None else False

    if query_use_raw:
        adata_query_use = sc.AnnData(adata_query.raw.X.copy())
        adata_query_use.var_names = adata_query.raw.var_names.values.astype(str)
    else:
        adata_query_use = sc.AnnData(adata_query.X.copy())
        adata_query_use.var_names = adata_query.var_names.values.astype(str)

    genes1 = adata_ref.uns['max_likelihood']['profiles'].columns.values.astype(str)
    genes2 = adata_query_use.var_names.values.astype(str)
    shared_genes = np.intersect1d(genes1, genes2)
    
    if len(shared_genes) >= min_shared_genes:
        adata_query_use = adata_query_use[:, shared_genes]

        profiles = adata_ref.uns['max_likelihood']['profiles'].loc[:, shared_genes].values.copy()
        profile_labels = adata_ref.uns['max_likelihood']['profiles'].index.values.astype(str)

        if scipy.sparse.issparse(profiles):
            profiles = profiles.toarray()
        
        profiles = profiles * postnorm_total / np.sum(profiles, axis=1)[:,None]
        profiles += pseudocount
        profiles = profiles / np.sum(profiles, axis=1)[:,None]
        profiles = np.log10(profiles)
        
        scores = adata_query_use.X.dot(profiles.T)
        adata_query.obs['max_likelihood_label'] = profile_labels[scores.argmax(1)]

    return


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


def get_pca(E, base_ix=[], numpc=50, keep_sparse=False, normalize=True, random_state=0):
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
            Z = sparse_rowwise_multiply(E.T, 1 / zstd).T
        else:
            Z = E
        pca = TruncatedSVD(n_components=numpc, random_state=random_state)

    else:
        if normalize:
            zmean = E[base_ix,:].mean(0)
            zstd = np.sqrt(sparse_var(E[base_ix,:]))
            Z = sparse_rowwise_multiply((E - zmean).T, 1/zstd).T
        else:
            Z = E
        pca = PCA(n_components=numpc, random_state=random_state)

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

def get_knn_graph(X, k=5, dist_metric='euclidean', approx=False, return_edges=True, random_seed=0):
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
        annoy_index.set_seed(random_seed)

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

def get_smooth_values(input_values, adjacency_matrix, beta=0.1, n_rounds=10):
    ''' Smooth values on a graph
    
    Example: 
    (Starting with adata containing PCA coordinates in `adata.obsm['X_pca']`)

    # Get knn edges
    edges = hf.get_knn_graph(adata.obsm['X_pca'], k=10, return_edges=True)[0]

    # Build adjacency matrix from edges
    adjacency_matrix = hf.build_adj_mat(edges, adata.shape[0])

    # Get raw data
    input_data = adata[:, 'Klf1'].X

    # Get smoothed values
    smoothed_value = get_smooth_values(input_data, adjacency_matrix, beta=0.1, n_rounds=10)
    '''

    if len(input_values.shape) == 1:
        smoothed_values = input_values[:, None]
    else:
        smoothed_values = input_values.copy()

    adjacency_matrix = sparse_rowwise_multiply(adjacency_matrix, 1 / adjacency_matrix.sum(1).A.squeeze())
    for iRound in range(n_rounds):
        smoothed_values = (beta * smoothed_values + ((1 - beta) * adjacency_matrix) * smoothed_values)

    return smoothed_values.squeeze()

########## CLUSTERING

def get_spectral_clusters(A, k):
    from sklearn.cluster import SpectralClustering
    spec = SpectralClustering(n_clusters=k, random_state = 0, affinity = 'precomputed', assign_labels = 'discretize')
    return spec.fit_predict(A)

def get_louvain_clusters(n_nodes, edges, resolution=1.0, random_seed=0):
    #https://github.com/theislab/scanpy/blob/master/scanpy/tools/louvain.py
    import louvain
    import igraph as ig
    g = ig.Graph()
    g.add_vertices(n_nodes)
    g.add_edges(list(edges))

    partition_type = louvain.RBConfigurationVertexPartition
    louvain.set_rng_seed(random_seed)
    part = louvain.find_partition(
        g, partition_type,
        resolution_parameter=resolution
    )

    return np.array(part.membership, dtype=int)

def get_louvain_clusters_old(nodes, edges):
    import networkx as nx
    import community
    
    G = nx.Graph()
    G.add_nodes_from(nodes)
    G.add_edges_from(edges)
    
    return np.array(list(community.best_partition(G).values()))

########## 2-D EMBEDDING
def get_tsne(X, angle=0.5, perplexity=30, verbose=False):
    from sklearn.manifold import TSNE
    return TSNE(angle=angle, perplexity=perplexity, verbose=verbose).fit_transform(X)

def get_force_layout(X, n_neighbors=5, approx_neighbors=False, n_iter=300, verbose=False):
    edges = get_knn_graph(X, k=n_neighbors, approx=approx_neighbors, return_edges=True)[0]
    return run_force_layout(edges, X.shape[0], n_iter=n_iter, verbose=verbose)

def run_force_layout(links, n_cells, n_iter=100, edgeWeightInfluence=1, barnesHutTheta=2, scalingRatio=1, gravity=0.05, jitterTolerance=1, verbose=False):
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

def get_umap(X, n_neighbors=10, min_dist=0.1, spread=1.0, metric='euclidean', random_state=1):
    import umap
    # UMAP: https://github.com/lmcinnes/umap
    embedding = umap.UMAP(
        n_neighbors=n_neighbors, 
        min_dist=min_dist, 
        spread=spread, 
        metric=metric, 
        random_state=random_state
        ).fit_transform(X)
    return embedding

########## SCANPY PREPROCESSING

def process_ad(adata, normalize=True, min_counts=3, min_cells=3, min_vscore_pctl=85, plot_vscore=False, 
               n_components=30, n_neighbors=15, umap_min_dist=0.3, leiden_res=1.0, sparse_pca=False,
               log_transform=False, scale=True, copy=False, random_state=0,
               exclude_gene_sets=None, exclude_gene_corr=0.2, exclude_gene_list=None, verbose=False, batch_base_mask=None, run_umap=True, run_leiden=True
              ):
    
    if copy:
        adata = adata.copy()
    
    if batch_base_mask is None:
        batch_base_mask = np.repeat(True, adata.shape[0])
    if normalize:
        if verbose:
            print('    normalizing...')
        adata.raw = adata
        adata.obs['n_counts'] = adata.X.sum(1).A.squeeze()
        if LooseVersion(sc.__version__) > LooseVersion('1.4.1'):
            sc.pp.normalize_total(adata, target_sum=adata.obs['n_counts'].values[batch_base_mask].mean())
        else:
            sc.pp.normalize_per_cell(adata, 
                               counts_per_cell_after=adata.obs['n_counts'].values[batch_base_mask].mean())
        
        
    if verbose:
        print('    filtering genes...')
    hvg = filter_genes(
        adata.X[batch_base_mask,:], 
        min_counts=min_counts, 
        min_cells=min_cells, 
        min_vscore_pctl=min_vscore_pctl, 
        show_vscore_plot=plot_vscore)
    adata.var['highly_variable'] = False
    adata.var.loc[adata.var_names[hvg], 'highly_variable'] = True
    if verbose:
        print(f'        {len(hvg)} highly variable genes')
    if exclude_gene_sets is not None:
        hvg = remove_corr_genes(
            adata.X[batch_base_mask,:], 
            adata.var_names.values.astype(str),
            exclude_gene_sets,
            hvg,
            min_corr=exclude_gene_corr
        )
        if verbose:
            print(f'        {len(hvg)} genes pass correlation filter')
    adata.var['highly_variable_prefilt'] = adata.var['highly_variable'].copy()
    adata.var['highly_variable'] = False
    adata.var.loc[adata.var_names[hvg], 'highly_variable'] = True
    if exclude_gene_list is not None:
        exclude_gene_list = [g for g in exclude_gene_list if g in adata.var_names]
        adata.var.loc[exclude_gene_list, 'highly_variable'] = False
    
    if normalize:
        if LooseVersion(sc.__version__) > LooseVersion('1.4.1'):
            sc.pp.normalize_total(adata, target_sum=1e4, key_added='n_counts2')
        else:
            sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4, key_n_counts='n_counts2')
    
    pca_input = adata.X[:, adata.var['highly_variable'].values].copy()
    if log_transform:
        pca_input.data = np.log10(1 + pca_input.data)
    if verbose:
        print('    running PCA...')
    adata.obsm['X_pca'] = get_pca(pca_input, 
                                     numpc=n_components, 
                                     keep_sparse=sparse_pca, 
                                     normalize=scale,
                                     random_state=random_state,
                                     base_ix=batch_base_mask
                                    )
    del pca_input
    if verbose:
        print('    finding neighbors...')
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, random_state=random_state)
    if run_umap:
        if verbose:
            print('    running UMAP...')
        sc.tl.umap(adata, min_dist=umap_min_dist, random_state=random_state)
    if run_leiden:
        if verbose:
            print('    running Leiden clustering...')
        sc.tl.leiden(adata, resolution=leiden_res)
    
    return adata if copy else None


########## GENE ENRICHMENT

def rank_enriched_genes(E, gene_list, cell_mask, min_counts=3, min_cells=3, verbose=False):
    gix = (E[cell_mask,:]>=min_counts).sum(0).A.squeeze() >= min_cells
    if verbose:
        print('%i cells in group' %(sum(cell_mask)))
        print('Considering %i genes' %(sum(gix)))
    
    gene_list = gene_list[gix]
    E = E.copy()[:, gix]

    means = E.mean(0).A.squeeze()
    stdevs = np.sqrt(sparse_var(E))
    scores = (E[cell_mask,:].mean(0).A.squeeze() - means) / stdevs
    
    o = np.argsort(-scores)
    
    return gene_list[o], scores[o]

def get_dge(ad, mask1, mask2, min_frac_expr=0.05, pseudocount=1):
    import statsmodels.sandbox.stats.multicomp
    import scipy.stats
    
    gene_mask = ((ad.X[mask1,:]>0).sum(0).A.squeeze()/mask1.sum() > min_frac_expr) | ((ad.X[mask2,:]>0).sum(0).A.squeeze()/mask2.sum() > min_frac_expr)
    print(gene_mask.sum())
    E1 = ad.X[mask1,:][:,gene_mask].toarray()
    E2 = ad.X[mask2,:][:,gene_mask].toarray()
    
    m1 = E1.mean(0) + pseudocount
    m2 = E2.mean(0) + pseudocount
    r = np.log2(m1 / m2)
    
    pv = np.zeros(gene_mask.sum())
    for ii,iG in enumerate(np.nonzero(gene_mask)[0]):
        pv[ii] = scipy.stats.ranksums(E1[:,ii], E2[:,ii])[1]
    pv = statsmodels.sandbox.stats.multicomp.multipletests(pv, alpha=0.05, method='fdr_bh',)[1]
    
    df = pd.DataFrame({
        'gene': ad.var_names.values.astype(str)[gene_mask],
        'pv': pv,
        'm1': m1 - pseudocount, 
        'm2': m2 - pseudocount, 
        'ratio': r
    })
    
    return df

def find_markers(norm_counts, gene_list, groups, groups_find=None, min_frac_expr=0.01, min_fold_change=2, pseudocount=0.1, max_p=0.05, verbose=False):
    '''
    norm_counts: normalized counts matrix (scipy.sparse with shape (n_cells, n_genes))
    gene_list: numpy array of gene names (length = n_genes)
    min_frac_expr: only test genes expressed in this fraction of cells in the cluster of interest
    min_fold_change: min fold-change between highest avg expression and 2nd-highest avg expression
    pseudocount: pseudocount to use when calculating fold-change
    max_p: maximum multiple hypothesis-corrected p-value
    '''
    import pandas as pd
    import numpy as np
    import statsmodels.sandbox.stats.multicomp
    import scipy.stats
    
    ##########
    min_fold_change = np.log2(min_fold_change)

    # cluster labels for each cell
    clusts_use = groups.copy() 
    if groups_find is None:
        groups_find = np.unique(clusts_use)

    # initialize results
    results = {'group': [], 'gene': [], 'p-value': [], 'log2_fold_change': []}

    for c1 in groups_find:

        #############    SETUP 

        # Find cells in test cluster
        i1 = clusts_use == c1
        i2 = ~i1

        # Find genes expressed in some percent of cells 
        # in test cluster
        gix = ((norm_counts[i1,:]>0).sum(0).A.squeeze()/i1.sum() > min_frac_expr)
        if verbose:
            print('Testing {} genes for cluster {}'.format(gix.sum(), c1))


        #############    P-VALS

        E1 = norm_counts[i1,:][:,gix].toarray()
        E2 = norm_counts[i2,:][:,gix].toarray()
        pv = np.zeros(gix.sum())
        for ii,iG in enumerate(np.nonzero(gix)[0]):
            pv[ii] = scipy.stats.ranksums(E1[:,ii], E2[:,ii])[1]
        pv = statsmodels.sandbox.stats.multicomp.multipletests(pv, alpha=0.05, method='fdr_bh',)[1]

        del E1,E2


        #############    HIGHEST VS. 2ND-HIGHEST MEANS

        m1 = norm_counts[i1,:][:,gix].mean(0).A.squeeze()
        m2 = np.zeros(gix.sum())
        for c in np.unique(clusts_use):
            if c != c1:
                i_test = clusts_use == c
                m_test = norm_counts[i_test,:][:,gix].mean(0).A.squeeze()
                gix_change = m_test-m2>0
                m2[gix_change] = m_test[gix_change]

        logFc = np.log2((m1+pseudocount)/(m2+pseudocount))

        #############    IDENTIFY SIG
        hi = (pv < max_p) & (logFc > min_fold_change)
        if hi.sum() > 0:
            results['group'].extend([c1 for i in range(hi.sum())])
            results['gene'].extend(list(gene_list[gix][hi]))
            results['p-value'].extend(list(pv[hi]))
            results['log2_fold_change'].extend(list(logFc[hi]))

        if verbose:
            print('cluster={}, n_cells={}, n_diff={}'.format(c1, i1.sum(), hi.sum()))

    results = pd.DataFrame(results)[['group', 'gene', 'p-value', 'log2_fold_change']]
    return results

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


def write_spring_main(directory, counts_norm, gene_list, total_counts):
    if not os.path.exists(directory):
        os.makedirs(directory)

    save_hdf5_cells(counts_norm, '{}/counts_norm_sparse_cells.hdf5'.format(directory))
    save_hdf5_genes(counts_norm, gene_list, '{}/counts_norm_sparse_genes.hdf5'.format(directory))
    save_sparse_npz(counts_norm, '{}/counts_norm.npz'.format(directory), compressed=True)
    np.savetxt('{}/genes.txt'.format(directory), gene_list, fmt='%s')
    np.savetxt('{}/total_counts.txt'.format(directory), total_counts, fmt='%.1f')

    return

def rescale_coordinates(coords):
    coords = coords - np.min(coords, axis=0) - np.ptp(coords, axis=0) / 2.0
    coords = coords / np.ptp(coords, axis=0) * 30 * np.sqrt(coords.shape[0])
    coords[:,0] = coords[:,0] + 750
    coords[:,1] = coords[:,1] + 250
    return coords

def write_spring_subplot(directory, counts_norm, gene_list, cell_subset, coordinates, edges, 
    total_counts, prin_comps, gene_filter, obs_df=None, categ_vars=None, contin_vars=None):

    if not os.path.exists(directory):
        os.makedirs(directory)

    # save custom colors
    custom_colors = {'Total counts': total_counts}
    custom_colors['Uniform'] = np.zeros(counts_norm.shape[0])
    if contin_vars is not None:
        for k in contin_vars:
            custom_colors[k] = obs_df[k].values
    write_color_tracks(custom_colors, '{}/color_data_gene_sets.csv'.format(directory))

    # create and save a dictionary of color profiles to be used by the visualizer
    color_stats = {}
    color_stats = get_color_stats_genes(color_stats, counts_norm, gene_list)
    color_stats = get_color_stats_custom(color_stats, custom_colors)
    save_color_stats('{}/color_stats.json'.format(directory), color_stats)

    # save cell labels
    cell_groupings = {}
    categorical_coloring_data = {}
    if categ_vars is not None:
        for k in categ_vars:
            cell_groupings[k] = list(obs_df[k].values.astype(str))
    categorical_coloring_data = build_categ_colors(categorical_coloring_data, cell_groupings)
    save_cell_groupings('{}/categorical_coloring_data.json'.format(directory), categorical_coloring_data)

    # write graph
    write_graph('{}/graph_data.json'.format(directory), counts_norm.shape[0], edges)
    write_edges('{}/edges.csv'.format(directory), edges)

    # save useful intermediates:
    np.savez_compressed('{}/intermediates.npz'.format(directory), Epca=prin_comps, gene_filter=gene_filter, total_counts=total_counts)

    np.save('{}/cell_filter.npy'.format(directory), cell_subset)
    np.savetxt('{}/cell_filter.txt'.format(directory), cell_subset, fmt='%i')

    coordinates = rescale_coordinates(coordinates)
    np.savetxt('{}/coordinates.txt'.format(directory),
                       np.hstack((np.arange(coordinates.shape[0])[:,None], coordinates)), fmt='%i,%.5f,%.5f')
    
    info_dict = {}
    info_dict['Date'] = '{}'.format(datetime.now())
    info_dict['Nodes'] = prin_comps.shape[0]
    info_dict['Filtered_Genes'] = len(gene_filter)
    info_dict['Gene_Var_Pctl'] = None
    info_dict['Min_Cells'] = None
    info_dict['Min_Counts'] = None
    info_dict['Num_Neighbors'] = None
    info_dict['Num_PCs'] = prin_comps.shape[1]
    info_dict['Num_Force_Iter'] = None
    with open('{}/run_info.json'.format(directory),'w') as f:
        f.write(json.dumps(info_dict, indent=4, sort_keys=True))

    return

   


#========================================================================================#

def make_spring_subplot(E, gene_list, save_path, base_ix=None, normalize=True, exclude_dominant_frac=1.0, min_counts=3, min_cells=5, min_vscore_pctl=85, show_vscore_plot=False, exclude_gene_names=None, num_pc=30, sparse_pca=False, pca_norm=True, k_neigh=4, cell_groupings={}, num_force_iter=100, output_spring=True, precomputed_pca=None, gene_filter=None, custom_colors={}, exclude_corr_genes_list=None, exclude_corr_genes_minCorr = 0.2, dist_metric = 'euclidean', use_approxnn=False, run_doub_detector = False, dd_k=50, dd_frac=5, dd_approx=True, tot_counts_final = None):
    
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
                gene_filter = filter_genes(E, base_ix, min_vscore_pctl=min_vscore_pctl, min_counts=min_counts, min_cells=min_cells,show_vscore_plot=show_vscore_plot)
            else:
                gene_filter = np.arange(E.shape[1])

            if len(gene_filter) == 0:
                print('Error: No genes passed filter')
                sys.exit(2)

            if not exclude_corr_genes_list is None:
                gene_filter = remove_corr_genes(E, gene_list, exclude_corr_genes_list, gene_filter, min_corr=exclude_corr_genes_minCorr)
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
    else:
        print('Using user-supplied PCA coordinates')
        Epca = precomputed_pca
        out['Epca'] = Epca
        if gene_filter is not None:
            out['gene_filter'] = gene_filter

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
        positions = run_force_layout(links, Epca.shape[0], n_iter=num_force_iter, 
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

def set_plot_defaults(fontsize=10):
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = 'Arial'
    plt.rc('font', size=fontsize)
    plt.rcParams['pdf.fonttype'] = 42
    return

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

def plot_groups(x, y, groups, buffer_pct=0.03, saving=False, fig_dir='./', fig_name='fig', res=300, close_after=False, title_size=12, point_size=3, ncol=5, fig_width=14, row_height=3):
    import matplotlib.pyplot as plt

    n_col = int(ncol)
    ngroup = len(np.unique(groups))
    nrow = int(np.ceil(ngroup / float(ncol)))
    fig = plt.figure(figsize = (fig_width, row_height * nrow))
    for ii, c in enumerate(np.unique(groups)):
        ax = plt.subplot(nrow, ncol, ii+1)
        ix = groups == c

        ax.scatter(x[~ix], y[~ix], s=point_size, c=[[.8,.8,.8]], edgecolors = '')
        ax.scatter(x[ix], y[ix], s=point_size, c=[[0,0,0]], edgecolors = '')
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlim(x.min()-x.ptp()*buffer_pct, x.max()+x.ptp()*buffer_pct)
        ax.set_ylim(y.min()-y.ptp()*buffer_pct, y.max()+y.ptp()*buffer_pct)

        ax.set_title(str(c), fontsize = title_size)

    fig.tight_layout()

    if saving:
        if not os.path.exists(fig_dir):
            os.makedirs(fig_dir)
        plt.savefig(fig_dir + '/' + fig_name + '.png', dpi=res)

    if close_after:
        plt.close()

def plot_one_gene(E, gene_list, gene_to_plot, x, y, normalize=False, ax=None, order_points=True, col_range=(0,100), buffer_pct=0.03, point_size=1, color_map=None, smooth_operator=None):
    if color_map is None:
        color_map = darken_cmap(plt.cm.Reds,.9)
    if ax is None:
        fig,ax=plt.subplots()
        
    if normalize:
        E = tot_counts_norm(E, target_mean=1e6)[0]
    
    k = list(gene_list).index(gene_to_plot)
    coldat = E[:,k].A
    
    if smooth_operator is None:
        coldat = coldat.squeeze()
    else:
        coldat = np.dot(smooth_operator, coldat).squeeze()
    
    if order_points:
        o = np.argsort(coldat)
    else:
        o = np.arange(len(coldat))
        
    vmin = np.percentile(coldat, col_range[0])
    vmax = np.percentile(coldat, col_range[1])
    if vmax==vmin:
        vmax = coldat.max()
        
    pp = ax.scatter(x[o], y[o], c=coldat[o], s=point_size, cmap=color_map,
               vmin=vmin, vmax=vmax)
    
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlim(x.min()-x.ptp()*buffer_pct, x.max()+x.ptp()*buffer_pct)
    ax.set_ylim(y.min()-y.ptp()*buffer_pct, y.max()+y.ptp()*buffer_pct)
    
    return pp

def plot_gene_list(adata, genes_plot, embedding='X_umap', save_filename=None, close_after_save=False, dpi=75, n_columns=6, vmax=99.6, **plot_kwargs):
    x=adata.obsm[embedding][:,0]
    y=adata.obsm[embedding][:,1]


    fig,nrow,ncol = start_subplot_figure(len(genes_plot), row_height=2.5, n_columns=n_columns, fig_width=16 * n_columns / 6, dpi=dpi)
    axs = []

    for iG,g in enumerate(genes_plot):
        ax = plt.subplot(nrow, ncol, iG+1)
        axs.append(ax)
        plot_one_gene(
            adata.X, 
            adata.var_names.values.astype(str), 
            g.split(' ')[0], x, y, 
            ax=ax, 
            col_range=(0, vmax), 
            **plot_kwargs
        )
        ax.set_title(g)

    fig.tight_layout()
    
    if save_filename is not None:
        save_dir='/'.join(save_filename.rstrip('/').split('/')[:-1])
        os.makedirs(save_dir, exist_ok=True)
        plt.savefig(save_filename, dpi=dpi)
    if close_after_save:
        plt.close()
    else:
        return fig, axs

def plot_continuous(adata, var_name, embedding='X_umap', ax=None, order_points=False, 
                    col_range=(0,100), buffer_pct=0.03, point_size=1, color_map=None, 
                    log_transform=False, pseudocount=1):

    x = adata.obsm[embedding][:,0]
    y = adata.obsm[embedding][:,1]

    if color_map is None:
        color_map = darken_cmap(plt.cm.Reds,.9)
    if ax is None:
        fig,ax=plt.subplots()

    coldat = adata.obs_vector(var_name)
    #if var_name in adata.var_names:
    #    coldat = adata[:, var_name].X.squeeze().copy()
    #else:
    #    coldat = adata.obs[var_name].values.copy()
        
    if log_transform:
        coldat = np.log10(pseudocount + coldat)
        
    if order_points:
        o = np.argsort(coldat)
    else:
        o = np.arange(len(coldat))

    vmin = np.percentile(coldat, col_range[0])
    vmax = np.percentile(coldat, col_range[1])
    if vmax==vmin:
        vmax = coldat.max()
    
    pp = ax.scatter(x[o], y[o], c=coldat[o], s=point_size, cmap=color_map,
               vmin=vmin, vmax=vmax)

    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlim(x.min()-x.ptp()*buffer_pct, x.max()+x.ptp()*buffer_pct)
    ax.set_ylim(y.min()-y.ptp()*buffer_pct, y.max()+y.ptp()*buffer_pct)

    return pp

    
def plot_categorical(x, y, data, ax=None, buffer_pct=0.03, point_size=5, color_map=None, 
                     show_centroids=False, centroid_label_size=12, centroid_label_weight='normal',
                     centroid_label_align='center'):
    if color_map is None:
        color_map = plt.cm.Paired
    if ax is None:
        fig,ax=plt.subplots(figsize=(4,4))

        
    group_strings = data.copy()
    group_labels = np.unique(group_strings)
    group_nums = np.zeros(len(group_strings), dtype=int)
    centroids = np.zeros((len(group_labels), 2))
    
    for iG, g in enumerate(group_labels):
        group_nums[group_strings == g] = iG
        centroids[iG,0] = np.mean(x[group_strings==g])
        centroids[iG,1] = np.mean(y[group_strings==g])

    ax.scatter(x, y, c=group_nums, s=point_size, cmap=color_map)
    
    if show_centroids:
        for iLab, lab in enumerate(group_labels):
            ax.text(centroids[iLab,0], centroids[iLab,1], lab, fontsize=centroid_label_size, fontweight=centroid_label_weight, ha=centroid_label_align)

    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlim(x.min()-x.ptp()*buffer_pct, x.max()+x.ptp()*buffer_pct)
    ax.set_ylim(y.min()-y.ptp()*buffer_pct, y.max()+y.ptp()*buffer_pct)

    return

def start_subplot_figure(n_subplots, n_columns=5, fig_width=14, row_height=3, dpi=75):
    n_rows = int(np.ceil(n_subplots / float(n_columns)))
    fig = plt.figure(figsize = (fig_width, n_rows * row_height), dpi=dpi)
    return fig, n_rows, n_columns

def total_counts_histogram(total_counts, log_x=True, counts_per_bin=False, min_bin=10, max_bin=10e5, n_bins=50, ax=None):
    '''histogram of counts per barcode
    If counts_per_bin is True, the histogram is weighted, i.e., the value of y(x)
    is the total number of UMIs from barcodes with x UMIs.
    '''
    if log_x:
        bins = np.logspace(np.log10(min_bin), np.log10(max_bin), n_bins)
        xscale = 'log'
    else:
        bins = np.linspace(min_bin, max_bin, n_bins)
        xscale = 'linear'

    if ax is None:
        fig, ax = plt.subplots()

    if counts_per_bin:
        from collections import defaultdict
        count_freq = defaultdict(int)
        for count in total_counts:
            count_freq[count] += 1
        x = np.array(list(count_freq.keys()))
        y = np.array(list(count_freq.values()))
        w = x*y
        ax.hist(x, bins=bins, weights=w)
        ax.set_ylabel('# counts coming from bin')
    else:
        ax.hist(total_counts, bins=bins)
        ax.set_ylabel('# barcodes')

    ax.set_xscale(xscale)
    ax.set_xlabel('Counts per barcode')

    return

def custom_violinplot(data, ax, showmeans=False, showmedians=False, show_xticks=False, show_yticks=False, colors=None):
    parts = ax.violinplot(
            data, showmeans=False, showmedians=False,
            showextrema=False)
    
    if showmeans:
        means = [np.mean(d) for d in data]
        ax.scatter(np.arange(1, len(means)+1), means, c='black', s=10,marker='o')
    if showmedians:
        medians = [np.median(d) for d in data]
        ax.scatter(np.arange(1, len(medians)+1), medians, c='black', s=100, marker='_')


    for iP,pc in enumerate(parts['bodies']):
        if colors is None:
            pc.set_facecolor([.7,.7,.7])
        else:
            pc.set_facecolor(colors[iP])
        pc.set_edgecolor('none')
        pc.set_alpha(1)
        
    if not show_xticks:
        ax.set_xticks([])
    if not show_yticks:
        ax.set_yticks([])
        
    return parts


def plot_categ_bar(data, main_group_label, sub_group_label=None, plot_pct=True, main_x_offset=0.3, bar_width=None, fig_size=(14, 4), show_legend=True):
    main_group_vals = data[main_group_label].values.astype(str)
    x_groups = np.unique(main_group_vals)
    
    if sub_group_label is None:
        sub_group_label = 'all cells'
        sub_group_vals = np.repeat(sub_group_label, data.shape[0])
        x_subgroups = np.unique(sub_group_vals)
        sub_x_offsets = [0]
    else:
        sub_group_vals = data[sub_group_label].values.astype(str)
        x_subgroups = np.unique(sub_group_vals)
        sub_x_offsets = np.linspace(-main_x_offset, main_x_offset, len(x_subgroups))
    
    if bar_width is None:
        bar_width = 2 * main_x_offset / len(x_subgroups)

    plot_dat = {s: [] for s in sub_group_vals}
    denom = 1
    for iGroup, group in enumerate(x_groups):
        for iSub,sub in enumerate(x_subgroups):
            x = iGroup + sub_x_offsets[iSub]
            if plot_pct:
                denom = (sub_group_vals == sub).sum() * 0.01
            y = (main_group_vals[sub_group_vals == sub] == group).sum() / denom
            plot_dat[sub].append([x, y])
    for k,v in plot_dat.items():
        plot_dat[k] = np.array(v)

    fig, ax = plt.subplots(figsize=fig_size)
    for s in plot_dat:
        ax.bar(plot_dat[s][:, 0], plot_dat[s][:, 1], bar_width, label=s)

    if show_legend:
        ax.legend();
    ax.set_xticks(np.arange(len(x_groups)));
    ax.set_xticklabels(x_groups);
    
    if plot_pct:
        ax.set_ylabel('Percent')
    else:
        ax.set_ylabel('Count')

    fig.tight_layout()
    
    return fig, ax


