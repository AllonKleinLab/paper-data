import os

import numpy as np
import scipy.stats
import scipy.sparse as sp
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

### loading counts

def file_opener(filename):
    '''Open file and return a file object, automatically decompressing zip and gzip 
    Arguments 
    
    S. Wolock
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

def load_annotated_text(file_data, delim='\t', read_row_labels=False, read_column_labels=False, transpose=False, chunk_size=2000):
    '''Load text file as scipy.sparse.csc_matrix, returning column and/or row labels if desired.
    Loads rows in chunks to ease memory demands.
    
    S. Wolock
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
                    E_chunks.append(sp.coo_matrix((X_data, (X_row, X_col)), dtype=float, shape=(nrow, ncol)))
                    X_data = []
                    X_row = []
                    X_col = []
                    nrow = 0

    if nrow > 0:
        E_chunks.append(sp.coo_matrix((X_data, (X_row, X_col)), dtype=float, shape=(nrow, ncol)))

    E = sp.vstack(E_chunks)
    if transpose: 
        E = E.T
    
    return E.tocsc(), np.array(row_labels), np.array(column_labels)

########## CALCULATE COMPOSITE GENE SCORES
def calc_score_denserank(adata, score_name, score_gene_list):
    ''' Calculate a composite score per cell from a set of marker genes
    Arguments
    
    H.M. Cambra
    
    - adata : annData object
        containing the single cell data structure
    - score_gene_list : list of str
        the gene symbols to be used for calculating a score
        
    This code adds (or updates) two dictionaries:
        adata.uns['scores'][score_name] contains an np array of size (num cells) x 1 with the scores
        adata.uns['score_genes'][score_name] contains a list of gene symbols used to generate the score
    
    The score is calculated as the (unweighted) average of normalized dense ranks of the genes across cells.
    Let x_{ij} be the expression of gene j in cell i
      and r_{ij} = denserank(x_ij) over cells
      then the score s_i = (1/nGenes) (sum_{j}  r_{ij}/max_i(r_{ij}))
    '''
    
    score = np.zeros(adata[:,1].shape) # Initialize the score for each cell to zero
    
    for gene in score_gene_list: # Iterate over all genes to form the score
    	ranks = scipy.stats.rankdata(scipy.sparse.csc_matrix(adata[:,gene].X).todense(),method='dense') # Calculate the dense rank
    	ranks = ranks.reshape([len(ranks),1])
    	score = score + ranks/np.max(ranks) # Normalize the ranks to max=1, and add to the score
    
    score = np.squeeze(score/(len(score_gene_list)));
    
    score_n = (score - np.mean(score))/(np.std(score)) # zscore the score

    # Save the score as a new unstructured annData dictionary    
    adata.obs[f'gene_score_z_{score_name}'] = score_n
    
    # Save the score as a new unstructured annData dictionary    
    adata.obs[f'gene_score_{score_name}'] = score
    
    if f'score_genes_{score_name}' not in adata.uns:
    	adata.uns[f'score_genes_{score_name}'] = {}
    	adata.uns[f'score_genes_{score_name}'][score_name] = score_gene_list

########## CLASSIFY CELLS BASED ON GENE SCORES

def classify_gene_score(adata, gene_score_prefix):
	'''
	Requires running 'calc_score_denserank' first to generate gene scores for each state based on the number of 
    marker genes chosen for that state.
    
    H.M. Cambra
    
    Arguments:
    ----------
    adata: AnnData object (scanpy)
        AnnData object to be classified
    gene_score_prefix: string
    	string identifier for the gene score to use to classify the data; should not include the number of genes used
    ngenes: int
        number of marker genes used to calculate the gene score for a state; should match already calculated gene 
        score. 
    Returns:
    --------
    adata: AnnData object (scanpy)
        modifided AnnData object with newly added obs DataFrame with classified states based on marker gene scores.
	
	'''
	gene_score_dict = {}
	
	for groups in adata.obs:
		if groups.startswith(f'{gene_score_prefix}'):
			gene_score_dict[groups.split('_')[4]] = adata.obs[groups]
	gene_score_df = pd.DataFrame.from_dict(gene_score_dict)
	assign_state = gene_score_df.idxmax(axis=1)
	adata.obs[f'classify_{gene_score_prefix}']=np.array(assign_state)
	
########## ABUNDANCE ANALYSIS

def abundance_df(adata, obs_group):
	
	'''
	Takes in AnnData object and creates abundance dataframe with a column for each obs_group
	and an abundance value for each state. 
	
	H.M. Cambra
	
	Arguments:
	----------
	adata: AnnData object (scanpy)
        Object containing counts matrix (AnnData.X), unnormalized counts (AnnData.raw.X), other matrices associated
        with data processing (i.e. umap coordinates, PCs, neighbors), 
        and metadata (i.e. gene lists, sample information) 
	
	obs_group: string
		string name of pandas series of obs dataframe of AnnData object
		specifies which column you'd like to subset data on for each state
	
	Returns:
	--------
	abund_df_ctype: pandas dataframe
		abundances of each state in each obs_group where all cells are represented in 
		obs_group and a subset in each state
	'''
	
	# generate dictionary of abundances by each category in obs group (cell state)
	abund_dict = {}
	for state_name in adata.obs[f'{obs_group}'].unique():
		abund_dict[state_name] = len(adata.obs.loc[adata.obs[f'{obs_group}'] == f'{state_name}'])
	
	# generate dataframe of abundances
	abund_df = pd.DataFrame.from_dict(abund_dict, orient='index')
	abund_df = abund_df.reset_index()
	abund_df = abund_df.rename(columns={'index' : 'state', 0 : f'{obs_group}_abundance'})
	
	# get normalized abundance levels for each state
	abund_df_ctype = abund_df
	total_df_ctype = abund_df_ctype[f'{obs_group}_abundance'].sum(axis=0)
	abund_df_ctype[f'{obs_group}_abundance'] = abund_df_ctype[f'{obs_group}_abundance']/total_df_ctype
	abund_df_ctype[f'{obs_group}_log2_abundance'] = np.log2(abund_df_ctype[f'{obs_group}_abundance'])
	
	return abund_df_ctype

########## PLOTTING

def set_plot_defaults(fontsize=10):
'''
Set default parameters for plotting.
'''
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = 'Arial'
    plt.rc('font', size=fontsize)
    plt.rcParams['pdf.fonttype'] = 42
    return

def start_subplot_figure(n_subplots, n_columns=5, fig_width=14, row_height=3, dpi=75):
    n_rows = int(np.ceil(n_subplots / float(n_columns)))
    fig = plt.figure(figsize = (fig_width, n_rows * row_height), dpi=dpi)
    return fig, n_rows, n_columns

def total_counts_histogram(total_counts, log_x=True, counts_per_bin=False, min_bin=10, max_bin=10e5, n_bins=50, ax=None, color='g'):
    '''histogram of counts per barcode
    If counts_per_bin is True, the histogram is weighted, i.e., the value of y(x)
    is the total number of UMIs from barcodes with x UMIs.
    
    S. Wolock
    
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
        ax.hist(x, bins=bins, weights=w, color=color)
        ax.set_ylabel('# counts coming from bin')
    else:
        ax.hist(total_counts, bins=bins, color=color)
        ax.set_ylabel('# barcodes')

    ax.set_xscale(xscale)
    ax.set_xlabel('Counts per barcode')

    return

def plot_abundance(x, y, group, corr_method='Pearson', show_corr=True):
    '''
    Plot the abundance values for each state in a merged dataframe of abundance levels for different sets of data. 
    
    H.M. Cambra
    
    Arguments:
    ----------
    x: col of dataframe
        set of abundance levels for x-axis of plot; represents one set of abundance levels
    y: col of dataframe
        set of abundance levels for y-axis of plot; represents another set of abundance levels
    group: col of dataframe
        abundance level group - usually the state variable.
    corr_method: string
        scipy.stats method to calculate the correlation of abundance between x levels and y levels
        Options: Pearson or Spearman
    show_corr: bool
    	if True, shows the correlation method on the plot
        
    Returns:
    --------
    plot: seaborn scatterplot object
        scatterplot of abundance levels for each group with best fit line and Pearson or Spearman correlation.
        
    '''
    
    fig = plt.figure(figsize=(15, 15))
    sns.set(style='white', font_scale=2.0, font='Arial')
    
    plot = sns.scatterplot(x=x, y=y, hue=group, s=120)
    
    plt.plot(np.unique(x.values), np.poly1d(np.polyfit(x.values, y.values, 1))(np.unique(x.values)))
    
    if corr_method == 'Pearson':
    	corr = scipy.stats.pearsonr(x, y)[0]
    elif corr_method == 'Spearman':
    	corr = scipy.stats.spearmanr(x, y)[0]
    	
    corr = np.around(corr, decimals=2)
    
    if show_corr:
    	plt.text(x.max() + y.min()/2, y.min(), f'{corr_method} correlation = {corr}')
    	
    plt.legend(loc='best', prop={'size': 8}, markerscale=1.0, fontsize=18)
    
    plt.tight_layout()
    
    return plot
    

def plot_groups(x, y, groups, buffer_pct=0.03, saving=False, fig_dir='./', fig_name='fig', res=300, close_after=False, title_size=12, point_size=3, ncol=5, fig_width=14, row_height=3):
    
    ''' 
    S. Wolock
    '''

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


def plot_processing(adata, coordinates, group, x=0, y=1, fig_dir='./', fig_name='group_plot', coordinates_pca=False, coordinates_umap=True, saving=False):
    '''
    
    Calls def plot_groups function with annotations in either PCA or UMAP space.  
    
    H.M. Cambra
    
    Arguments:
    ----------
    adata: AnnData object (scanpy)
        Object containing counts matrix (AnnData.X), unnormalized counts (AnnData.raw.X), other matrices associated
        with data processing (e.g. umap coordinates, PCs, neighbors), 
        and metadata (e.g. gene lists, sample information)
    coordinates: numpy array
        matrix of reduced dimensions of 2D plot of shape ncells x 2
    x: int
        first axis to plot (e.g. first dimension of UMAP coodinates)
    y: int
        second axis to plot (e.g. second dimension of UMAP coodinates)
    group: string
        string name of pandas series representing each cell along the same index as coordinates matrix
        (e.g. Series in pandas dataframe given by the AnnData object's obs structure)
        
    Returns:
    --------
    plot: matplolib plot object
        plots group associated of cells in 2D coordinates given 
    
    '''
    
    if coordinates_umap:
        fig_prefix = 'UMAP_coord_{}'.format(group)
        fig = plot_groups(
        coordinates[:, x],
        coordinates[:, y],
        adata.obs[f'{group}'].values.astype(str),
        saving=False, 
        fig_dir=fig_dir, 
        fig_name=fig_prefix+fig_name)
    
    elif coordinates_pca:
        fig_prefix = '{}_PCs_PC{}_by_PC{}_{}'.format(adata.obsm['X_pca'].shape[1], x, y, group)
        plot_groups(
            coordinates[:, x],
            coordinates[:, y],
            adata.obs[f'{group}'].values.astype(str),
            saving=False, 
            fig_dir=fig_dir, 
            fig_name=fig_prefix+fig_name)
    
    plt.suptitle('Group plot by {}\nNumber of neighbors: {}\nNumber of PCs: {}\nVisualization: {}'.format(group, list(adata.uns['neighbors']['params'].values())[0], 
                              adata.obsm['X_pca'].shape[1], 
                              list(adata.uns['neighbors']['params'].values())[1]), y=0.01)
    plt.show()
    