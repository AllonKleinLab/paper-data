import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import scrublet as scr
import os
import sys
import time
from statistics import median
from sklearn.neighbors import KNeighborsClassifier
from scipy.cluster.hierarchy import dendrogram, linkage

import helper_functions as hf

# Set random seed
np.random.seed(seed = 0)

# Set output path
out_path = 'figure2_plots_output/'
if not os.path.exists(out_path): os.mkdir(out_path)

# ============================================================================
# SET UP FOR SCRNA PLOTS

title_list = {
    'GSM8869530_Cr_blood_asw': 'dCMF-ASW',
    'GSM8869531_Cr_blood_mannitol': 'PBS-M'
}
color_list = {
    'GSM8869530_Cr_blood_asw': '#147dbe',
    'GSM8869531_Cr_blood_mannitol': '#ab004d'
}

# ------------------------------------
# Arguments for this filtering

# SET FILTERING PARAMETERS
# v1: UMI/bc thresholds used in preprint
version = 'v1'
# v2: UMI/bc thresholds used in Figure S1
#version = 'v2'
# v3: using 10X Cell Ranger's automatic barcode filtering
#version = 'v3'

if version == 'v1':
    arg_dict = {
        'version': 'v1',
        'aligned_genome': 'HT2019_KY21_with_Ens_mito',
        'max_mito_pct': [10, 10],   # order: dCMF-ASW, PBS-M
        'min_num_UMI': [2000, 800], # order: dCMF-ASW, PBS-M
        'min_num_genes': [0, 0],    # order: dCMF-ASW, PBS-M
        'run_scrublet': True
    }
elif version == 'v2':
    arg_dict = {
        'version': 'v2',
        'aligned_genome': 'HT2019_KY21_with_Ens_mito',
        'max_mito_pct': [10, 10],   # order: dCMF-ASW, PBS-M
        'min_num_UMI': [800, 800],  # order: dCMF-ASW, PBS-M
        'min_num_genes': [0, 0],    # order: dCMF-ASW, PBS-M
        'run_scrublet': True
    }
elif version == 'v3':
    arg_dict = {
        'version': 'v3',
        'aligned_genome': 'HT2019_KY21_with_Ens_mito',
        'max_mito_pct': [10, 10],   # order: dCMF-ASW, PBS-M
        'min_num_UMI': ['10x_list', '10x_list'],    # order: dCMF-ASW, PBS-M
        'min_num_genes': [0, 0],    # order: dCMF-ASW, PBS-M
        'run_scrublet': True
    }

# Update out_path
out_path = out_path + arg_dict['version'] + '/'
if not os.path.exists(out_path): os.mkdir(out_path)

# ------------------------------------
# Import count matrices

t = time.time()
print('Importing count matrices...')

# Get paths to folders, create folder for saving figures
data_path = 'raw_unfiltered_data/'
adict = {}

# Import counts matrix
filename_end = '_raw_unfiltered_matrix.h5'
libs = [d.split(filename_end)[0] for d in os.listdir(data_path)
        if d.endswith('_raw_unfiltered_matrix.h5')]
for lib in libs:
    adict[lib] = sc.read_10x_h5(data_path + lib + filename_end)

# Reverse order of libs so dCMF-ASW is first
libs = libs[::-1]

# ------------------------------------
# Mitochondrial gene labeling

for lib in libs:
    # Add a column to adict[lib].var that is a boolean for mitochondrial genes
    adict[lib].var['mt'] = adict[lib].var.index.str.contains('ENSCING')

# ------------------------------------
# Calculate summary statistics: number of UMIs per cell/barcode ("count
# depth"), proportion of genes/counts from mtDNA, and total genes per barcode.
for lib in libs:
    sc.pp.calculate_qc_metrics(adict[lib], qc_vars=['mt'], percent_top=None,
                               log1p=False, inplace=True)

# Before starting the filtration: remove all barcodes with 0 total counts (so
# we can visualize this in log space)
for lib in libs:
    adict[lib] = adict[lib][adict[lib].obs['total_counts'] > 0, :]

# ============================================================================
# 2d - UMIs/barcode

print('\nPlotting UMIs/barcode')
print('`````````````````````')

with plt.style.context('tal_paper'):
    # Set count thresholds (change in arg_dict above based on plots)
    for i in range(len(libs)):
        adict[libs[i]].uns['min_num_UMI'] = arg_dict['min_num_UMI'][::-1][i]

    ncol = 1
    nrow = len(adict)

    fig = plt.figure(figsize = (ncol * 2.5, nrow * 1))
    for i, lib in enumerate(libs[::-1]):
        ax0 = plt.subplot(nrow, ncol, i + 1)
        
        min_num_UMI = adict[lib].uns['min_num_UMI']
        
        # Plot histogram
        (freq, bins) = np.histogram(adict[lib].obs['total_counts'],
                                    np.logspace(1, 4.5, 50))
        ax0.bar(bins[:-1], freq*bins[:-1],
                width=0.9*np.diff(bins),
                color=color_list[lib], alpha=0.65)
        ax0.set_xscale('log')
        # ax0.set_ylim(0, 10000)

        # Include lines for mean UMIs/barcode in empty droplets, print values
        if isinstance(arg_dict['min_num_UMI'][i], int):
            barcode_passlist = (adict[lib].obs['total_counts']
                                > adict[lib].uns['min_num_UMI'])
        else:
            barcode_passlist_df = pd.read_csv(data_path + libs[i] + '_barcodes.tsv',
                                header=None)
            barcode_passlist = list(barcode_passlist_df.iloc[:, 0])
            barcode_passlist = adict[lib].obs.index.isin(barcode_passlist)
        empty = adict[lib].obs.loc[~barcode_passlist, 'total_counts']
        th = 10; empty = empty[empty > th]
        val_to_plot = median(empty)
        print(f'{lib} mean UMI/bc for empty droplets (excluding UMI/bc={th})'
              + f': {val_to_plot}')
        ax0.axvline(x=val_to_plot, linestyle='--', linewidth=1.5,
                    color='k')
        
        ax0.set_xlim(40,)
        # ax0.set_title(title_list[lib])
    ax0.set_xlabel('Num. unique transcripts per barcode\n'
                   + '(UMI-filtered mapped reads)')

    fig.tight_layout()
    plt.savefig(out_path + '2d.pdf')
    plt.close()

with plt.style.context('tal_paper'):
    # % reads in cells (reads which are confidently mapped to transcriptome)
    # Read in from the output of data_quality_metrics.py
    x = ['\n'.join(title_list[l].split(' ')) for l in title_list]
    y = []
    for lib in title_list:
        if isinstance(arg_dict['min_num_UMI'][i], int):
            barcode_passlist = (adict[lib].obs['total_counts']
                                > adict[lib].uns['min_num_UMI'])
        else:
            barcode_passlist_df = pd.read_csv(data_path + lib + '_barcodes.tsv',
                                              header=None)
            barcode_passlist = list(barcode_passlist_df.iloc[:, 0])
            barcode_passlist = adict[lib].obs.index.isin(barcode_passlist)
        empty = adict[lib].obs.loc[~barcode_passlist, 'total_counts']
        th = 10; empty = empty[empty > th]
        val_to_plot = median(empty)
        y.append(val_to_plot)
    # y = 100 * np.array(y)

    f, ax = plt.subplots(1, 1, figsize=(1.35, 2.25))
    # plt.bar(x, y, color='#999999')
    for i in range(len(y)):
        plt.bar(i, y[i], color=color_list[libs[::-1][i]], alpha=0.65, width=0.7)
    plt.xticks(ticks=np.arange(len(x)), labels=x, rotation=90)
    plt.xlim(plt.xlim()[0] - 0.15, plt.xlim()[1] + 0.15)
    plt.ylabel('Median UMI per barcode\nin empty droplets')

    plt.tight_layout()
    plt.savefig(out_path + 'empty_umi_per_bc.pdf')
    plt.close()

# ============================================================================
# 2e - Mitochondrial fraction (only barcodes which are cells)

print('\nPlotting mitochondrial fraction')
print('```````````````````````````````')

with plt.style.context('tal_paper'):
    ncol = 1
    nrow = 1#len(adict)

    if arg_dict['version'] == 'v1':
        fig = plt.figure(figsize = (ncol * 3.15, nrow * 2))
    else:
        fig = plt.figure(figsize = (ncol * 2.3, nrow * 2))
    for i, lib in enumerate(libs[::-1]):

        # Remove droplets passed the UMI threshold
        if isinstance(arg_dict['min_num_UMI'][i], int):
            barcode_passlist = (adict[lib].obs['total_counts']
                                > adict[lib].uns['min_num_UMI'])
        else:
            barcode_passlist_df = pd.read_csv(data_path + lib + '_barcodes.tsv',
                                              header=None)
            barcode_passlist = list(barcode_passlist_df.iloc[:, 0])
            barcode_passlist = adict[lib].obs.index.isin(barcode_passlist)
        to_plot = adict[lib][barcode_passlist]

        # Print number of cells with >20% mitochondrial fraction
        th = 20
        numer = np.sum(to_plot.obs["pct_counts_mt"] > th)
        denom = to_plot.shape[0]
        print(f'{lib}: {100*numer/denom:.1f}% ({numer}/{denom}) cells with '
              + f'mito fraction > {th}%')

        ax = plt.subplot(nrow, ncol, 1)#i + 1)

        # plot a histogram of the % of counts from mtDNA
        (freq, bins) = np.histogram(to_plot.obs['pct_counts_mt'], bins=50)
        freq = freq / np.sum(freq)
        ax.bar(bins[:-1], freq, width=np.diff(bins), color=color_list[lib],
            alpha=0.65, label=title_list[lib])

        ax.set_xlim(0,)
        ax.set_xlabel('Mitochondrial fraction (%)')
        ax.set_ylabel('Fraction of cells')
        ax.set_yscale('log')

    if arg_dict['version'] == 'v1':
        ax.legend(loc='center right')
    else:
        ax.legend(loc='upper right')
    fig.tight_layout()
    plt.savefig(out_path + '2e.pdf')
    plt.close()

# ============================================================================
# 2c - Fraction of reads from droplets with cells

with plt.style.context('tal_paper'):
    # % reads in cells (reads which are confidently mapped to transcriptome)
    # Read in from the output of data_quality_metrics.py
    x = ['\n'.join(title_list[l].split(' ')) for l in title_list]
    y = []
    for l in title_list:
        df = pd.read_csv(f'data_quality_metrics/{arg_dict["version"]}/'
                         + f'{l}_quality_stats.csv')
        num = df.loc['Reads confidentally mapped to transcriptome in droplets with cells', 'MAPPING']
        num = int(num.split(' ')[0])
        denom = df.loc['Valid reads mapped confidently to transcriptome', 'MAPPING']
        denom = int(denom)
        y.append(100 - (100 * num / denom))
    # y = 100 * np.array(y)

    f, ax = plt.subplots(1, 1, figsize=(1.35, 2.25))
    # plt.bar(x, y, color='#999999')
    for i in range(len(y)):
        plt.bar(i, y[i], color=color_list[libs[::-1][i]], alpha=0.65, width=0.7)
        plt.text(i, y[i]+2, "{:.1f}%".format(y[i]),
                 horizontalalignment='center', fontweight='bold')
    plt.xticks(ticks=np.arange(len(x)), labels=x, rotation=90)
    plt.xlim(plt.xlim()[0] - 0.15, plt.xlim()[1] + 0.15)
    plt.ylabel('Fraction of reads in\nempty droplets (%)')
    plt.ylim(0, 100)

    plt.tight_layout()
    plt.savefig(out_path + '2c.pdf')
    plt.close()


# ============================================================================
# IMPORT PROCESSED DATA

# Import filtered & preprocessed adata
version = arg_dict['version']
adata_path = (f'filtering_and_preprocessing_output/{version}/')
adata = sc.read_h5ad(adata_path + 'adata.h5ad')

# ============================================================================
# 2b - Number of cells recovered

with plt.style.context('tal_paper'):
    # % reads in cells (reads which are confidently mapped to transcriptome)
    x = ['\n'.join(title_list[l].split(' ')) for l in title_list]
    y = [
        np.sum(adata.obs['sample'] == 'GSM8869530_Cr_blood_asw'),
        np.sum(adata.obs['sample'] == 'GSM8869531_Cr_blood_mannitol')
    ]

    f, ax = plt.subplots(1, 1, figsize=(1.25, 2.15))
    for i in range(len(y)):
        plt.bar(i, y[i], color=color_list[libs[::-1][i]], alpha=0.65, width=0.7)
    plt.axhline(y = 6000, color = 'k', alpha=0.5, linestyle = '--',
                label='Expected cell number')
    plt.xticks(ticks=np.arange(len(x)), labels=x, rotation=90)
    plt.xlim(plt.xlim()[0] - 0.15, plt.xlim()[1] + 0.15)
    plt.ylabel('Number of cells recovered')

    plt.tight_layout()
    plt.savefig(out_path + '2b.pdf')
    plt.close()

# ============================================================================
# 2f - UMAP with samples highlighted

print('\nUMAP plots')
print('``````````')

with plt.style.context('tal_paper_spine'):
    hf.pl_umap_separate(adata, 'sample', alpha=0.6,
                        # title=[title_list[s] for s in adata.obs['sample'].cat.categories],
                        title=['' for s in adata.obs['sample'].cat.categories],
                        palette=[color_list[s] for s in adata.obs['sample'].cat.categories])
    plt.gcf().set_size_inches(4, 2)
    plt.tight_layout()
    plt.savefig(out_path + '2f.png', dpi=300)
    plt.close()

# ============================================================================
# 2g - UMAP density plot

num_neighbors = 200

with plt.style.context('tal_paper_spine'):
    f, ax = plt.subplots(1, 1)
    f = hf.sample_density(
        adata,
        'sample',
        [(
            'GSM8869531_Cr_blood_mannitol',
            'GSM8869530_Cr_blood_asw'
        )],
        num_neighbors=num_neighbors,
        cmap=hf.cmap_pink2blue_r,
        sort_order=False,
        log_offset=1,
    )
    plt.close()

    # Plot with more control
    obs = 'GSM8869531_Cr_blood_mannitol_over_GSM8869530_Cr_blood_asw_density'
    vmax = max(
        abs(np.percentile(adata.obs[obs], 99)),
        abs(np.percentile(adata.obs[obs], 1))
    )
    sc.pl.umap(adata, color=obs, cmap=hf.cmap_pink2blue_r, sort_order=False,
               show=False, vmax=vmax, vmin=-vmax)
    plt.title('')
    plt.gcf().set_size_inches(2.25, 2)
    plt.tight_layout()
    plt.savefig(out_path + '2g.png', dpi=300)
    plt.close()

# ============================================================================
# 2h,i - Which cell states are enriched/depleted in PBS-M vs. dCMF-ASW?

# ------------------------------------
# Import cell state labels for PBS-M sample
# Annotated dataset from Scully et al. 2025, bioRxiv
#   Paper: https://www.biorxiv.org/content/10.1101/2025.05.20.655184v1
#   GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE296253
# Note that the combined and processed data h5ad file must be downloaded
# to this repo, stored in the directory Scully_Ciona_blood_2025/data/
# See Scully_Ciona_blood_2025/data/README.md for more info.
# Alternatively, download to a different folder and adjust the path below.

# Import annotated atlas
repo_dir = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))))
full_atlas_dataset = os.path.join(
    repo_dir,
    'Scully_Ciona_blood_2025',
    'data',
    'GSE296253_combined_processed_data.h5ad'
)
ad_ref = sc.read_h5ad(full_atlas_dataset)

# Add cell type information from full dataset
annotated_pbsm_bc = ad_ref.obs.loc[ad_ref.obs['library']=='220428_1', 'cell_type']
annotated_pbsm_bc.index = [x.replace('220428_1', 'GSM8869531_Cr_blood_mannitol')
                        for x in annotated_pbsm_bc.index]
adata.obs.loc[annotated_pbsm_bc.index, 'cell_type'] = annotated_pbsm_bc.astype('str')
cell_type_colors = ad_ref.uns['cell_type_colors']
del ad_ref

# Using Table S6 of Scully et al. to annotate HA-3/SRC doublets
with open('Scully_et_al_2025_atlas/table_s6_doublet_bc_list.csv', 'r') as f:
    doublet_bc_list = []
    l = f.readline()
    while l != '':
        l = f.readline()
        if '220428_1' in l:
            bc = l.strip().replace('220428_1', 'GSM8869531_Cr_blood_mannitol')
            doublet_bc_list.append(bc)
adata.obs['cell_type'] = adata.obs['cell_type'].astype('str')
adata.obs.loc[doublet_bc_list, 'cell_type'] = 'HA-3/SRC doublet'

# ------------------------------------
# Use annotations from PBS-M sample to annotate other cells
# (written and adjusted from ChatGPT)

# Mask for labeled and unlabeled cells
labels = adata.obs['cell_type'].copy()
labeled_mask = labels != 'nan'
unlabeled_mask = ~labeled_mask

X = adata.obsm["X_pca"]
X_labeled = X[labeled_mask]
y_labeled = labels[labeled_mask]
X_unlabeled = X[unlabeled_mask]

# Fit KNN classifier on labeled cells
knn = KNeighborsClassifier(n_neighbors=10)
knn.fit(X_labeled, y_labeled)

# Predict labels for unlabeled cells
y_pred = knn.predict(X_unlabeled)

# Assign predictions back into obs
adata.obs.loc[unlabeled_mask, 'cell_type'] = y_pred

# Set coloring to match Scully et al. 2025
adata.obs['cell_type'] = adata.obs['cell_type'].astype('str').astype('category')
adata.uns['cell_type_colors'] = (
    list(cell_type_colors[:5])
    + ["#caa5c9",]
    + list(cell_type_colors[5:])
)

# Plot
with plt.style.context('tal_paper_spine'):
    sc.pl.umap(adata, color='cell_type', show=False)
    plt.title('')

    plt.gcf().set_size_inches(2, 2)
    plt.tight_layout()
    plt.savefig(out_path + '2h.pdf', dpi=300)
    plt.close()

# ------------------------------------
# Fraction of cells in each sample

df = pd.DataFrame(
    index=[c for c in adata.obs['cell_type'].cat.categories
           if c!='HA-3/SRC doublet'],
)
adata_dasw = adata[adata.obs['sample'] == 'GSM8869530_Cr_blood_asw']
adata_pbsm = adata[adata.obs['sample'] == 'GSM8869531_Cr_blood_mannitol']
# Exclude HA-3/SRC doublet cluster for this analysis
adata_dasw = adata_dasw[adata_dasw.obs['cell_type'] != 'HA-3/SRC doublet']
adata_pbsm = adata_pbsm[adata_pbsm.obs['cell_type'] != 'HA-3/SRC doublet']

for c in df.index:
    df.loc[c, 'p_dasw'] = (np.sum(adata_dasw.obs['cell_type'] == c)
                           / adata_dasw.shape[0])
    df.loc[c, 'p_pbsm'] = (np.sum(adata_pbsm.obs['cell_type'] == c)
                           / adata_pbsm.shape[0])
    df.loc[c, 'cell_state_size'] = np.sum(adata.obs['cell_type'] == c)
df = df.astype('float')

# Ignore cell states with <5 cells total
df = df[df['cell_state_size'] >= 5]

# Log foldchange of proportions
df['logfc'] = np.log2(df['p_pbsm'] / df['p_dasw'])
df = df.sort_values(by='logfc')

# Error bars, assuming binomial distribution
"""
HOW ERROR BARS WERE CALCULATED

Variable definitions:
p_pbsm = proportion of PBS-M cells in this cell state
p_dasw = proportion of dCMF-ASW cells in this cell state
N_pbsm = total # cells in PBS-M
N_dasw = total # cells in dCMF-ASW

Standard error of a proportion:
sig_p_pbsm = sqrt(p_pbsm * (1-p_pbsm) / N_pbsm)
sig_p_dasw = sqrt(p_dasw * (1-p_dasw) / N_dasw)

Log2 fold change:
f(p_pbsm, p_dasw) = log2(p_pbsm / p_dasw)
df/dp_pbsm = 1 / (ln(2) * p_pbsm)
df/dp_dasw = 1 / (ln(2) * p_dasw)
sig_f = ln(2) * sqrt((sig_p_pbsm^2 / p_pbsm^2) + (sig_p_dasw^2 / p_dasw^2))
"""
N_pbsm = adata_pbsm.shape[0]
N_dasw = adata_dasw.shape[0]
df['p_pbsm_err'] = np.sqrt(df['p_pbsm'] * (1-df['p_pbsm']) / N_pbsm)
df['p_dasw_err'] = np.sqrt(df['p_dasw'] * (1-df['p_dasw']) / N_dasw)
df['logfc_err'] = np.log(2) * np.sqrt(
    (df['p_pbsm_err']**2 / df['p_pbsm']**2)
    + (df['p_dasw_err']**2 / df['p_dasw']**2)
)

# Statistical significance (two-tailed Z-test)
from scipy.stats import norm
df['z_score'] = (df['logfc'] - 0) / df['logfc_err']
df['p_val'] = 2 * (1 - norm.cdf(abs(df['z_score'])))
# Bonferroni correction
df['p_val_adj'] = df['p_val'] / df.shape[0]
df['null_rejected'] = (df['p_val_adj'] < 0.01)

# Plot
with plt.style.context('tal_paper_spine'):
    plt.figure(figsize=(3.75, 2))
    plt.bar(
        x=np.arange(df.shape[0]),
        height=df['logfc'],
        yerr=df['logfc_err'] * 1.96,
        color="#888888",
        capsize=1,
        error_kw={'elinewidth': 1}
    )

    plt.axhline(y=0, color='k', linewidth=0.75)
    plt.xticks(ticks=np.arange(df.shape[0]), labels=df.index, rotation=90)
    # plt.grid(axis='x')
    plt.xlabel('Cell state')
    plt.ylabel('Log2 fold change\n(PBS-M / dCMF-ASW)')
    plt.xlim(-1, df.shape[0])
    
    plt.tight_layout()
    plt.savefig(out_path + '2i.pdf')
    plt.close()

# ============================================================================
# Transcriptional differences between dCMF-ASW and PBS-M

adata.obs['cell_type_sample'] = [
    adata.obs.loc[i, 'cell_type'] + ' ' + adata.obs.loc[i, 'sample']
    for i in adata.obs.index
]
markers_dict = {}

for this_state in adata.obs['cell_type'].cat.categories:
    # this_state = 'HA-3'
    pbsm_state = f'{this_state} GSM8869531_Cr_blood_mannitol'
    dasw_state = f'{this_state} GSM8869530_Cr_blood_asw'

    if dasw_state in adata.obs['cell_type_sample'].values:

        # Marker genes
        sc.tl.rank_genes_groups(adata, groupby='cell_type_sample', method='wilcoxon',
                                groups=[pbsm_state], reference=dasw_state)
        df = sc.get.rank_genes_groups_df(adata, pbsm_state)
        markers = df.loc[
            (np.abs(df['logfoldchanges']) > 1) * (df['pvals_adj'] < 0.05),
            'names'
        ].values

        # Centroids of clusters
        centr = pd.DataFrame(index=adata.obs['cell_type_sample'].cat.categories,
                            columns=markers)
        for state in centr.index:
            these_cells = adata.obs['cell_type_sample'] == state
            centr.loc[state] = np.mean(adata[these_cells, markers].X, axis=0)
        centr = (centr - centr.mean()) / centr.std()
        centr = centr.loc[[pbsm_state, dasw_state], :]
        centr = centr.astype('float')

        # ------------------------------------
        # Order genes etc. for plotting

        centr = centr.sort_values(axis=1, by=dasw_state, ascending=False)

        # # Hierarchically cluster genes
        # with plt.style.context('tal_paper'):
        #     f = plt.figure(figsize=(3, 3))
        #     ax = plt.subplot(1, 1, 1)
        #     linkage_data = linkage(centr[markers].T, method='single',
        #                            metric='cosine', optimal_ordering=True)
        #     dendr = dendrogram(linkage_data, labels=centr[markers].columns,
        #                        count_sort='descending', color_threshold=0.00532)
        #     # plt.tight_layout()
        #     # plt.savefig(out_path + f'{"-".join(cl_ordered)}_gene_dendrogram.pdf')
        #     plt.close()
        # gene_order = np.array(dendr['ivl'][::-1])
        # cluster_arr = np.array(dendr['leaves_color_list'])
        # # gene_order_cluster_specific = gene_order[cluster_arr != 'C1']

        # ------------------------------------
        # Plot expression for each cluster and animal

        with plt.style.context('tal_paper_spine'):
            f, ax = plt.subplots(1, 1, figsize=(1.5, 3.5))

            # f = plt.figure(figsize=(1, 1.75))#clustermtx.shape[1]/35))
            # ax = plt.subplot(1, 1, 1)
            sns.heatmap(centr.T, cmap='RdBu_r', center=0, vmax=2, vmin=-2,
                        cbar_kws={'label': 'Z-score across\ncell types'},
                        xticklabels=True, yticklabels=False, ax=ax)
            plt.ylabel('Marker genes')
            ax.spines.top.set_visible(True)
            ax.spines.right.set_visible(True)
            ax.spines.left.set_visible(True)
            ax.spines.bottom.set_visible(True)

            # Plot dendrogram
            # from scipy.cluster import hierarchy
            # Z = hierarchy.linkage(clustermtx, 'single')
            # dn = hierarchy.dendrogram(Z)

            # Save
            plt.tight_layout()
            plt.savefig(out_path + f'expr_heatmap_{this_state.replace("/", "-")}.pdf')
            plt.close()

        markers_dict[this_state] = markers

# Markers shared across multiple cell states
markers_dict = {c: list(markers_dict[c]) for c in markers_dict}

# Count how many lists each gene appears in
from collections import Counter
counts = Counter()
for genes in markers_dict.values():
    counts.update(set(genes))  # use set() so duplicates within a list are not double counted

# Threshold: at least half the lists
threshold = len(markers_dict) / 3

# Get genes meeting the threshold
result = [gene for gene, c in counts.items() if c >= threshold]

for g in result:
    g = g.replace('KY21:', '')
    if g in hf.ciona2human:
        print(g, hf.ciona2human[g])

centr = pd.DataFrame(index=adata.obs['cell_type_sample'].cat.categories,
                     columns=result)
for state in centr.index:
    these_cells = adata.obs['cell_type_sample'] == state
    centr.loc[state] = np.mean(adata[these_cells, result].X, axis=0)
centr = (centr - centr.mean()) / centr.std()
centr = centr.astype('float')
