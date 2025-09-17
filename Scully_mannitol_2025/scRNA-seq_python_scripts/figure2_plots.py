import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
import scrublet as scr
import os
import sys
import time
from statistics import median
from tqdm import tqdm

import helper_functions as hf

# Change this path to point to folder containing tal_helper_functions.py
path_to_dropbox = os.environ['PATH_TO_DROPBOX']

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

# Filtration conditions
arg_dict = {
    'version': 'vx',
    'max_mito_pct': '10,10',
    'min_num_UMI': '2000, 800',
    'min_num_genes': '0,0',
    'run_scrublet': 'True'
}

# Format arguments with commas into a list of ints
for a in arg_dict:
    if ',' in arg_dict[a]:
        arg_dict[a] = [int(x) for x in arg_dict[a].split(',')]
    elif 'True' == arg_dict[a] or 'False' == arg_dict[a]:
        arg_dict[a] = (arg_dict[a] == 'True')

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
        if d.endswith('.h5')]
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
        empty = adict[lib].obs.loc[
            adict[lib].obs['total_counts'] < adict[lib].uns['min_num_UMI'],
            'total_counts']
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

# ============================================================================
# 2e - Mitochondrial fraction (only barcodes which are cells)

print('\nPlotting mitochondrial fraction')
print('```````````````````````````````')

with plt.style.context('tal_paper'):
    N = {}
    for i in range(len(libs)):
        N[libs[i]] = arg_dict['min_num_UMI'][::-1][i]

    ncol = 1
    nrow = 1#len(adict)

    fig = plt.figure(figsize = (ncol * 3.15, nrow * 2))
    for i, lib in enumerate(libs[::-1]):

        # Remove droplets passed the UMI threshold
        to_plot = adict[lib][adict[lib].obs.total_counts > N[lib]]

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

    ax.legend(loc='center right')
    fig.tight_layout()
    plt.savefig(out_path + '2e.pdf')
    plt.close()

# ============================================================================
# IMPORT PROCESSED DATA

# Import filtered & preprocessed adata
version = 'v1'
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
# 2c - Fraction of reads from droplets with cells

with plt.style.context('tal_paper'):
    # % reads in cells (reads which are confidently mapped to transcriptome)
    x = ['\n'.join(title_list[l].split(' ')) for l in title_list]
    y = [49306418/145747789, 108362126/145988546]
    y = 100 * np.array(y)

    f, ax = plt.subplots(1, 1, figsize=(1.35, 2.25))
    # plt.bar(x, y, color='#999999')
    for i in range(len(y)):
        plt.bar(i, y[i], color=color_list[libs[::-1][i]], alpha=0.65, width=0.7)
        plt.text(i, y[i]+2, "{:.1f}%".format(y[i]),
                 horizontalalignment='center', fontweight='bold')
    plt.xticks(ticks=np.arange(len(x)), labels=x, rotation=90)
    plt.xlim(plt.xlim()[0] - 0.15, plt.xlim()[1] + 0.15)
    plt.ylabel('Fraction of reads in\ndroplets with cells (%)')
    plt.ylim(0, 100)

    plt.tight_layout()
    plt.savefig(out_path + '2c.pdf')
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

try:
    adata.obs['cell_type'] = pd.read_csv(out_path + 'cell_type_annotation.csv',
                                         index_col=0)
    with open(out_path + 'cell_type_colors.txt', 'r') as f:
        adata.uns['cell_type_colors'] = [x.strip() for x in f.readlines()]

except:
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
    adata.obs.loc[annotated_pbsm_bc.index, 'cell_type'] = annotated_pbsm_bc
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

    # Look at 20 nearest neighbors
    sc.pp.neighbors(adata, n_neighbors=20, n_pcs=50, use_rep='X_pca')

    unlabeled_bc = adata.obs.index[adata.obs['cell_type'] == 'nan']
    count = 1

    while len(unlabeled_bc) > 0:
        print(f'Round {count}: {len(unlabeled_bc)} unlabeled barcodes')

        # Count cell state annotations in neighborhood
        neighborhood_counts = hf.neighborhood_composition(
            adata,
            'cell_type',
            use_existing_neighbor_graph=True,
            cell_subset=unlabeled_bc
        )
        del neighborhood_counts['nan']  # don't look at unlabeled neighbors
        bc_has_labeled_neighbors = (neighborhood_counts.sum(axis=1) > 0)

        # Label cell states based on highest number of neighbor labels
        for bc in unlabeled_bc:
            these_neighbors = neighborhood_counts.loc[bc, :]
            if bc_has_labeled_neighbors[bc]:
                adata.obs.loc[bc, 'cell_type'] = these_neighbors.idxmax()
        
        count += 1
        unlabeled_bc = adata.obs.index[adata.obs['cell_type'] == 'nan']

    # Save
    adata.obs['cell_type'].to_csv(out_path + 'cell_type_annotation.csv')
    with open(out_path + 'cell_type_colors.txt', 'w') as f:
        for c in adata.uns['cell_type_colors']: f.write(c + '\n')


# Plot
with plt.style.context('tal_paper_spine'):
    sc.pl.umap(adata, color='cell_type', show=False)
    plt.title('')

    plt.gcf().set_size_inches(2, 2)
    plt.tight_layout()
    plt.savefig(out_path + '2h.pdf', dpi=300)
    plt.close()

# ------------------------------------
# Fraction of cells in dCMF-ASW sample

df = pd.DataFrame(
    index=[c for c in adata.obs['cell_type'].cat.categories
           if c!='HA-3/SRC doublet'],
    columns=['num_dasw', 'num_pbsm', 'num_cells']
)
for c in df.index:
    samples = adata.obs.loc[adata.obs['cell_type'] == c, 'sample']
    df.loc[c, 'num_dasw'] = np.sum(samples.values=='GSM8869530_Cr_blood_asw')
    df.loc[c, 'num_pbsm'] = np.sum(samples.values=='GSM8869531_Cr_blood_mannitol')
    df.loc[c, 'num_cells'] = len(samples)
df = df.astype('int')

# Ignore cell states with <5 cells total
df = df[df['num_cells'] >= 5]
import numpy as np
df['logfc'] = np.log2(df['num_pbsm'] / df['num_dasw'])
df = df.sort_values(by='logfc')

# Error bars, assuming binomial distribution
# n = num PBS-M cells
# N = total num cells
# (N-n = num dCMF/ASW cells)
# p = n / N
# f = log2(n / (N-n)) = log2(p / (1-p))
# sigma_f = (1 / ln(2)) * sqrt(1 / (N * p * (1-p)))
overall_n = df['num_pbsm'].sum()
overall_N = df['num_pbsm'].sum() + df['num_dasw'].sum()
overall_p = overall_n / overall_N
overall_fc = np.log2(overall_p / (1-overall_p))
overall_fc_err = (
    (1 / np.log(2)) 
    * np.sqrt(1 / (overall_N * overall_p * (1-overall_p)))
)

df['n'] = df['num_pbsm']
df['N'] = df['num_pbsm'] + df['num_dasw']
df['p'] = df['n'] / df['N']
df['logfc_err'] = (
    (1 / np.log(2)) 
    * np.sqrt(1 / (df['N'] * df['p'] * (1-df['p'])))
)

# # Statistical significance (Welch's t-test)
# # https://en.wikipedia.org/wiki/Welch%27s_t-test
# df['t'] = ((df['logfc'] - overall_fc)
#            / np.sqrt(df['logfc_err']**2 + overall_fc_err**2))

# Plot
with plt.style.context('tal_paper_spine'):
    plt.figure(figsize=(3.75, 2))
    # plt.scatter(
    #     x=np.arange(df.shape[0]),
    #     y=df['logfc'],
    #     s=df['num_cells'] / 3,
    #     c="#717171"
    # )
    plt.bar(
        x=np.arange(df.shape[0]),
        height=df['logfc'],
        yerr=df['logfc_err'] * 1.98,
        color="#888888",
        capsize=1,
        error_kw={'elinewidth': 1}
    )
    plt.axhline(y=overall_fc, color="#1a8300", linewidth=1.5, alpha=0.8)

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
