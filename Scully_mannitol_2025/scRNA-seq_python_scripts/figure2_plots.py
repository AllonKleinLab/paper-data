import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
import scrublet as scr
import os
import sys
import time
from statistics import median

import helper_functions as hf

# Change this path to point to folder containing tal_helper_functions.py
path_to_dropbox = os.environ['PATH_TO_DROPBOX']

# Set random seed
np.random.seed(seed = 0)

# Set output path
out_path = 'figure2_plots/'
if not os.path.exists(out_path): os.mkdir(out_path)

# ============================================================================
# SET UP FOR SCRNA PLOTS

title_list = ['dCMF-ASW', 'PBS-M']
color_list = ['#147dbe', '#ab004d']#np.array(hf.palette1)[[0, 2]].tolist()

# ------------------------------------
# Arguments for this filtering

# Filtration conditions
arg_dict = {
    'version': 'vx',
    'max_mito_pct': '10,10',
    'min_num_UMI': '800,2000',
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
exp_path = 'klein_lab/evolution/ciona/experiments/2022/220428_10x_mannitol/'
data_path = (path_to_dropbox + exp_path
             + 'data_processing/data/HT2019_KY21_with_Ens_mito_alignment/')

adict = {}

# Import counts matrix
libs = [d for d in os.listdir(data_path) if os.path.isdir(data_path + d)]
for lib in libs:
    adict[lib] = sc.read_10x_h5(data_path + lib + '/raw_feature_bc_matrix.h5')

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

with plt.style.context('tal_light'):
    # Set count thresholds (change in arg_dict above based on plots)
    for i in range(len(libs)):
        adict[libs[i]].uns['min_num_UMI'] = arg_dict['min_num_UMI'][::-1][i]

    ncol = 1
    nrow = len(adict)

    fig = plt.figure(figsize = (ncol * 7, nrow * 2.5))
    for i, lib in enumerate(libs):
        ax0 = plt.subplot(nrow, ncol, i + 1)
        
        min_num_UMI = adict[lib].uns['min_num_UMI']
        
        # Plot histogram
        (freq, bins) = np.histogram(adict[lib].obs['total_counts'],
                                    np.logspace(1, 4.5, 50))
        ax0.bar(bins[:-1], freq*bins[:-1],
                width=0.9*np.diff(bins),
                color=color_list[i], alpha=0.65)
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
        ax0.axvline(x=val_to_plot, linestyle='--', linewidth=2,
                    color='k')
        
        ax0.set_title(lib)
        ax0.grid(False)

        plt.xticks(fontsize=16)
        plt.yticks(fontsize=16)

        xl = np.array(ax0.get_xlim())
        yl = np.array(ax0.get_ylim())
        ax0.set_title(title_list[i], fontsize=20)
    ax0.set_xlabel('Number of unique mRNA transcripts per barcode\n'
                   + '(UMI-filtered mapped reads)',
                   fontsize=18)

    fig.tight_layout()
    plt.savefig(out_path + '2d.pdf')
    plt.close()

# ============================================================================
# 2e - Mitochondrial fraction (only barcodes which are cells)

print('\nPlotting mitochondrial fraction')
print('```````````````````````````````')

with plt.style.context('tal_light'):
    N = {}
    for i in range(len(libs)):
        N[libs[i]] = arg_dict['min_num_UMI'][::-1][i]

    ncol = 1
    nrow = 1#len(adict)

    fig = plt.figure(figsize = (ncol * 7.5, nrow * 5))
    for i, lib in enumerate(libs):

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
        ax.bar(bins[:-1], freq, width=np.diff(bins), color=color_list[i],
            alpha=0.65, label=title_list[i])#color='#999999')

        ax.set_xlim(0,)
        ax.set_xlabel('Mitochondrial fraction (%)', fontsize=18)
        ax.set_ylabel('Number of cells', fontsize=18)
        ax.set_yscale('log')

        plt.xticks(fontsize=16)
        plt.yticks(fontsize=16)

    ax.legend(loc='center right', fontsize=17)
    fig.tight_layout()
    plt.savefig(out_path + '2e.pdf')
    plt.close()

# ============================================================================
# IMPORT PROCESSED DATA

# Import filtered & preprocessed adata
adata_path = (path_to_dropbox + exp_path + 'data_processing/'
            + '2_preprocessing_output/v5.1/all/')
adata = sc.read_h5ad(adata_path + 'adata.h5ad')

# ============================================================================
# 2b - Number of cells recovered

with plt.style.context('tal_light'):
    # % reads in cells (reads which are confidently mapped to transcriptome)
    x = ['\n'.join(t.split(' ')) for t in title_list]
    y = [
        np.sum(adata.obs['sample'] == 'Cr_blood_asw'),
        np.sum(adata.obs['sample'] == 'Cr_blood_mannitol')
    ]

    f, ax = plt.subplots(1, 1, figsize=(3.4, 4))
    for i in range(len(y)):
        plt.bar(i, y[i], color=color_list[i], alpha=0.65, width=0.7)
    plt.axhline(y = 6000, color = 'k', alpha=0.5, linestyle = '--',
                label='Expected cell number')
    plt.xticks(ticks=np.arange(len(x)), labels=x, fontsize=15)
    plt.xlim(plt.xlim()[0] - 0.15, plt.xlim()[1] + 0.15)
    plt.ylabel('Number of cells recovered', fontsize=16)

    plt.tight_layout()
    plt.savefig(out_path + '2b.pdf')
    plt.close()

# ============================================================================
# 2f - UMAP with samples highlighted

print('\nUMAP plots')
print('``````````')

with plt.style.context('tal_light_spine'):
    hf.pl_umap_separate(adata, 'sample', s=30, title=title_list,
                        palette=color_list[::-1], alpha=0.6)
    plt.gcf().set_size_inches(8, 4)
    plt.tight_layout()
    plt.savefig(out_path + '2f.pdf', dpi=300)
    plt.close()

# ============================================================================
# 2g - UMAP density plot

num_neighbors = 200

with plt.style.context('tal_light_spine'):
    f, ax = plt.subplots(1, 1)
    hf.sample_density(
        adata,
        'sample',
        [('Cr_blood_mannitol', 'Cr_blood_asw')],
        num_neighbors=num_neighbors,
        cmap=hf.cmap_pink2blue_r,
        sort_order=False,
        log_offset=1,
    )
    plt.title('Relative Sample Density')

    plt.gcf().set_size_inches(4, 4)
    plt.tight_layout()
    plt.savefig(out_path + '2g.pdf', dpi=300)
    plt.close()
