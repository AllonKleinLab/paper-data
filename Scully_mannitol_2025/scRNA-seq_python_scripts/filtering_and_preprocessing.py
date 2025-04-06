import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
import scrublet as scr
import os
import sys
import time
from tqdm import tqdm

import helper_functions as hf

# Plotting settings
plt.style.use('tal_light_spine')

# Set random seed
np.random.seed(seed = 0)

# ============================================================================
# ARGUMENTS/VARIABLES FOR THIS FILTERING

# Or input arguments here
arg_dict = {
    'version': 'v1',
    'aligned_genome': 'HT2019_KY21_with_Ens_mito',
    'max_mito_pct': [10, 10],   # order: dCMF-ASW, PBS-M
    'min_num_UMI': [2000, 800], # order: dCMF-ASW, PBS-M
    'min_num_genes': [0, 0],    # order: dCMF-ASW, PBS-M
    'run_scrublet': True    # include iff scrublet shouldn't be run
}

# ============================================================================
# IMPORT COUNT MATRICES

t = time.time()
print('Importing count matrices...')

# Get paths to folders, create folder for saving figures
data_path = 'raw_unfiltered_data/'
out_path = 'filtering_and_preprocessing_output/'
if not os.path.exists(out_path): os.mkdir(out_path)
out_path += arg_dict['version'] + '/'
if not os.path.exists(out_path): os.mkdir(out_path)

adict = {}

libs = ['_'.join(f.split('_')[:4]) for f in os.listdir(data_path)
        if f[-3:] == '.h5']
for lib in libs:
    # Import counts matrix
    adict[lib] = sc.read_10x_h5(data_path + lib + '_raw_unfiltered_matrix.h5')

    # Update genome name to be more useful than "10x_genome"
    # adict[lib].var['genome'] = 'HT2019 assembly KY2021 + Ensembl mito'
    adict[lib].var['genome'] = arg_dict['aligned_genome']

print(f'DONE ({time.time()-t} sec)\n')

# ============================================================================
# IDENTIFY MITOCHONDRIAL GENES

t = time.time()
print('Identifying mitochondrial genes...')

for lib in libs:
    # Add a column to adict[lib].var that is a boolean for mitochondrial genes
    if arg_dict['aligned_genome'] == 'HT2019_KY21_with_Ens_mito':
        adict[lib].var['mt'] = adict[lib].var.index.str.contains('ENSCING')

print(f'DONE ({time.time()-t} sec)\n')

# ============================================================================
# DATA QUALITY CONTROL

# Calculate summary statistics: number of UMIs per cell/barcode ("count
# depth"), proportion of genes/counts from mtDNA, and total genes per barcode.
for lib in libs:
    sc.pp.calculate_qc_metrics(adict[lib], qc_vars=['mt'], percent_top=None,
                               log1p=False, inplace=True)


# Before starting the filtration: remove all barcodes with 0 total counts (so
# we can visualize this in log space)
with open(out_path + 'filtering_summary.csv', 'w') as f:
    f.write(f'Genome used for alignment: {arg_dict["aligned_genome"]}\n\n')

    f.write('Removing barcodes with 0 UMIs\n')

    f.write('Sample,Filter,# barcodes pre-filter,# barcodes post-filter\n')
    for lib in libs:
        n_orig = adict[lib].shape[0]
        adict[lib] = adict[lib][adict[lib].obs['total_counts'] > 0, :]
        n_pass = adict[lib].shape[0]
        f.write('{},>0,{},{}\n'.format(lib, n_orig, n_pass))

# ------------------------------------
# QC1. Mitochondrial fraction
# ````````````````````````````
# Cells can lyse during RNA-Seq/encapsulation; the mRNA in the cytoplasm might
# leak out, while mtDNA, which is in the mitochondrial  matrix, may remain. So
# high mtDNA proportions may indicate cell death. We'll filter out barcodes
# ~>=20% of mtDNA (roughly--depending on distribution).

t = time.time()
print('Plotting QC thresholds')
print('  Mitochondrial fraction')

# Set mito thresholds (change in arg_dict above based on plots)
for i in range(len(libs)):
    adict[libs[i]].uns['max_mito_pct'] = arg_dict['max_mito_pct'][i]

ncol = min(3, len(adict))
nrow = int(np.ceil(len(adict) / ncol))

fig = plt.figure(figsize = (ncol * 5, nrow * 4))
for i, lib in enumerate(libs):
    ax = plt.subplot(nrow, ncol, i + 1)

    # plot a histogram of the % of counts from mtDNA
    (freq, bins) = np.histogram(adict[lib].obs['pct_counts_mt'], bins=50)
    ax.bar(bins[:-1], freq, width=np.diff(bins), color='#999999')

    # plot the mtDNA threshold as a red, dotted line
    ax.axvline(x = adict[lib].uns['max_mito_pct'], color = '#ff1f5b', 
        linestyle = '--')

    ax.set_ylim(0, 25000)
    ax.set_xlim(0, 100)
    ax.set_xlabel('mtDNA percentage')
    ax.set_ylabel('Number of cells')
    ax.set_title(lib)

    ntot = len(adict[lib].obs['pct_counts_mt'])
    npass = sum(
        adict[lib].obs['pct_counts_mt'] <= adict[lib].uns['max_mito_pct'])
    ax.text(0.98, 0.98, f'Threshold: {adict[lib].uns["max_mito_pct"]}%',
            ha="right", va="top", transform=ax.transAxes)
    ax.text(0.98, 0.85, f'{npass}/{ntot}', ha="right", va="top",
            transform=ax.transAxes)
    ax.text(0.98, 0.72, '{:.1f}%'.format(npass/ntot*100), ha="right",
            va="top", transform=ax.transAxes)

fig.tight_layout()
plt.savefig(out_path + 'qc1_mt_fraction.pdf')
plt.close()

# ------------------------------------
# QC2. UMIs per barcode
# ``````````````````````

print('  UMI/barcode')

# Set count thresholds (change in arg_dict above based on plots)
for i in range(len(libs)):
    adict[libs[i]].uns['min_num_UMI'] = arg_dict['min_num_UMI'][i]

ncol = 3
nrow = len(adict)

fig = plt.figure(figsize = (ncol * 5, nrow * 4))
for i, lib in enumerate(libs):
    ax0 = plt.subplot(nrow, ncol, 3*i + 1)
    ax1 = plt.subplot(nrow, ncol, 3*i + 2)
    ax2 = plt.subplot(nrow, ncol, 3*i + 3)
    
    min_num_UMI = adict[lib].uns['min_num_UMI']
    
    # 1. histogram
    (freq, bins) = np.histogram(adict[lib].obs['total_counts'],
                                np.logspace(1, 4.5, 50))
    
    ax0.bar(bins[:-1], freq, width=np.diff(bins), color='#999999')
    ax0.set_xscale('log')
    ax0.axvline(x = adict[lib].uns['min_num_UMI'], color = '#ff1f5b',
                linestyle = '--')
    
    ax0.set_title(lib)
    ax0.grid(False)
    ax0.set_xlabel('UMIs/barcode')
    ax0.set_ylabel('Frequency')
    
    # if s == 'ci_bl': ax0.set_ylim(0, 1000)
    # elif s == 'ci_dev': ax0.set_ylim(0, 5000)

    ntot = len(adict[lib].obs['total_counts'])
    npass = sum(adict[lib].obs['total_counts'] >= min_num_UMI)
    
    xl = np.array(ax0.get_xlim())
    yl = np.array(ax0.get_ylim())
    
    ax0.text(xl[0] + 10**(np.log10(xl.ptp())*0.05), yl[0] + yl.ptp() * 0.9,
             'Threshold: {} UMI/bc'.format(adict[lib].uns['min_num_UMI']),
             fontsize=14)
    ax0.text(xl[0] + 10**(np.log10(xl.ptp())*0.05), yl[0] + yl.ptp() * 0.8,
             '{}/{}'.format(npass, ntot), fontsize=14)
    ax0.text(xl[0] + 10**(np.log10(xl.ptp())*0.05), yl[0] + yl.ptp() * 0.7,
             '{:.1f}%'.format(npass/ntot*100), fontsize=14)

    # 2. scaled histogram
    ax1.bar(bins[:-1], freq*bins[:-1], width=np.diff(bins), color='#999999')
    ax1.set_xscale('log')
    ax1.axvline(x = adict[lib].uns['min_num_UMI'], color='#ff1f5b',
                linestyle='--')
    
    ax1.set_title(lib)
    ax1.grid(False)
    ax1.set_xlabel('UMIs/barcode')
    ax1.set_ylabel('Frequency * UMIs/barcode')
    
    #if s == 'ci_bl': ax1.set_ylim(0, 1000)
    #elif s == 'ci_dev': ax1.set_ylim(0, 5000)

    ntot = len(adict[lib].obs['total_counts'])
    npass = sum(adict[lib].obs['total_counts'] >= min_num_UMI)
    
    xl = np.array(ax1.get_xlim())
    yl = np.array(ax1.get_ylim())
    
    ax1.text(xl[0] + 10**(np.log10(xl.ptp())*0.05), yl[0] + yl.ptp() * 0.9,
             'Threshold: {} UMI/bc'.format(adict[lib].uns['min_num_UMI']),
             fontsize=14)
    ax1.text(xl[0] + 10**(np.log10(xl.ptp())*0.05), yl[0] + yl.ptp() * 0.8,
             '{}/{}'.format(npass, ntot), fontsize=14)
    ax1.text(xl[0] + 10**(np.log10(xl.ptp())*0.05), yl[0] + yl.ptp() * 0.7,
             '{:.1f}%'.format(npass/ntot*100), fontsize=14)

    # 3. elbow plot
    countranks = list(adict[lib].obs['total_counts'].sort_values(
        ascending=False))
    ax2.plot(countranks, color = '#888888')
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax2.axhline(y = adict[lib].uns['min_num_UMI'], color='#ff1f5b',
                linestyle='--')

    #ax2.set_xlim(0, 20000)
    #ax2.set_ylim(10, 10000)

    ax2.set_xlabel('Barcode rank')
    ax2.set_ylabel('UMIs per barcode')
    ax2.set_title(lib)
    
fig.tight_layout()
plt.savefig(out_path + 'qc2_umi_per_bc.pdf')
plt.close()

# ------------------------------------
# QC3. Total number of genes
# ```````````````````````````
# Does a particular barcode have unusually few genes associated with it? Then
# the cell may have lysed, or the replication step went wrong, or something
# similar.

print('  Genes/barcode')

# Set gene thresholds (change in arg_dict above based on plots)
for i in range(len(libs)):
    adict[libs[i]].uns['min_num_genes'] = arg_dict['min_num_genes'][i]

ncol = min(3, len(adict))
nrow = int(np.ceil(len(adict) / ncol))

fig = plt.figure(figsize = (ncol * 5, nrow * 4))
for i, lib in enumerate(libs):
    ax = plt.subplot(nrow, ncol, i + 1)

    (freq, bins) = np.histogram(adict[lib].obs['n_genes_by_counts'], bins=50)
    ax.bar(bins[:-1], freq, width=np.diff(bins), color='#999999')
    ax.axvline(adict[lib].uns['min_num_genes'], color='#ff1f5b',
               linestyle='--')

    #ax.set_xlim(0, 5000)
    ax.set_ylim(0, 50000)

    ntot = len(adict[lib].obs['n_genes_by_counts'])
    npass = sum(adict[lib].obs['n_genes_by_counts'] >= adict[lib].uns[
        'min_num_genes'])

    ax.text(0.98, 0.98, (f'Threshold: {adict[lib].uns["min_num_genes"]} ' +
            'genes/bc'), ha="right", va="top", transform=ax.transAxes)
    ax.text(0.98, 0.85, f'{npass}/{ntot}', ha="right", va="top",
            transform=ax.transAxes)
    ax.text(0.98, 0.72, '{:.1f}%'.format(npass/ntot*100), ha="right",
            va="top", transform=ax.transAxes)

    ax.set_xlabel('Number of genes')
    ax.set_ylabel('Barcodes')
    ax.set_title(lib)

fig.tight_layout()
plt.savefig(out_path + 'qc3_gene_per_bc.pdf')
plt.close()

#-------------------------------------
# QC4. Plot the thresholds together
# ``````````````````````````````````

print('  All together')

ncol = 4
nrow = len(adict)

fig = plt.figure(figsize = (ncol * 5, nrow * 5))
for i, lib in enumerate(libs):
    ax1 = plt.subplot(nrow, ncol, 4*i + 1)
    ax2 = plt.subplot(nrow, ncol, 4*i + 2)
    ax3 = plt.subplot(nrow, ncol, 4*i + 3)
    ax4 = plt.subplot(nrow, ncol, 4*i + 4)

    ntot = adict[lib].shape[0]
    pass1 = adict[lib].obs['pct_counts_mt'] <= adict[lib].uns['max_mito_pct']
    pass2 = adict[lib].obs['total_counts'] >= adict[lib].uns['min_num_UMI']
    pass3 = adict[lib].obs['n_genes_by_counts'] >= adict[lib].uns[
        'min_num_genes']
    npass = sum(pass1 & pass2 & pass3)

    # --------------------------------
    # 1. UMIs/bc vs. mito fraction

    # plot a histogram of the % of counts from mtDNA
    ax1.scatter(
        adict[lib].obs['total_counts'],
        adict[lib].obs['pct_counts_mt'],
        color='#000000', alpha=0.05, s=1)
    ax1.set_xscale('log')
    
    ax1.axvline(x = adict[lib].uns['min_num_UMI'], color='#ff1f5b',
                linestyle='--')
    ax1.axhline(y = adict[lib].uns['max_mito_pct'], color='#ff1f5b',
                linestyle='--')

    #ax.set_ylim(0, 200)
    ax1.set_xlabel('UMIs/barcode')
    ax1.set_ylabel('mtDNA percentage')
    ax1.set_title(lib)

    # --------------------------------
    # 2. UMIs/bc vs. mito fraction, add coloring based on gene/bc

    # Categorial shading of genes/barcode (above-below threshold)
    gene_mask = adict[lib].obs['n_genes_by_counts'] >= adict[lib].uns[
        'min_num_genes']
    points_passing = adict[lib][gene_mask]
    points_failing = adict[lib][np.invert(gene_mask)]
    
    # plot a histogram of the % of counts from mtDNA
    ax2.scatter(
        points_passing.obs['total_counts'],
        points_passing.obs['pct_counts_mt'],
        color='#000000', alpha=0.8, s=1)#, label='passing (genes/bc)')
    ax2.scatter(
        points_failing.obs['total_counts'],
        points_failing.obs['pct_counts_mt'],
        color='#d46a00', alpha=0.8, s=1, label='failing (genes/bc)')
    ax2.set_xscale('log')
    
    ax2.axvline(x = adict[lib].uns['min_num_UMI'], color='#ff1f5b',
                linestyle='--')
    ax2.axhline(y = adict[lib].uns['max_mito_pct'], color='#ff1f5b',
                linestyle='--')

    #ax.set_ylim(0, 200)
    ax2.set_xlabel('UMIs/barcode')
    ax2.set_ylabel('mtDNA percentage')
    ax2.set_title(lib)
    ax2.legend(loc='upper left')

    # --------------------------------
    # 3. UMIs/bc vs. gene/bc

    ax3.scatter(
        adict[lib].obs['total_counts'],
        adict[lib].obs['n_genes_by_counts'],
        color='#000000', alpha=0.05, s=2)
    ax3.set_xscale('log')
    
    ax3.axvline(x = adict[lib].uns['min_num_UMI'], color='#ff1f5b',
                linestyle='--')
    ax3.axhline(y = adict[lib].uns['min_num_genes'], color='#ff1f5b',
                linestyle='--')
    
    ax3.set_xlabel('UMIs/barcode')
    ax3.set_ylabel('Number of genes')
    ax3.set_title(lib)
    
    # --------------------------------
    # 4. UMIs/bc vs. gene/bc, add coloring based on mito fraction
    
    # Categorial shading of mitochondria fractions (above-below threshold)
    mito_mask = adict[lib].obs['pct_counts_mt'] <= adict[lib].uns[
        'max_mito_pct']
    points_passing = adict[lib][mito_mask]
    points_failing = adict[lib][np.invert(mito_mask)]

    ax4.scatter(
        points_passing.obs['total_counts'],
        points_passing.obs['n_genes_by_counts'],
        color='#000000', alpha=0.8, s=2)#, label='passing (% mito)')
    ax4.scatter(
        points_failing.obs['total_counts'],
        points_failing.obs['n_genes_by_counts'],
        color='#d46a00', alpha=0.8, s=2, label='failing (% mito)')
    ax4.set_xscale('log')
    
    ax4.axvline(x = adict[lib].uns['min_num_UMI'], color='#ff1f5b',
                linestyle='--')
    ax4.axhline(y = adict[lib].uns['min_num_genes'], color='#ff1f5b',
                linestyle='--')
    
    ax4.set_xlabel('UMIs/barcode')
    ax4.set_ylabel('Number of genes')
    ax4.set_title(lib)
    ax4.legend(loc='upper left')

fig.tight_layout()
plt.savefig(out_path + 'qc4_all.png', dpi=300)
plt.close()

print(f'DONE ({time.time()-t} sec)\n')

#-------------------------------------
# Now actually filter
# ````````````````````

print('Filtering')

with open(out_path + 'filtering_summary.csv', 'a') as f:

    # QC1. UMIs / barcode
    f.write('\nUMI PER BARCODE\n')
    f.write('Sample,Filter,# barcodes pre-filter,# barcodes post-filter\n')
    for lib in libs:
        n_orig = adict[lib].shape[0]
        adict[lib] = (adict[lib][adict[lib].obs['total_counts'] >=
                      adict[lib].uns['min_num_UMI']])
        n_pass = adict[lib].shape[0]
        f.write('{},>={},{},{}\n'.format(lib, adict[lib].uns['min_num_UMI'],
                n_orig, n_pass))

    # QC2. Mitochondrial fraction
    f.write('\nMITOCHONDRIAL FRACTION\n')
    f.write('Sample,Filter,# barcodes pre-filter,# barcodes post-filter\n')
    for lib in libs:
        n_orig = adict[lib].shape[0]
        adict[lib] = (adict[lib][adict[lib].obs['pct_counts_mt'] <=
                      adict[lib].uns['max_mito_pct']])
        n_pass = adict[lib].shape[0]
        f.write('{},<={},{},{}\n'.format(lib, adict[lib].uns['max_mito_pct'],
                n_orig, n_pass))

    # QC3. Gene counts
    f.write('\nGENE COUNTS\n')
    f.write('Sample,Filter,# barcodes pre-filter,# barcodes post-filter\n')
    for lib in libs:
        n_orig = adict[lib].shape[0]
        adict[lib] = (adict[lib][adict[lib].obs['n_genes_by_counts'] >=
                      adict[lib].uns['min_num_genes']])
        n_pass = adict[lib].shape[0]
        f.write('{},>={},{},{}\n'.format(lib, adict[lib].uns['min_num_genes'],
                n_orig, n_pass))

print('DONE')

# ============================================================================
# DOUBLET DETECTION

# ------------------------------------
# Detect doublets

if 'run_scrublet' in arg_dict and arg_dict['run_scrublet']:

    print('\nSCRUBLET')

    doublet_scores = {}
    predicted_doublets = {}

    for lib in libs:
        print('\nSAMPLE:', lib)

        # Run scrublet
        scrub = scr.Scrublet(adict[lib].X, expected_doublet_rate = 0.06)
        doublet_scores[lib], predicted_doublets[lib] = scrub.scrub_doublets(
            min_counts = 2, min_cells = 3)

        # Histogram of doublet scores
        scrub.plot_histogram()
        plt.savefig(out_path + f'scrublet_{lib}_hist.png', dpi=300)

        # UMAP of doublet scores
        scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10,
                            min_dist=0.5))
        scrub.plot_embedding('UMAP', order_points = True)
        plt.savefig(out_path + f'scrublet_{lib}_umap.png', dpi=300)

        # Save doublet scores to adata
        adict[lib].obs['doublet_scores'] = doublet_scores[lib]

    # ------------------------------------
    # Filter doublets

    print('\nFiltering out doublets')

    with open(out_path + 'filtering_summary.csv', 'a') as f:

        f.write('\nDOUBLET FILTERING WITH SCRUBLET\n')

        f.write('Sample,Filter,# barcodes pre-filter,# barcodes post-filter')
        f.write('\n')
        for lib in libs:
            n_orig = adict[lib].shape[0]
            adict[lib] = adict[lib][np.invert(predicted_doublets[lib]), :]
            n_pass = adict[lib].shape[0]
            f.write('{},doublets,{},{}\n'.format(lib, n_orig, n_pass))
            adict[lib].uns['scrublet'] = ('Ran scrublet, filtered out cells,'
                + ' doublet scores saved in obs.doublet_scores')

    print('DONE')

# ============================================================================
# MERGE

libs = [l for l in adict]

# Merge
adata = sc.AnnData.concatenate(*[adict[lib] for lib in libs], 
                               batch_categories=libs,
                               batch_key='library')

# Label samples as appropriate
adata.obs['sample'] = adata.obs['library']

# Transfer adict[lib].uns to adata
adata.uns = {}
for lib in libs:
    adata.uns.update({x+'-'+lib:adict[lib].uns[x] for x in adict[lib].uns})

# Set any NaN in donor column to 'unassigned'
if 'donor' in adata.obs.columns:
    adata.obs['donor'] = adata.obs['donor'].fillna('unassigned')

# ============================================================================
# NORMALIZATION & GENE FILTRATION

# Save a copy of raw (unnormalized data)
adata.layers['raw_unnorm_expression'] = adata.X

print('\nNormalization and gene filtration:')

# Total count normalization (without normalizing layers['raw_unnorm_expression'])
norm = sc.pp.normalize_total(adata, exclude_highly_expressed=False,
                             inplace=False)
adata.X = norm['X']
del norm

# Logarithmize
sc.pp.log1p(adata)

# Before filtering, set `adata.raw` to normalized & logarithmized data
adata.raw = adata

# Remove lowly expressed genes
sc.pp.filter_genes(adata, min_cells=5)

# Identify highly variable genes
sc.pp.highly_variable_genes(adata)
print('{} genes passing filter'.format(np.sum(adata.var['highly_variable'])))
sc.pl.highly_variable_genes(adata, show=False)
plt.savefig(out_path + "highly_variable_genes.png", dpi=300)
plt.close()

# Scale data
sc.pp.scale(adata)

# ============================================================================
# DIMENSIONALITY REDUCTION

print('\nDimensionality reduction:')

# Run PCA on z-scored data
sc.tl.pca(adata, use_highly_variable=True, n_comps=50)
sc.pl.pca(adata, color='total_counts', show=False,
          cmap=hf.cmap(['#f7eef0', '#ea707e', '#8a0015']))
plt.savefig(out_path + "pca_pc1_vs_pc2.png", dpi=300)
plt.close()
sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True, show=False)
plt.savefig(out_path + "pca_var_ratio.png", dpi=300)
plt.close()

# Identify k nearest neighbors (Euclidean distance in PCA space)
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=50, use_rep='X_pca')

# Get UMAP embedding
sc.tl.umap(adata, random_state=0)

# Cluster data (leiden)
# sc.tl.leiden(adata)

# ============================================================================
# PLOTTING

print('Plotting')

# Leiden clusters
hf.pl_umap_separate(adata, color='sample')
plt.gcf().set_size_inches(8, 4)
plt.tight_layout()
plt.savefig(out_path + "umap.png", dpi=300)
plt.close()

# ============================================================================
# SAVE ADATA

adata.write(out_path + f"adata.h5ad")
