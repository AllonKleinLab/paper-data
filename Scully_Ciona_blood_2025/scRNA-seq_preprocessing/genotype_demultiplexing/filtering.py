import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
import scrublet as scr
import os
import sys
import time
from tqdm import tqdm

# Change this path to point to folder containing tal_helper_functions.py
path_to_dropbox = os.environ['PATH_TO_DROPBOX']
sys.path.append(path_to_dropbox + 'klein_lab/resources/helper_functions/')
import scrna_helper_functions as hf

# Scanpy settings
# verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.verbosity = 2

# Plotting settings
plt.style.use('tal_light_spine')

# Set random seed
np.random.seed(seed = 0)

# ============================================================================
# ARGUMENTS/VARIABLES FOR THIS FILTERING

# Parse arguments from command line, if any
args = [
    'version',
    'aligned_genome',   # which genome was data aligned to?
    'max_mito_pct',     # for all libraries, list separated by ","
    'min_num_UMI',      # for all libraries, list separated by ","
    'min_num_genes',    # for all libraries, list separated by ","
    'run_scrublet',
]
arg_dict = hf.get_args([], args)

# Or input arguments here
if arg_dict == {}:
    arg_dict = {
        'version': 'v1',
        'aligned_genome': 'HT2019_KY21_wth_Ens_mito',
        'max_mito_pct': '10,10',
        'min_num_UMI': '800,2000',
        'min_num_genes': '180,180',
        'run_scrublet': 'True' # include iff scrublet shouldn't be run
    }

# Format arguments with commas into a list of ints
for a in arg_dict:
    if ',' in arg_dict[a]:
        arg_dict[a] = [int(x) for x in arg_dict[a].split(',')]
    elif 'True' == arg_dict[a] or 'False' == arg_dict[a]:
        arg_dict[a] = (arg_dict[a] == 'True')

# ============================================================================
# IMPORT COUNT MATRICES

t = time.time()
print('Importing count matrices...')

# Get paths to folders, create folder for saving figures
data_path = f'data/{arg_dict["aligned_genome"]}_alignment/'
out_path = 'filtering_output/'
out_path += arg_dict['version'] + '/'
if not os.path.exists(out_path): os.mkdir(out_path)

adict = {}

libs = [d for d in os.listdir(data_path) if os.path.isdir(data_path + d)]
for lib in libs:
    # Import counts matrix
    adict[lib] = sc.read_10x_h5(data_path + lib + '/raw_feature_bc_matrix.h5')

    # Update genome name to be more useful than "10x_genome"
    # adict[lib].var['genome'] = 'HT2019 assembly KY2021 + Ensembl mito'
    adict[lib].var['genome'] = arg_dict['aligned_genome']

print(f'DONE ({time.time()-t} sec)\n')

# ============================================================================
# DONOR LABELING FROM GENOTYPE DEMULTIPLEXING
# See folder genotype_demultiplexing for how these labels were generated.

t = time.time()
print('Labeling barcodes by donor (only if output files from genotype ' + \
    'demultiplexing are available)...')

donor_ids = {}

for lib in libs:

    # Import donor id labels
    # donor_ids_file = f'genotype_demultiplexing/from_o2/demultiplex_{lib}/donor_ids.tsv'
    
    # !! USING CORRECTED DONOR LABELS
    # see genotype_demultiplexing/demultiplexing_output/summary.key
    donor_ids_file = ('genotype_demultiplexing/demultiplexing_output/' + 
                      f'{lib}_{arg_dict["version"]}.1/donor_reassign.tsv')

    if os.path.exists(donor_ids_file):
        donor_ids[lib] = pd.read_csv(donor_ids_file, sep='\t')

        # Initialize new column in adata.obs
        adict[lib].obs['donor'] = ''

        # Loop through cells
        for i in donor_ids[lib].index:
            bc = donor_ids[lib].loc[i, 'cell']
            this_donor = donor_ids[lib].loc[i, 'donor_id']

            # Relabel as 'animal1' rather than 'donor1'
            # if 'donor' in this_animal:
            #     animal_ind = this_animal.split('donor')[1]
            #     this_animal = 'animal' + animal_ind

            # Add annotation to adata
            if bc in adict[lib].obs.index:
                adict[lib].obs.loc[bc, 'donor'] = this_donor

            # print(f'{i}/{donor_ids[lib].shape[0]} barcodes labeled', end='\r')

        # Add information to adata.uns
        adict[lib].uns['donor_demultiplexing'] = ('cellSNP-lite and vireo ' +
            'used for genotype demultiplexing of individual donor animals.')
        print(f'Donor labels added for {lib}')

del donor_ids
print(f'DONE ({time.time()-t} sec)\n')

# ============================================================================
# HELPFUL GENE LABELING
# 
# Each gene is named according to the HT assembly, KY2021 gene model
# (http://ghost.zool.kyoto-u.ac.jp/download_ht.html) notation as follows: 
# KY21.[chromosome or contig #].[id #]
# (note this leaves out transcript-specific information, since there are often
# >1 transcript version per gene). E.g. KY21.Chr1.1, or KY21.UAContig32.22.
# 
# For each gene I add the following information for each gene as new columns
# in `adata.var`:
#   - Listed human ortholog from KY2021 gene model, if applicable
#   - Whether it is a mitochondrial gene (I will do this in a later section)

t = time.time()
print('Formatting gene labels...')

for lib in libs:

    # ------------------------------------
    # 1. Add human homolog and any gene name to `adata.var`:

    print(f'  {lib}: Add human homologs')

    adict[lib].var['gene_names'] = ''
    adict[lib].var['human_homologs'] = ''

    for iGene in adict[lib].var.index:

        # Remove "KY21:" from the beginning of the id
        if 'KY21:' in iGene:
            iGene_short = iGene.split(':')[1]

        # Add human homologs to adict[lib].var
        if iGene_short in hf.ciona2human:
            # ciona2human outputs a list, so make it a comma separated string
            adict[lib].var['human_homologs'][iGene] = ', '.join(
                hf.ciona2human[iGene_short])

    # ------------------------------------
    # 2. Add `adata.var.id_w_human` to be a shortened HT gene ID with any
    # human homolog listed in parentheses:

    print(f'  {lib}: Add adata.var.id_w_human')

    # Make list of genes formated as desired
    adict[lib].var['id_w_human'] = ''

    for iGene in adict[lib].var.index:

        # Remove "KY21:" from the beginning of the id
        if 'KY21:' in iGene:
            iGene_short = iGene.split(':')[1]

            if adict[lib].var['human_homologs'][iGene] != '':
                adict[lib].var['id_w_human'][iGene] = (f'{iGene_short} (' + 
                    f'{adict[lib].var["human_homologs"][iGene]})')
            else:
                adict[lib].var['id_w_human'][iGene] = iGene_short
        
        # For non-KY21 genes, just use ID
        else:
            adict[lib].var['id_w_human'][iGene] = iGene

    # --------------------------------
    # 3. Make shortened ID the new adata.var.index

    # Make id_w_human new index
    adict[lib].var['short_gene_id'] = [g.split(':')[-1] for g in
                                       adict[lib].var.index.values]
    adict[lib].var.set_index('short_gene_id', inplace=True, drop=True)
    adict[lib].var.head(10)

    # ------------------------------------
    # 3. Idnetifying mitochondrial genes

    print(f'  {lib}: Identify mitochondrial genes')

    # Add a column to adict[lib].var that is a boolean for mitochondrial genes
    if arg_dict['aligned_genome'] == 'HT2019_KY21_with_Ens_mito':
        adict[lib].var['mt'] = adict[lib].var.index.str.contains('ENSCING')
    elif arg_dict['aligned_genome'] == 'HT2019_KY21_with_Griggio_mito':
        adict[lib].var['mt'] = adict[lib].var.index.str.contains('CAD')


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
# SAVE ADATA (each library separately)

for lib in libs:
    # Name library in adata.obs
    adict[lib].obs['library'] = lib
    adict[lib].obs['sample'] = lib

    # Save
    adict[lib].write(out_path + f'filtered_{lib}.h5ad')
