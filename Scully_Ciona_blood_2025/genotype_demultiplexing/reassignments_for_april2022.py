import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
import os, sys
from scipy.io import mmread

# Change this path to point to folder containing tal_helper_functions.py
path_to_dropbox = os.environ['PATH_TO_DROPBOX']
sys.path.append(path_to_dropbox + 'klein_lab/resources/helper_functions')
import scrna_helper_functions as hf

# Scanpy settings
sc.settings.verbosity = 2

# ============================================================================
# FUNCTIONS

def preprocess(adata, num_pcs, out_path):
    # Logarithmize
    sc.pp.log1p(adata)

    # Get highly variable SNPs
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3,
                                min_disp=0.5)
    sc.pl.highly_variable_genes(adata, show=False)
    plt.savefig(out_path + f'highly_variable_SNPs.png')
    plt.close()

    # PCA
    sc.pp.pca(adata, n_comps=30, use_highly_variable=True)
    sc.pl.pca_variance_ratio(adata, log=True, show=False)
    plt.savefig(out_path + f'pca_var_ratio.pdf')
    plt.close()

    # Nearest neighbors
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=num_pcs)

    # UMAP
    sc.tl.umap(adata)

    # Leiden clustering
    sc.tl.leiden(adata, resolution=0.1)


def enriched_snps(adata_SNP, cluster, cluster_type, verbose=True):
    """
    Get list of genes with score above a threshold
    Want upregulated genes: z-score is > 0 and p-value is significant

    cluster = str, name of cluster
    cluster_type = str, either 'snp_leiden' or 'transcriptome_leiden'
    """
    uns_str = 'rank_genes_groups_' + cluster_type
    snp_markers = pd.DataFrame(adata_SNP.uns[uns_str]['names'])
    snp_scores = pd.DataFrame(adata_SNP.uns[uns_str]['scores'])
    snp_pvals = pd.DataFrame(adata_SNP.uns[uns_str]['pvals_adj'])

    upreg_bool = (snp_pvals[cluster] < 0.01) * (snp_scores[cluster] > 0)

    snp_list = snp_markers[cluster][upreg_bool]
    score_list = snp_scores[cluster][upreg_bool]

    if verbose:
        print(f'Cluster {cluster}: {len(snp_list)} genes pass threshold.')
    return [snp_list, score_list]


# ============================================================================
# SET PATHS AND VARIABLES

lib = 'Cr_blood_mannitol'
version = 'v5.1'

# Get paths to folders, create folder for saving figures
exp_path = 'klein_lab/evolution/ciona/experiments/2022/220428_10x_mannitol/'
demux_path = (path_to_dropbox + exp_path + 'data_processing/' +
              'genotype_demultiplexing/')
adata_path = (path_to_dropbox + exp_path + 'data_processing/' +
              f'2_preprocessing_output/{version}/{lib}/')
snp_data_path = demux_path + f'from_o2/snp_matrices_{lib}/'

# Define output path for figure, etc.
out_path = demux_path + 'demultiplexing_output/'
if not os.path.exists(out_path): os.mkdir(out_path)
out_path += f'{lib}_{version}/'
if not os.path.exists(out_path): os.mkdir(out_path)
out_path_paper_figs = out_path + 'plots_for_paper/'
if not os.path.exists(out_path_paper_figs): os.mkdir(out_path_paper_figs)

# ============================================================================
# IMPORT DATA

print('Importing SNP data')

# Import SNP matrix
AD = mmread(snp_data_path + 'cellSNP.tag.AD.mtx').tocsc()

# Adata input with cell-by-gene expression
adata = sc.read(adata_path + f'adata.h5ad')

# List of barcodes
# ...from cellSNP output
df_CB = pd.read_csv(snp_data_path + 'cellSNP.samples.tsv', sep='\t',
                    header=None)
df_CB = df_CB.rename(columns={0:'cell'})
df_CB['SNP_cell_ID'] = df_CB['cell'].apply(lambda x: x[:-2])
# ...from adata
adata_CB=[x.split('_')[-1]  for x in adata.obs_names ]
df_adata=pd.DataFrame({
    'cell':adata.obs_names,
    'adata_index':np.arange(len(adata_CB))
})
# ...and merged
df_CB['SNP_index'] = np.arange(len(df_CB))
df_adata = df_adata.merge(df_CB,on='cell',how='left')

# ------------------------------------
# Format SNP matrix

# List of SNP indices in adata
sel_idx=df_adata['SNP_index'].values

min_cells_per_SNP = 50
# Get adata bc rows in SNP matrix, keep SNPs present in >min_cells_per_SNP
# cells
var_idx = AD[:,sel_idx].sum(1).A.flatten() > min_cells_per_SNP
SNP_matrix = AD[:,sel_idx][var_idx].transpose()

# Save SNP matrix to adata copy
adata_new = adata[df_adata['adata_index'].values]
adata_new.obsm['X_snp'] = SNP_matrix

# Save as new scanpy AnnData object
adata_SNP = sc.AnnData(SNP_matrix)
adata_SNP.obs.index = adata_new.obs.index
# adata_SNP.obs['donor'] = adata_new.obs['donor']
adata_SNP.obs['transcriptome_leiden'] = adata_new.obs['leiden']
adata_SNP.obs['SNP_count'] = np.sum(adata_SNP.X, axis=1).flatten()\
                                .tolist()[0]

# Get assigned donor labels from Vireo
assigned_donors = pd.read_csv(f'{demux_path}from_o2/demultiplex_{lib}/'
                              + 'donor_ids.tsv', sep='\t', index_col='cell')
adata_SNP.obs['donor_vireo'] = assigned_donors.loc[adata_SNP.obs.index,
                                                   'donor_id']
adata.obs['donor_vireo'] = assigned_donors.loc[adata.obs.index,
                                               'donor_id']


# ============================================================================ 
# Visualize percent vireo doublets per cluster

res = 0.5
sc.tl.leiden(adata, resolution=res, key_added=f'leiden_res{res}')

# UMAPs
with plt.style.context('tal_paper_spine'):
    hf.pl_umap(adata, color='donor_vireo', groups=['doublet'],
               width_scale=2.5, height_scale=2.5, palette='k')
    plt.legend().remove()
    plt.tight_layout()
    plt.savefig(out_path_paper_figs + 'umap_vireo_doublets.png', dpi=300)
    plt.close()

# Bar plot
doublet_frac = {}
for cl in adata.obs[f'leiden_res{res}'].cat.categories:
    adata_sub = adata[adata.obs[f'leiden_res{res}'] == cl]
    cl_size = adata_sub.shape[0]
    doublet_count = np.sum(adata_sub.obs['donor_vireo'] == 'doublet')
    doublet_frac[cl] = doublet_count / cl_size
doublet_frac = pd.Series(doublet_frac)
doublet_frac = doublet_frac.sort_values()[::-1]
with plt.style.context('tal_paper'):
    f = plt.figure(figsize=(2.3, 1.5))
    plt.bar(doublet_frac.index, doublet_frac.values, color='#888888')
    plt.xlabel(f'{lib} Leiden cluster')
    plt.ylabel('Fraction labeled doublets')
    plt.axhline(y=0.5, linestyle='--', color='k', linewidth=0.75)
    plt.xticks(ticks=[], labels=[])
    plt.xlim(-1, len(doublet_frac))
    plt.ylim(0, 1)
    plt.tight_layout()
    plt.savefig(out_path_paper_figs + 'vireo_doublet_frac.pdf')
    plt.close()

# ============================================================================
# PROCESS AND CLUSTER SNP MATRIX
# Similar to pre-processing on scRNA-seq data

print('Processing full SNP data')

num_pcs = 9
preprocess(adata_SNP, num_pcs, out_path)

# Plot
hf.pl_umap(adata_SNP, color=['leiden', 'transcriptome_leiden', 'donor_vireo'],
           ncol=1, width_scale=7, height_scale=4)
plt.savefig(out_path + f'umap_{num_pcs}_pcs.pdf')
plt.close()

pc1 = 1
pc2 = 2
pcs = ','.join([str(pc1), str(pc2)])
f, ax = plt.subplots(3, 1, figsize=[7, 12])
colors = ['leiden', 'transcriptome_leiden', 'donor_vireo']
for i, c in enumerate(colors):
    sc.pl.pca(adata_SNP, color=c, ax=ax[i], show=False,
              components=pcs)
plt.tight_layout()
plt.savefig(out_path + f'pca_{pc1}v{pc2}.pdf')
plt.close()

# Save summary with variables
with open(out_path + 'dimensionality_reduction.csv', 'w') as f:
    f.write(f'Min cells per SNP,{min_cells_per_SNP}\n')
    f.write(f'Number of PCs used,{num_pcs}\n\n')

# ============================================================================
# FOCUS ON CLUSTER SUBSET

# ------------------------------------
# Filter data

# Transcriptomic cluster(s) to keep
clust = '2,11,16,20'

# Filter dataset to include only the cluster or clusters
if ',' in clust:
    bool_overall = np.zeros((adata_SNP.shape[0],))
    for c in clust.split(','):
        bool_cl = (adata_SNP.obs['transcriptome_leiden'] == c)
        bool_overall += bool_cl.astype(int)
    adata_SNP_sub = adata_SNP[bool_overall.astype(bool), :]
else:
    adata_SNP_sub = adata_SNP[adata_SNP.obs['transcriptome_leiden'] == \
        clust, :]

out_path3 = out_path + f'SNP_subset_cluster_{clust}/'
if not os.path.exists(out_path3): os.mkdir(out_path3)

# ------------------------------------
# Re-do processing and clustering of SNP matrix

print(f'Processing SNP data (transcriptomic cluster {clust} only)')

num_pcs = 5
preprocess(adata_SNP_sub, num_pcs, out_path3)

# Plotting
hf.pl_umap(adata_SNP_sub, color=['leiden', 'transcriptome_leiden', 'donor_vireo'],
        ncol=1, width_scale=7, height_scale=4)
plt.savefig(out_path3 + f'umap_{num_pcs}pcs.pdf')
plt.close()

pc1 = 1
pc2 = 2
pcs = ','.join([str(pc1), str(pc2)])
f, ax = plt.subplots(3, 1, figsize=[7, 12])
colors = ['leiden', 'transcriptome_leiden', 'donor_vireo']
for i, c in enumerate(colors):
    sc.pl.pca(adata_SNP_sub, color=c, ax=ax[i], show=False,
              components=pcs)
plt.tight_layout()
plt.savefig(out_path3 + f'pca_{pc1}v{pc2}.pdf')
plt.close()

# Save summary with variables
with open(out_path3 + 'dimensionality_reduction.csv', 'w') as f:
    f.write(f'Number of PCs used,{num_pcs}\n\n')

# ============================================================================
# PROJECT FULL DATA INTO THIS SUBSET'S PC SPACE

out_path4 = out_path3 + '2_pca_project/'
if not os.path.exists(out_path4): os.mkdir(out_path4)

# Project into PC space
adata_SNP_new = adata_SNP.copy()
# adata_SNP_new = adata_SNP_new[adata_SNP.obs['transcriptome_leiden'] == '2']
adata_SNP_new.varm['PCs'] = adata_SNP_sub.varm['PCs']
adata_SNP_new.obsm['X_pca'] = adata_SNP_new.X @ adata_SNP_new.varm['PCs']

# Complete preprocessing with this PC space
sc.pp.neighbors(adata_SNP_new, n_neighbors=10, n_pcs=num_pcs)
sc.tl.umap(adata_SNP_new)#, min_dist=0.5, spread=1.0)
sc.tl.leiden(adata_SNP_new, resolution=0.1)

# Plotting
hf.pl_umap(adata_SNP_new, color=['leiden', 'transcriptome_leiden', 'donor_vireo'],
        ncol=1, width_scale=7, height_scale=4)
plt.savefig(out_path4 + f'umap_projected_pc_space.pdf')
plt.close()

pc1 = 1
pc2 = 5
pcs = ','.join([str(pc1), str(pc2)])
f, ax = plt.subplots(3, 1, figsize=[7, 12])
colors = ['leiden', 'transcriptome_leiden', 'donor_vireo']
for i, c in enumerate(colors):
    sc.pl.pca(adata_SNP_new, color=c, ax=ax[i], show=False,
              components=pcs)
    # ax[i].axvline(x=0.5, linestyle='--', color='k')
plt.tight_layout()
plt.savefig(out_path4 + f'pca_{pc1}v{pc2}.pdf')
plt.close()

# Plot SNPs shown to be enriched in each donor for this cluster
snp_dict = {
    'donor4': '1059',
    'donor2': '25334',
    'donor0': '10120',
    'donor1': '18456',
}
nrow, ncol = hf.get_nrow_ncol(len(snp_dict))
f, ax = plt.subplots(nrow, ncol, figsize=(5*ncol, 5*nrow))
for i, d in enumerate(snp_dict):
    color = snp_dict[d]
    sc.pl.pca(adata_SNP_new, color=color, cmap=hf.cmap_pinks, ax=ax[i],
              show=False, title=f'{d} marker: {color}')
plt.tight_layout()
plt.savefig(out_path4 + f'pca_marker_snps.png', dpi=150)
plt.close()

# ============================================================================
# Donor re-assignement based on PC values

adata_SNP_new.obs['donor_reassign'] = 'unassigned'

# ------------------------------------
# Donor4 thresholds
this_donor = 'donor4'

# Visualize threshold on PC plot
x = 0.3
# m1 = 4; b1 = -2
# m2 = -3; b2 = 2
pc1 = 1; pc2 = 2; pcs = ','.join([str(pc1), str(pc2)])
pc1_data = adata_SNP_new.obsm['X_pca'][:, pc1-1]
pc2_data = adata_SNP_new.obsm['X_pca'][:, pc2-1]

f, ax = plt.subplots(1, 2, figsize=[12, 4])
for iAx in ax:
    sc.pl.pca(adata_SNP_new, color='donor_vireo', ax=iAx, show=False, components=pcs,
            title='donor_vireo (PC space, zoomed in)')#, groups='donor4')
    xmin, xmax = iAx.get_xlim()
    ymin, ymax = iAx.get_ylim()

    iAx.axvline(x=x, linestyle='--', color='k')
    x_arr = np.linspace(iAx.get_xlim()[0], min(x, iAx.get_xlim()[1]), 200)
    iAx.fill_between(x_arr, iAx.get_ylim()[0], iAx.get_ylim()[1],
                    interpolate=True, color='#ceedff', zorder=0)
    # iAx.plot(x_arr, m1*x_arr + b1, 'k--')
    # iAx.plot(x_arr, m2*x_arr + b2, 'k--')
    # iAx.fill_between(x_arr, m1*x_arr + b1, iAx.get_ylim()[1],
    #                 interpolate=True, color='#ceedff', zorder=0)
    # iAx.fill_between(x_arr, m2*x_arr + b2, iAx.get_ylim()[0],
    #                 interpolate=True, color='#ceedff', zorder=0)

ax[0].set_title('donor (PC space)')
ax[-1].text(1.05, 0.98, f'PC{pc1} < {x}', ha="left", va="top",
            transform=ax[-1].transAxes)
# ax[-1].text(1.05, 0.90, f'PC{pc2} > {m1}*PC{pc1} + {b1}', ha="left",
#         va="top", transform=ax[-1].transAxes)
# ax[-1].text(1.05, 0.82, f'PC{pc2} < {m2}*PC{pc1} + {b2}', ha="left",
#         va="top", transform=ax[-1].transAxes)

# Adjust zoom for different axes
ax[0].set_xlim(xmin, xmax)
ax[0].set_ylim(ymin, ymax)
ax[1].set_xlim(-6, 8)
ax[1].set_ylim(-10, 13)

plt.tight_layout()
plt.savefig(out_path4 + f'reassign_{this_donor}_thresholds.pdf')
plt.close()
# plt.show()

# Re-assign donor labels
pc1_data = adata_SNP_new.obsm['X_pca'][:, pc1-1]
pc2_data = adata_SNP_new.obsm['X_pca'][:, pc2-1]
adata_SNP_new.obs['donor_reassign'] = adata_SNP_new.obs['donor_reassign'] \
    .cat.add_categories([this_donor])
adata_SNP_new.obs.loc[
    # row (boolean array)
    (
        (pc1_data < x).astype(int)# *
        # (
        #     ((m1*pc1_data + b1) < pc2_data).astype(int) +
        #     ((m2*pc1_data + b2) > pc2_data).astype(int)
        # )
    ).astype(bool),
    # column name
    'donor_reassign'
] = this_donor

# ------------------------------------
# Not donor4 thresholds
this_donor = 'not_donor4'

# Visualize threshold on PC plot
x = 0.7
pc1 = 1; pc2 = 2; pcs = ','.join([str(pc1), str(pc2)])
pc1_data = adata_SNP_new.obsm['X_pca'][:, pc1-1]
pc2_data = adata_SNP_new.obsm['X_pca'][:, pc2-1]

f, ax = plt.subplots(1, 2, figsize=[12, 4])
for iAx in ax:
    sc.pl.pca(adata_SNP_new, color='donor_vireo', ax=iAx, show=False, components=pcs,
            title='donor_vireo (PC space, zoomed in)')#, groups='donor4')
    xmin, xmax = iAx.get_xlim()
    ymin, ymax = iAx.get_ylim()

    iAx.axvline(x=x, linestyle='--', color='k')
    
    x_arr = np.linspace(max(x, iAx.get_xlim()[0]), iAx.get_xlim()[1], 200)
    iAx.fill_between(x_arr, iAx.get_ylim()[0], iAx.get_ylim()[1],
                    interpolate=True, color='#ceedff', zorder=0)

ax[0].set_title('donor_vireo (PC space)')
ax[-1].text(1.05, 0.98, f'PC{pc1} > {x}', ha="left", va="top",
            transform=ax[-1].transAxes)

# Adjust zoom for different axes
ax[0].set_xlim(xmin, xmax)
ax[0].set_ylim(ymin, ymax)
ax[1].set_xlim(-6, 8)
ax[1].set_ylim(-10, 13)

plt.tight_layout()
plt.savefig(out_path4 + f'reassign_{this_donor}_thresholds.pdf')
plt.close()
# plt.show()

# Re-assign donor labels
pc1_data = adata_SNP_new.obsm['X_pca'][:, pc1-1]
pc2_data = adata_SNP_new.obsm['X_pca'][:, pc2-1]
adata_SNP_new.obs['donor_reassign'] = adata_SNP_new.obs['donor_reassign'] \
    .cat.add_categories([this_donor])
adata_SNP_new.obs.loc[
    # row (boolean array)
    (
        (pc1_data > x).astype(int)
    ).astype(bool),
    # column name
    'donor_reassign'
] = this_donor

# ============================================================================
# Visualize on transcriptomic data

# Rename donor4 to donor1
adata_SNP_new.obs['donor_reassign'] = \
    adata_SNP_new.obs['donor_reassign'].cat.rename_categories({
        'donor4': 'donor1',
        'not_donor4': 'not_donor1'
    })

adata.obs['donor_reassign'] = adata_SNP_new.obs['donor_reassign']

# Original donor labels
title_list = [f'{d} ({np.sum(adata.obs["donor"] == d)} bc)' for d in
              np.unique(adata.obs['donor'])]
hf.pl_umap_separate(adata, 'donor', title=title_list, palette=['k'])
plt.tight_layout()
plt.savefig(out_path4 + f"donors_orig.png", dpi=150)
plt.close()

# New donor labels
title_list = [f'{d} ({np.sum(adata.obs["donor_reassign"] == d)} bc)' for d in
              np.unique(adata.obs['donor_reassign'])]
hf.pl_umap_separate(adata, 'donor_reassign', title=title_list, palette=['k'])
plt.tight_layout()
plt.savefig(out_path4 + f"donors_reassign.png", dpi=150)
plt.close()

colors = ['donor_vireo', 'donor_reassign']
f, ax = plt.subplots(1, len(colors), figsize=(3.5*len(colors), 4))
for i, c in enumerate(colors):
    sc.pl.pca(adata_SNP_new, color=c, show=False, ax=ax[i])
    if (i + 1) != len(color):
        ax[i].legend().remove()
plt.tight_layout()
plt.savefig(out_path4 + f"donors_pca_labels.png", dpi=300)
plt.close()

# ============================================================================
# SAVE DONOR REASSIGNMENT

# Save donor assignment list
save_donors = adata.obs['donor_reassign'].copy()
save_donors = save_donors.reset_index()
save_donors = save_donors.rename(columns = {'index': 'cell',
                                            'donor_reassign': 'donor_id'})

save_donors.to_csv(out_path + 'donor_reassign.tsv', sep='\t')
with open(out_path + 'donor_reassign_info.txt', 'w') as f:
    f.write('HOW DONOR LABELS WERE REASSIGNED\n')
    f.write('1. Worked with the SNP count matrix from cellSNP-lite\n')
    f.write('2. Preprocessed (log transform, highly variable SNPs) and did ' +
            'PCA on SNP matrix, only including barcodes from transcriptomic' +
            f' cluster {clust}*\n')
    f.write('3. Projected the whole data (all barcodes) into the cluster(s)' +
            f' {clust} PC space\n')
    f.write('4. Set thresholds in this PC space to re-label donors**\n') 
    f.write('** See thresholds plotted - files starting with "reassign_" ' +
            f'in {out_path4}\n')

# ============================================================================
# Figures for paper

if not os.path.exists(out_path_paper_figs): os.mkdir(out_path_paper_figs)
with plt.style.context('tal_paper_spine'):
    box_x = [-3, 4]
    box_y = [-6, 6]
    f, ax = plt.subplots(1, 1, figsize=(3, 2.1))
    sc.pl.pca(adata_SNP_new, color='donor_vireo', show=False, ax=ax, s=5)
    plt.title('')
    plt.tight_layout()
    plt.savefig(out_path_paper_figs + 'SNP_pc_vireo.png', dpi=300)
    plt.axvline(x=box_x[0], linewidth=0.75, color='k')
    plt.axvline(x=box_x[1], linewidth=0.75, color='k')
    plt.axhline(y=box_y[0], linewidth=0.75, color='k')
    plt.axhline(y=box_y[1], linewidth=0.75, color='k')
    plt.savefig(out_path_paper_figs + 'SNP_pc_vireo_with_box.png', dpi=300)
    plt.close()

    f, ax = plt.subplots(1, 1, figsize=(3, 2.1))
    sc.pl.pca(adata_SNP_new, color='donor_vireo', show=False, ax=ax)
    plt.title('')
    plt.xlim(box_x[0], box_x[1])
    plt.ylim(box_y[0], box_y[1])
    plt.axvline(x=0.3, linestyle='--', linewidth=0.75, color='k')
    plt.axvline(x=0.7, linestyle='--', linewidth=0.75, color='k')
    plt.tight_layout()
    plt.subplots_adjust(top=0.88)
    plt.savefig(out_path_paper_figs + 'SNP_pc_vireo_zoomed_in.png', dpi=300)
    plt.close()
