import os, sys
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import pickle
from scipy.stats import pearsonr
from statsmodels.stats.multitest import fdrcorrection

# Change this path to point to folder containing gene_hf.py
# This imports dictionaries and functions for easily converting gene ids
path_to_repo_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(os.path.join(path_to_repo_dir, 'helper_functions'))
import scrna_helper_functions as hf
import gillis_style_coexpression_hf as coexpr

# Plot style
plt.style.use('tal_paper')

# ============================================================================
# OUTPUT PATH

out_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        'tf_coexpr_Hs_vs_Dr_output')
if not os.path.exists(out_path): os.mkdir(out_path)

# ============================================================================
# THESE SPECIES

species1 = hf.Species('Danio', 'rerio')
species2 = hf.Species('Homo', 'sapiens')

# ------------------------------------
# Zfish dataset

print('Importing zebrafish dataset...')
adata_dr = sc.read_h5ad(os.path.join(hf.path.external_datasets,
                                   'Dr_zfish_kidney_marrow_Wang.h5ad'))

# Save a copy of raw (unnormalized data)
adata_dr.layers['raw_unnorm_expression'] = adata_dr.X

# Total count normalization
# (without normalizing layers['raw_unnorm_expression'])
norm = sc.pp.normalize_total(adata_dr, exclude_highly_expressed=False,
                             inplace=False, target_sum=10**4)
adata_dr.X = norm['X']
del norm

# Logarithmize
sc.pp.log1p(adata_dr, base=10)

# Set `adata.raw` to normalized & logarithmized data
adata_dr.raw = adata_dr

# PCA and nearest neighbor graph
sc.pp.pca(adata_dr)
sc.pp.neighbors(adata_dr)

# ------------------------------------
# Human dataset

print('Importing human dataset...')
adata_hs = sc.read_h5ad(os.path.join(hf.path.external_datasets,
                                     'Hs_Tabula_Sapiens_blood_bone_marrow.h5ad'))

# Save a copy of raw (unnormalized data)
adata_hs.layers['raw_unnorm_expression'] = adata_hs.X

# Total count normalization
# (without normalizing layers['raw_unnorm_expression'])
norm = sc.pp.normalize_total(adata_hs, exclude_highly_expressed=False,
                            inplace=False, target_sum=10**4)
adata_hs.X = norm['X']
del norm

# Identify highly variable genes using Klein et al. 2015 method
# (assumes normalized counts but NOT log transformed)
hvg_list = hf.find_hvgs(adata_hs, min_cells=5)
plt.title(f'{len(hvg_list)} highly variable genes')
plt.close()
adata_hs.var['highly_variable'] = False
for g_idx in hvg_list:
    g = adata_hs.var.index[g_idx]
    adata_hs.var.loc[g, 'highly_variable'] = True
print(f'highly variable genes: {len(hvg_list)} passed')

# Logarithmize
sc.pp.log1p(adata_hs, base=10)

# Before filtering, set `adata_hs.raw` to normalized & logarithmized data
adata_hs.raw = adata_hs

# Scale data
sc.pp.scale(adata_hs)

# Run PCA on z-scored data
sc.tl.pca(adata_hs, use_highly_variable=True, n_comps=50)

# Identify k nearest neighbors (Euclidean distance in PCA space)
sc.pp.neighbors(adata_hs, n_neighbors=15, n_pcs=50, use_rep='X_pca')

# ------------------------------------
# Orthology dict

print('\nGetting OrthoFinder orthologs for all gene pairs...')
path_to_dict = os.path.join(path_to_repo_dir, 'helper_functions', 'gene_dictionaries')
with open(os.path.join(path_to_dict, 'hsap2drer.pickle'), 'rb') as fp:
    hs2dr = pickle.load(fp)
    dr2hs = pickle.load(fp)

# ------------------------------------
# Coexpression conservation

print('Importing coexpression conservation scores...')
path_to_coexpr_data = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                   'coexpr_cons_Hs_vs_Dr_output')
coexpr_conservation = pd.read_csv(
    os.path.join(path_to_coexpr_data, 'coexpr_conservation_all.csv'),
    index_col=0)

# ============================================================================
# Heatmap plotting function
# (Copied and adjusted from SAMap heatmap plotting function)

def plot_together(mtx_hs, mtx_dr, save_as='heatmaps', cmap=None):
    if cmap is None:
        cmap = sns.diverging_palette(230, 10, s=100, l=35, as_cmap=True)

    # Set figure size and aspect ratios
    num_human_genes = mtx_hs.shape[1]
    num_zfish_genes = mtx_dr.shape[1]
    total_gene_rows = max(num_human_genes, num_zfish_genes)
    plot_width = 6.5
    plot_height_adjust = 1.4335
    f, axs = plt.subplots(
        1, 2,
        figsize=(plot_width, 0.08618*total_gene_rows + plot_height_adjust),
        gridspec_kw={'width_ratios': [mtx_hs.shape[0], mtx_dr.shape[0]]}
    )
    f_cbar, ax_cbar = plt.subplots(1, 1, figsize=(0.75, 1.2))
    yticks = [[], []]

    # Get cbar range for all heatmaps
    vmax = max(mtx_hs.max().max(), mtx_dr.max().max())
    vmin = min(mtx_hs.min().min(), mtx_dr.min().min())

    # Plot human expression
    sns.heatmap(
        mtx_hs.T,
        cmap=cmap,
        center=0, vmax=vmax, vmin=vmin,
        xticklabels=True, yticklabels=False,
        ax=axs[0], cbar=False,
    )
    yticks[0] = mtx_hs.columns

    # Plot ciona expression
    sns.heatmap(
        mtx_dr.T,
        cmap=cmap,
        center=0, vmax=vmax, vmin=vmin,
        xticklabels=True, yticklabels=False,
        ax=axs[1],
        cbar_kws={'label': 'Expression z-score'}, cbar_ax=ax_cbar
    )
    yticks[1] = mtx_dr.columns

    # Set gene labels after tight_layout() to get consistent heatmap widths
    f.tight_layout()
    for i, mtx in enumerate([mtx_hs, mtx_dr]):
        axs[i].set_yticks(np.arange(len(yticks[i])) + 0.5)
        axs[i].set_yticklabels(yticks[i])
    f.subplots_adjust(left=0.12, wspace=0.6)

    f.savefig(os.path.join(out_path, f'{save_as}.pdf'))
    plt.close()

    f_cbar.tight_layout()
    f_cbar.savefig(os.path.join(out_path, f'{save_as}_cbar.pdf'))
    plt.close()

# ============================================================================
# Get cluster centroids

centr_hs = hf.cluster_centroids(adata_hs, groupby='cell_type_coarse')
centr_dr = hf.cluster_centroids(adata_dr, groupby='cell_type')

# Adjust human cell type names
human_cell_types = [x.replace('-positive', '+').replace(',', '')\
                    .replace('Kuppfer', 'Kupffer cell')
                    for x in centr_hs.index]
centr_hs.index = human_cell_types

centr_hs_norm = centr_hs / centr_hs.max()
centr_dr_norm = centr_dr / centr_dr.max()
centr_z_hs = (centr_hs - centr_hs.mean()) / centr_hs.std()
centr_z_dr = (centr_dr - centr_dr.mean()) / centr_dr.std()

# ============================================================================

# Number of top coexpression partners in human to look at
n=30

# Import coexpression networks themselves
network1 = np.load(os.path.join(path_to_coexpr_data, 'network1_weights.npy'))
network2 = np.load(os.path.join(path_to_coexpr_data, 'network2_weights.npy'))
gene_list1 = np.load(os.path.join(path_to_coexpr_data, 'network1_gene_list.npy'))
gene_list2 = np.load(os.path.join(path_to_coexpr_data, 'network2_gene_list.npy'))
network1 = pd.DataFrame(network1, index=gene_list1, columns=gene_list1)
network2 = pd.DataFrame(network2, index=gene_list2, columns=gene_list2)

top_n_genes1 = coexpr.get_top_n(np.array(network1), n=n)
top_n_genes2 = coexpr.get_top_n(np.array(network2), n=n)
top_n_genes1 = pd.DataFrame(top_n_genes1, index=gene_list1, columns=gene_list1)
top_n_genes2 = pd.DataFrame(top_n_genes2, index=gene_list2, columns=gene_list2)

del gene_list1, gene_list2

# ------------------------------------
# Coexpression conservation example

goi_list = []

# Genes involved in human/vertebrate hematopoiesis
hs_gene_list = [
    # High coexpression conservation
    'E2F8',

    # Genes related to hematopoiesis
    'SPI1',
    'CEBPA',
    'CEBPB',
    # 'CEBPE',  # no Zebrafish ortholog
    'IRF4',
    'IRF8',
    'JUN',
    'LEF1',
    'GFI1',
    'GATA2',
    'MITF',
    'STAT5A',
    'STAT5B',
    # 'KLF4',   # Low expression
    'RUNX1',
    'TAL1',
    'NFKB1',
    # 'STAT3',  # no Ciona ortholog
    'EGR1',
    'EGR2',
    'EGR3',
    'MAFB',
    'MAF',
]

# For each human TF get homolog in Ciona
for g_hs in hs_gene_list:
    g_dr_list = hs2dr[g_hs]
    for g_dr in g_dr_list:
        if g_dr in centr_dr.columns:
            g_dr_print = g_dr
            goi_list.append([g_hs, g_dr, g_dr_print])

out_path2 = os.path.join(out_path, 'example_coepxression')
if not os.path.exists(out_path2): os.mkdir(out_path2)

score_per_tf = {}
for goi_hs, goi_dr, goi_dr_print in goi_list:

    # Get top n genes in human, sort by decreasing correlation rank
    top_coexpr_partners2 = top_n_genes2.index[top_n_genes2.loc[goi_hs, :]]

    # Get orthologs of these genes in Ciona
    orthologs_sp1 = [goi_dr]
    for g in list(top_coexpr_partners2):#[goi_hs] + list(top_coexpr_partners2):
        if g in hs2dr:
            for g_dr in hs2dr[g]:
                if (g_dr not in orthologs_sp1) and (g_dr in centr_dr.columns):
                    orthologs_sp1.append(g_dr)

    # Ciona matrix to plot
    mtx_dr = centr_z_dr[orthologs_sp1]
    # Remove NaNs
    mtx_dr = mtx_dr.loc[:, mtx_dr.isnull().sum() == 0]
    # Order of genes
    gene_order_dr = np.argsort(mtx_dr.corrwith(mtx_dr.iloc[:, 0]))[::-1]
    mtx_dr = mtx_dr.iloc[:, gene_order_dr]
    # Order of cell types
    mtx_dr = mtx_dr.sort_values(by=mtx_dr.columns[0])[::-1]
    # Gene labeling
    gene_labels = []
    for g_dr in mtx_dr.columns:
        this_str = g_dr.replace('KY21.Chr', 'C')\
            .replace('KY21.UAContig', 'UAC') + ' ('
        g_hs_list = [g_hs for g_hs in dr2hs[g_dr]
                    if (g_hs in top_coexpr_partners2) or (g_hs == goi_hs)]
        this_str += '/'.join(g_hs_list) + ')'
        gene_labels.append(this_str)
    mtx_dr.columns = gene_labels

    # Human matrix to plot
    mtx_hs = centr_z_hs[[goi_hs] + list(top_coexpr_partners2)]
    # Order of genes
    gene_ordered_list_hs = []
    for g_dr in mtx_dr.columns:
        g_hs_str = g_dr.split('(')[1][:-1]
        g_hs_list = g_hs_str.split('/')
        for g_hs in g_hs_list:
            if g_hs not in gene_ordered_list_hs:
                gene_ordered_list_hs.append(g_hs)
    mtx_hs = mtx_hs.loc[:, gene_ordered_list_hs]
    # Order of cell types
    mtx_hs = mtx_hs.sort_values(by=goi_hs)[::-1]

    plot_together(mtx_hs, mtx_dr,
                  save_as=os.path.join('example_coepxression',
                                       f'{goi_hs}_{goi_dr}'),
                  cmap='seismic')
    
    # --------------------------------
    # Convert to score

    # Get correlation within Ciona genes
    corr_values = []
    for g in mtx_dr.columns[1:]:
        corr, pval = pearsonr(mtx_dr.iloc[:, 0], mtx_dr[g])
        corr_values.append([g, corr, pval])
        if pval > 0.05: break
    corr_values = pd.DataFrame(corr_values, columns=['gene', 'corr', 'pval'])
    passed, fdr = fdrcorrection(corr_values['pval'], alpha=0.05)
    corr_values['fdr'] = fdr

    # Get genes above threshold
    threshold_param = 'corr'
    threshold = 0.5
    corr_values = corr_values[corr_values[threshold_param] > threshold]

    # Save score
    score_per_tf[f'{goi_hs} vs. {goi_dr_print}'] = len(corr_values)

# ------------------------------------
# Summary plots

# Bar graph summarizing values for several TFs
tf_list = np.array([tf for tf in score_per_tf])
tf_scores = np.array([score_per_tf[tf] for tf in score_per_tf])
reorder = np.argsort(tf_scores)
f = plt.figure(figsize=(1.75, 3))
plt.barh(tf_list[reorder[:-1]], tf_scores[reorder[:-1]], color='#cc9200', zorder=3)
plt.barh(tf_list[reorder[-1]], tf_scores[reorder[-1]], color='#888888', zorder=3)
plt.ylabel('Cross-species TF pair')
plt.xlabel('Number of\nshared coexpr.\npartners')
plt.xlim(0, 30.5)
plt.ylim(-0.75, len(tf_list)-0.25)
plt.xticks(ticks=[0, 10, 20, 30], labels=[0, 10, 20, 30], fontsize=6.5)
plt.yticks(fontsize=6.5)
plt.grid(color='#d0d0d0', which='major', axis='x', linestyle='-',
         linewidth=0.75, zorder=-1)
plt.tight_layout()
plt.savefig(os.path.join(out_path, 'summary_bar_graph.pdf'))
plt.close()

# Bar graph with AUROCs instead
tf_list = []
tf_scores = []
for tf_hs, tf_dr, tf_dr_print in goi_list:
    foo = coexpr_conservation.loc[(coexpr_conservation['Dr_gene'] == tf_dr)
                                  * (coexpr_conservation['Hs_gene'] == tf_hs)]
    if len(foo) > 0:
        tf_scores.append(foo['raw_auroc'].values[0])
        tf_list.append(f'{tf_hs} vs. {tf_dr_print}')
tf_list = np.array(tf_list)
tf_scores = np.array(tf_scores)
reorder = np.argsort(tf_scores)#[::-1]
f = plt.figure(figsize=(1.75, 3))
plt.barh(tf_list[reorder[:-1]], tf_scores[reorder[:-1]], color='#cc9200', zorder=3)
plt.barh(tf_list[reorder[-1]], tf_scores[reorder[-1]], color='#888888', zorder=3)
plt.axvline(x=0.5, linestyle='--', linewidth=0.75, color='k', zorder=4)
plt.ylabel('Cross-species TF pair')
plt.xlabel('Coexpression\nconservation\nscore (AUROC)')
plt.xlim(0, 1.01)
plt.ylim(-1, len(tf_list))
plt.xticks(ticks=[0, 0.5, 1], labels=[0, 0.5, 1], fontsize=6.5)
plt.yticks(fontsize=6.5)#, rotation=90)
plt.grid(color='#d0d0d0', which='major', axis='x', linestyle='-',
         linewidth=0.75, zorder=-1)
plt.tight_layout()
plt.savefig(os.path.join(out_path, 'summary_bar_graph_AUROCs.pdf'))
plt.close()
