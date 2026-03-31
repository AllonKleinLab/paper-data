import os, sys
import pickle
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

# Change this path to point to folder containing helper functions scripts
path_to_repo_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(os.path.join(path_to_repo_dir, 'helper_functions'))
import scrna_helper_functions as hf
import gillis_style_coexpression_hf as coexpr

# ============================================================================
# OUTPUT PATH

out_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        'coexpr_cons_Hs_vs_Dr_output')
if not os.path.exists(out_path): os.mkdir(out_path)

# ============================================================================
# THESE SPECIES

species1 = hf.Species('Danio', 'rerio')
species2 = hf.Species('Homo', 'sapiens')

# ------------------------------------
# Zfish dataset

print('Importing zebrafish dataset...')
adata1 = sc.read_h5ad(os.path.join(hf.path.external_datasets,
                                   'Dr_zfish_kidney_marrow_Wang.h5ad'))

# Save a copy of raw (unnormalized data)
adata1.layers['raw_unnorm_expression'] = adata1.X

# Total count normalization
# (without normalizing layers['raw_unnorm_expression'])
norm = sc.pp.normalize_total(adata1, exclude_highly_expressed=False,
                             inplace=False, target_sum=10**4)
adata1.X = norm['X']
del norm

# Logarithmize
sc.pp.log1p(adata1, base=10)

# Set `adata.raw` to normalized & logarithmized data
adata1.raw = adata1

# PCA and nearest neighbor graph
sc.pp.pca(adata1)
sc.pp.neighbors(adata1)

# ------------------------------------
# Human dataset

print('Importing human dataset...')
adata2 = sc.read_h5ad(os.path.join(hf.path.external_datasets,
                                   'Hs_Tabula_Sapiens_blood_bone_marrow.h5ad'))

# Save a copy of raw (unnormalized data)
adata2.layers['raw_unnorm_expression'] = adata2.X

# Total count normalization
# (without normalizing layers['raw_unnorm_expression'])
norm = sc.pp.normalize_total(adata2, exclude_highly_expressed=False,
                            inplace=False, target_sum=10**4)
adata2.X = norm['X']
del norm

# Identify highly variable genes using Klein et al. 2015 method
# (assumes normalized counts but NOT log transformed)
hvg_list = hf.find_hvgs(adata2, min_cells=5)
plt.title(f'{len(hvg_list)} highly variable genes')
plt.close()
adata2.var['highly_variable'] = False
for g_idx in hvg_list:
    g = adata2.var.index[g_idx]
    adata2.var.loc[g, 'highly_variable'] = True
print(f'highly variable genes: {len(hvg_list)} passed')

# Logarithmize
sc.pp.log1p(adata2, base=10)

# Before filtering, set `adata2.raw` to normalized & logarithmized data
adata2.raw = adata2

# Scale data
sc.pp.scale(adata2)

# Run PCA on z-scored data
sc.tl.pca(adata2, use_highly_variable=True, n_comps=50)

# Identify k nearest neighbors (Euclidean distance in PCA space)
sc.pp.neighbors(adata2, n_neighbors=15, n_pcs=50, use_rep='X_pca')

# ============================================================================
# Orthology dict

print('\nGetting OrthoFinder orthologs for all gene pairs...')
path_to_dict = os.path.join(path_to_repo_dir, 'helper_functions', 'gene_dictionaries')
with open(os.path.join(path_to_dict, 'hsap2drer.pickle'), 'rb') as fp:
    hs2dr = pickle.load(fp)
    dr2hs = pickle.load(fp)

# ============================================================================
# Get list of genes for coexpression network

# Only consider genes detected in >20 cells
gene_bool1 = np.array([np.sum(adata1.raw.X != 0, axis=0) > 20]).flatten()
gene_bool2 = np.array([np.sum(adata2.raw.X != 0, axis=0) > 20]).flatten()
gene_list1 = adata1.var.index[gene_bool1]
gene_list2 = adata2.var.index[gene_bool2]

# Filter to known orthology gene lists (based on OrthoFinder)
gene_list1 = np.array([g for g in gene_list1 if g in dr2hs])
gene_list2 = np.array([g for g in gene_list2 if g in hs2dr])

# ============================================================================
# Get co-expression ranked networks

# Generate gene coexpression ranked networks
print('Coexpression ranked networks')
print('  species 1')
rank1 = coexpr.coexpr_network(adata1, gene_list1)
print('  species 2')
rank2 = coexpr.coexpr_network(adata2, gene_list2)

# Save networks
np.save(os.path.join(out_path, 'network1_weights'), rank1)
np.save(os.path.join(out_path, 'network2_weights'), rank2)
np.save(os.path.join(out_path, 'network1_gene_list'), gene_list1)
np.save(os.path.join(out_path, 'network2_gene_list'), gene_list2)

# Coexpression conservation scores
print('Coexpression conservation scores')
example_gene_to_plot = [
    ('gfi1aa', 'GFI1'),
    ('gfi1ab', 'GFI1'),
    ('cebpa', 'CEBPA'),
    ('spi1', 'SPI1'),
    ('stat5a', 'STAT5A'),
    ('stat5b', 'STAT5A'),
    ('mki67', 'MKI67'),
]
results = coexpr.calculate_coexpression_conservation_complex_orthology(
    [rank1, rank2],
    gene_list1, gene_list2,
    orthology_1to2=dr2hs,
    orthology_2to1=hs2dr,
    n=20,
    example_gene_to_plot=example_gene_to_plot,
    save_as=os.path.join(out_path, 'example_'),
    num_random_pairs=1, # Set to a large number for Suresh et al.-like score
                        # (see README.md in this directory for explanation)
    orthology_dict=dr2hs,
)

orthologs_coexpr = pd.DataFrame(index=results['gene_pairs'])
orthologs_coexpr[f'{species1.Xx}_gene'] = [
    gene_list1[pair[0]] for pair in results['gene_pairs']
]
orthologs_coexpr[f'{species2.Xx}_gene'] = [
    gene_list2[pair[1]] for pair in results['gene_pairs']
]
orthologs_coexpr['coexpr_conservation'] = results['coexpr_conservation']
orthologs_coexpr['raw_auroc'] = results['auroc_orthologs']

# Save scores
filename = 'coexpr_conservation_all'
orthologs_coexpr.to_csv(os.path.join(out_path, filename + '.csv'))

# Save non-orthologous genes distribution of AUROCs
np.save(os.path.join(out_path, 'non_orthologs_aurocs'), results['auroc_non_orthologs'])
