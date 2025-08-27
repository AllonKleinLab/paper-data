import os, sys
import pickle
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm
from Bio.KEGG import REST
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import fdrcorrection

# Change this path to point to folder containing helper functions scripts
path_to_repo_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(os.path.join(path_to_repo_dir, 'helper_functions'))
import scrna_helper_functions as hf

# ============================================================================
# OUTPUT PATH

out_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        'kegg_pathway_fishers_output')
if not os.path.exists(out_path): os.mkdir(out_path)

# ============================================================================
# THESE SPECIES

species1 = hf.Species('Ciona', 'robusta')
species2 = hf.Species('Homo', 'sapiens')

# ------------------------------------
# Ciona dataset

# Import counts matrix
adata_cr = sc.read_h5ad(hf.path.Crob_adata_file)

# ------------------------------------
# Human dataset

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

# Logarithmize
sc.pp.log1p(adata_hs, base=10)

# Set `adata.raw` to normalized & logarithmized data
adata_hs.raw = adata_hs

# ============================================================================
# DEGs for each species

try:
    with open(os.path.join(out_path, 'degs_dict.p'), 'rb') as f:
        degs_cr = pickle.load(f)
        degs_hs = pickle.load(f)

except:
    print(f'DEGs for {species1.name}')
    sc.tl.rank_genes_groups(adata_cr, groupby='cell_type', method='wilcoxon')
    degs_cr = {}
    for cell_type in tqdm(adata_cr.obs['cell_type'].cat.categories):
        df = sc.get.rank_genes_groups_df(adata_cr, cell_type)
        these_genes = set(df['names'][(df['logfoldchanges'] > 2)
                                    * (df['pvals_adj'] < 0.05)])
        degs_cr[cell_type] = these_genes

    print(f'DEGs for {species2.name}')
    sc.tl.rank_genes_groups(adata_hs, groupby='cell_type_coarse',
                            method='wilcoxon')
    degs_hs = {}
    for cell_type in tqdm(adata_hs.obs['cell_type_coarse'].cat.categories):
        df = sc.get.rank_genes_groups_df(adata_hs, cell_type)
        these_genes = set(df['names'][(df['logfoldchanges'] > 2)
                                    * (df['pvals_adj'] < 0.05)])
        degs_hs[cell_type] = these_genes

    with open(os.path.join(out_path, 'degs_dict.p'), 'wb') as f:
        pickle.dump(degs_cr, f)
        pickle.dump(degs_hs, f)

# ============================================================================
# IMPORT KEGG PATHWAYS

try:
    with open(os.path.join(out_path, 'kegg_pathways_dict.p'), 'rb') as f:
        pathway_names = pickle.load(f)
        pathway_dict = pickle.load(f)

except:
    # Get pathway list for Homo sapiens
    response = REST.kegg_list("pathway", "hsa")
    pathways = response.read()
    pathways = pd.DataFrame(
        [l.split('\t')[-1].replace(' - Homo sapiens (human)', '')
         for l in pathways.split('\n')],
        index=[l.split('\t')[0] for l in pathways.split('\n')],
        columns=['name'])
    pathways = pathways.loc[[p for p in pathways.index if p!='']]

    pathway_dict = {}
    for pathway_id in tqdm(pathways.index):
        # Query this pathway
        response = REST.kegg_get(pathway_id)
        details = response.read()

        # Ignore if Human Diseases or Drug Development
        pathway_class = details.split('CLASS')[-1].strip().split(';')[0]
        if pathway_class not in ['Human Diseases', 'Drug Development']:
            # Extract gene list from output
            pathway_dict[pathway_id] = []
            lines = details.split('\n')
            these_lines_are_genes = False
            count = 0
            while lines[count] != '///':

                # Determine if these lines contain gene symbols
                if lines[count].startswith('GENE'):
                    these_lines_are_genes = True
                elif lines[count].startswith('COMPOUND'):
                    these_lines_are_genes = False
                
                # Save gene symbols
                if these_lines_are_genes:
                    if len(lines[count].split(';')) > 1:
                        gene = lines[count].split(';')[0].split(' ')[-1]
                        if (gene in adata_hs.var.index) or (gene in hf.human2ciona):
                            pathway_dict[pathway_id].append(gene)
                        else:
                            print(f'{pathway_id}: {gene} not in '
                                + 'adata_hs.var.index or hf.human2ciona')
                
                # Read next line
                count += 1

    pathway_names = {x: pathways.loc[x].iloc[0] for x in pathways.index}
    del pathways
    with open(os.path.join(out_path, 'kegg_pathways_dict.p'), 'wb') as f:
        pickle.dump(pathway_names, f)
        pickle.dump(pathway_dict, f)

# ============================================================================
# PATHWAYS WITH ENOUGH CIONA HOMOLOGS TO INCLUDE

pathway_list = []
pathway_dict_ciona = {}
for term in tqdm(pathway_dict):
    pathway_dict_ciona[term] = []
    this_gene_set = pathway_dict[term]
    
    # Track how many human genes have orthologs
    count_human_genes_with_orthologs = 0

    for human_g in this_gene_set:
        # If there is >1 ortholog...
        if human_g in hf.human2ciona:
            # ...increment count of human genes with orthologs
            count_human_genes_with_orthologs += 1

            # ...and add each ciona ortholog to the gene set (without repeats)
            ciona_g_list = hf.human2ciona[human_g]
            for ciona_g in ciona_g_list:
                if ciona_g not in pathway_dict_ciona[term]:
                    pathway_dict_ciona[term].append(ciona_g)
    
    # If >=10 human genes have orthologs
    if count_human_genes_with_orthologs >= 10:
        # AND if >=10 ciona genes in final list
        if len(pathway_dict_ciona[term]) >= 10:
            pathway_list.append(term)


# ============================================================================
# HYPERGEOMETRIC TEST / FISHER'S EXACT TEST

species_data_list = [
    ('cr', adata_cr, degs_cr, pathway_dict_ciona),
    ('hs', adata_hs, degs_hs, pathway_dict)
]

for species, adata, degs, sp_pathway_dict in species_data_list:
    out_path2 = os.path.join(out_path, species)
    if not os.path.exists(out_path2): os.mkdir(out_path2)

    try:
        fisher_results = pd.read_csv(os.path.join(out_path2, 'KEGG_PATHWAY_fisher_results.tsv'),
                                        sep='\t', index_col=0)
        fisher_results['Cell type'] = fisher_results['Cell type'].astype('category')
        fisher_results['KEGG pathway ID'] = fisher_results['KEGG pathway ID'].astype('category')
        fisher_results['KEGG pathway name'] = fisher_results['KEGG pathway name'].astype('category')

    except:
        # Cell types to look at
        cell_type_list = [c for c in degs]

        # Define the world (all genes which could be in a set)
        gene_bool = np.array(np.sum(adata.raw.X != 0, axis=0) > 20).flatten()
        # gene_bool = adata.var['highly_variable']
        full_set = adata.var.index[gene_bool]

        # Initialize pandas dataframe to save info
        fisher_results = {
            'Cell type': [],
            'KEGG pathway ID': [],
            'KEGG pathway name': [],
            'Enrichment': [],
            'p-value': [],
            'num_genes_overlap': [],
            'num_genes_in_gene_set': [],
        }

        # Loop through every pair of clusters / cell types
        count = 0
        for i, cell_type in tqdm(enumerate(cell_type_list)):
            # print(f'{cell_type} ({i}/{len(cell_type_list)})')
            for kegg_term in pathway_list:
                kegg_name = pathway_names[kegg_term]

                # Get lists of genes
                genes_this_set = set(sp_pathway_dict[kegg_term])

                # Get list of genes for this cell type
                genes_c = degs[cell_type]

                # Filter so gene sets not in "world" are ignored
                genes_this_set = {g for g in genes_this_set if g in full_set}
                genes_c = {g for g in genes_c if g in full_set}

                # Create contingency table
                # Yes Ciona cell type and yes GO term
                a = len(genes_c.intersection(genes_this_set))
                # Yes Ciona cell type and not GO term
                b = len(genes_c.difference(genes_this_set))
                # Not Ciona cell type and yes GO term
                c = len(genes_this_set.difference(genes_c))
                # Not Ciona cell type and not GO term
                # a + b + c + d = len(full_set)
                d = len(full_set) - a - b - c

                contingency_table = np.array([[a, b], [c, d]])

                # Perform Fisher's exact test
                [enrichment, pval] = fisher_exact(contingency_table,
                                                alternative='greater')
                fisher_results['Cell type'].append(cell_type)
                fisher_results['KEGG pathway ID'].append(kegg_term)
                fisher_results['KEGG pathway name'].append(kegg_name)
                fisher_results['Enrichment'].append(enrichment)
                fisher_results['p-value'].append(pval)
                fisher_results['num_genes_overlap'].append(a)
                fisher_results['num_genes_in_gene_set'].append(a + c)

                count += 1

        fisher_results = pd.DataFrame(fisher_results)
        fisher_results['Cell type'] = fisher_results['Cell type'].astype('category')
        fisher_results['KEGG pathway ID'] = fisher_results['KEGG pathway ID'].astype('category')
        fisher_results['KEGG pathway name'] = fisher_results['KEGG pathway name'].astype('category')

        # --------------------------------
        # MULTIPLE HYPOTHESIS CORRECTION

        # Run FDR Benjamini-Hochberg correction
        [rejected, pval_fdr] = fdrcorrection(fisher_results['p-value'],
                                            alpha=0.05)

        fisher_results['null rejected'] = rejected
        fisher_results['FDR'] = pval_fdr

        # Save results
        fisher_results[fisher_results['null rejected']].to_csv(
            os.path.join(out_path2, f'KEGG_PATHWAY_fisher_results.tsv'), sep='\t'
        )

# # ============================================================================
# # OVERVIEW MATRIX

# # Ciona results
# fisher_cr = pd.read_csv(
#     os.path.join(out_path, 'cr', 'KEGG_PATHWAY_fisher_results.tsv'),
#     sep='\t',
#     index_col=0
# )
# fisher_cr['Cell type'] = fisher_cr['Cell type'].astype('category')
# fisher_cr['KEGG pathway ID'] = fisher_cr['KEGG pathway ID'].astype('category')
# fisher_cr['KEGG pathway name'] = fisher_cr['KEGG pathway name'].astype('category')

# # Human results
# fisher_hs = pd.read_csv(
#     os.path.join(out_path, 'hs', 'KEGG_PATHWAY_fisher_results.tsv'),
#     sep='\t',
#     index_col=0
# )
# fisher_hs['Cell type'] = fisher_hs['Cell type'].astype('category')
# fisher_hs['KEGG pathway ID'] = fisher_hs['KEGG pathway ID'].astype('category')
# fisher_hs['KEGG pathway name'] = fisher_hs['KEGG pathway name'].astype('category')

# # Get list of KEGG pathway gene sets to plot
# sets_to_plot = set({})
# for cell_type in fisher_cr['Cell type'].cat.categories:
#     these_terms = fisher_cr.loc[fisher_cr['Cell type'] == cell_type]
#     these_terms = these_terms.sort_values(by='Enrichment', ascending=False)
#     sets_to_plot = sets_to_plot.union(set(these_terms['KEGG pathway ID']))

# # Matrices to plot
# mtx_cr = pd.DataFrame(0, index=fisher_cr['Cell type'].cat.categories,
#                       columns=list(sets_to_plot))
# mtx_hs = pd.DataFrame(0, index=fisher_hs['Cell type'].cat.categories,
#                       columns=list(sets_to_plot))

# # Populate with enrichment if statistically significant
# for kegg_id in sets_to_plot:
#     cr_hits_df = fisher_cr[fisher_cr['KEGG pathway ID'] == kegg_id]
#     cr_hits_df = cr_hits_df.set_index('Cell type')
#     mtx_cr.loc[cr_hits_df.index, kegg_id] = cr_hits_df['Enrichment']

#     hs_hits_df = fisher_hs[fisher_hs['KEGG pathway ID'] == kegg_id]
#     hs_hits_df = hs_hits_df.set_index('Cell type')
#     mtx_hs.loc[hs_hits_df.index, kegg_id] = hs_hits_df['Enrichment']

# # Order rows and columns
# # ...Ciona cell types
# # cell_type_order_idx = np.argsort(np.sum(mtx_cr > 0, axis=1))[::-1]
# cell_type_order_idx = hf.hierarchical_clustering(mtx_cr)['leaves'][::-1]
# cell_type_order_cr = mtx_cr.index[cell_type_order_idx]
# # ...human cell types
# # cell_type_order_idx = np.argsort(np.sum(mtx_hs > 0, axis=1))[::-1]
# cell_type_order_idx = hf.hierarchical_clustering(
#     mtx_hs.loc[mtx_hs.sum(axis=1)!=0])['leaves'][::-1]
# cell_type_order_idx += list(np.where(mtx_hs.sum(axis=1)==0)[0])
# cell_type_order_hs = mtx_hs.index[cell_type_order_idx]
# # ...KEGG pathways based on Ciona dataset
# kegg_max_match_to_cr = mtx_cr.idxmax(axis=0)
# kegg_id_order = []
# for cr_cell in cell_type_order_cr:
#     these_kegg = kegg_max_match_to_cr.index[
#         np.where(kegg_max_match_to_cr == cr_cell)[0]]
#     marginals = mtx_cr[these_kegg].max(axis=0)
#     these_kegg_ordered = these_kegg[np.argsort(marginals)[::-1]]
#     kegg_id_order += list(these_kegg_ordered)
# # Actually reorder matrices
# mtx_cr = mtx_cr.loc[cell_type_order_cr, kegg_id_order]
# mtx_hs = mtx_hs.loc[cell_type_order_hs, kegg_id_order]

# # Relabel KEGG IDs by their name
# mtx_cr.columns = [pathway_names[kegg_id] for kegg_id in mtx_cr.columns]
# mtx_hs.columns = [pathway_names[kegg_id] for kegg_id in mtx_hs.columns]

# # CMAP for GO terms with changed bottom color
# from matplotlib.colors import LinearSegmentedColormap
# original_cmap = sns.color_palette("crest_r", as_cmap=True)
# # Define the new color for the bottom value
# new_color = (1.0, 1.0, 1.0, 1.0)
# new_colors = original_cmap.reversed()(np.linspace(0, 1, 256))
# new_colors[0] = new_color
# # Create the new colormap
# cmap_for_go = LinearSegmentedColormap.from_list("new_crest", new_colors)

# # Plot
# f, axs = plt.subplots(1, 2, figsize=(6, 6.5),
#                       width_ratios=[mtx_cr.shape[0], mtx_hs.shape[0]])
# f_cbar, ax_cbar = plt.subplots(1, 1, figsize=(0.75, 1.5))
# vmax = mtx_cr.max().max()
# sns.heatmap(mtx_cr.T, ax=axs[0], xticklabels=True, yticklabels=True,
#             cmap=cmap_for_go, linewidths=0.05, linecolor='#888888',
#             cbar=False, vmax=vmax)
# sns.heatmap(mtx_hs.T, ax=axs[1], xticklabels=True, yticklabels=False,
#             cmap=cmap_for_go, linewidths=0.05, linecolor='#888888',
#             cbar_ax=ax_cbar, cbar_kws={'label': 'Enrichment'}, vmax=vmax)
# # axs[0].set_xlabel('KEGG PATHWAY')
# axs[0].set_xlabel('C. robusta cell state')
# axs[1].set_xlabel('H. sapiens cell state')
# axs[0].xaxis.set_label_position('top')
# axs[1].xaxis.set_label_position('top')
# axs[0].tick_params(axis='x', which='major', rotation=90,
#                    labelbottom=False, bottom=False, top=True, labeltop=True)
# axs[1].tick_params(axis='x', which='major', rotation=90,
#                    labelbottom=False, bottom=False, top=True, labeltop=True)
# for ax in axs:
#     for _, spine in ax.spines.items():
#         spine.set_visible(True)
# f.tight_layout()
# f.savefig(out_path + 'overview_heatmap.pdf')
# f_cbar.tight_layout()
# f_cbar.savefig(out_path + 'overview_heatmap_cbar.pdf')
# plt.close('all')

# # ============================================================================
# # PLOT SUBSET FOR PAPER FIGURE

# pathway_subset_met = [
#     'Phenylalanine metabolism',
#     'Glycosaminoglycan degradation',
#     'Mannose type O-glycan biosynthesis',
#     'Pentose phosphate pathway',
# ]
# pathway_subset_imm = [
#     'Fc gamma R-mediated phagocytosis',
#     'Leukocyte transendothelial migration',
#     'Natural killer cell mediated cytotoxicity',
#     'C-type lectin receptor signaling pathway',
#     'RIG-I-like receptor signaling pathway',
#     'Chemokine signaling pathway',
# ]
# pathway_subset_pro = [
#     'Apoptosis',
#     'Necroptosis',
#     'Lysosome',
#     'Endocytosis',
#     'Focal adhesion',
# ]
# pathway_subset_sig = [
#     'TNF signaling pathway',
#     'NF-kappa B signaling pathway',
#     'VEGF signaling pathway',
#     'MAPK signaling pathway',
#     # 'cAMP signaling pathway',
#     # 'Phospholipase D signaling pathway',
#     'Oxytocin signaling pathway',
# ]
# # pathway_subset_ids = [x for x in pathway_names
# #                       if pathway_names[x] in pathway_subset]
# pathway_subset_list = [
#     pathway_subset_imm,
#     pathway_subset_pro,
#     pathway_subset_sig,
#     pathway_subset_met
# ]
# pathway_subset_together = []
# for pw in pathway_subset_list:
#     pathway_subset_together += pw

# # ------------------------------------
# # Plot heatmap

# mtx_cr = mtx_cr[pathway_subset_together]
# mtx_hs = mtx_hs[pathway_subset_together]
# mtx_cr = mtx_cr[mtx_cr.sum(axis=1) > 0]
# mtx_hs = mtx_hs[mtx_hs.sum(axis=1) > 0]
# mtx_cr = mtx_cr.loc[[x for x in adata_cr.obs['cell_type'].cat.categories
#                      if x in mtx_cr.index]]
# mtx_hs = mtx_hs.loc[[x for x in adata_hs.obs['cell_type_coarse'].cat.categories
#                      if x in mtx_hs.index]]


# def order_rows_by_col_max(mtx):
#     marginals = mtx.max(axis=0)
#     col_order = mtx.columns
#     row_max_match_to_col = mtx.idxmax(axis=1)
#     row_order = []
#     for iCol in col_order:
#         these_rows = row_max_match_to_col.index[
#             np.where(row_max_match_to_col == iCol)[0]]
#         marginals = mtx.loc[these_rows].max(axis=1)
#         these_rows_ordered = these_rows[np.argsort(marginals)[::-1]]
#         row_order += list(these_rows_ordered)
#     return row_order


# # Order rows and columns
# # ...Ciona cell types
# cell_type_order_cr = order_rows_by_col_max(mtx_cr)
# # ...human cell types
# cell_type_order_hs = order_rows_by_col_max(mtx_hs)
# # Actually reorder matrices
# mtx_cr = mtx_cr.loc[cell_type_order_cr, :]
# mtx_hs = mtx_hs.loc[cell_type_order_hs, :]

# # Plot
# f, axs = plt.subplots(len(pathway_subset_list), 2, figsize=(4, 3),
#                       width_ratios=[mtx_cr.shape[0], mtx_hs.shape[0]],
#                       height_ratios=[len(pw) for pw in pathway_subset_list])
# f_cbar, ax_cbar = plt.subplots(1, 1, figsize=(0.75, 1.5))
# vmax = mtx_cr.max().max()
# for i, pathway_subset in enumerate(pathway_subset_list):
#     mtx_to_plot_cr = mtx_cr[pathway_subset]
#     mtx_to_plot_hs = mtx_hs[pathway_subset]

#     sns.heatmap(mtx_to_plot_cr.T, ax=axs[i, 0],
#                 xticklabels=True, yticklabels=True,
#                 cmap=cmap_for_go, linewidths=0.05, linecolor='#888888',
#                 cbar=False, vmax=vmax)
#     sns.heatmap(mtx_to_plot_hs.T, ax=axs[i, 1],
#                 xticklabels=True, yticklabels=False,
#                 cmap=cmap_for_go, linewidths=0.05, linecolor='#888888',
#                 cbar_ax=ax_cbar, cbar_kws={'label': 'Enrichment'}, vmax=vmax)
#     # axs[0].set_xlabel('KEGG PATHWAY')
#     if i == 0:
#         # axs[i, 0].set_xlabel('C. robusta cell state')
#         # axs[i, 1].set_xlabel('H. sapiens cell state')
#         axs[i, 0].xaxis.set_label_position('top')
#         axs[i, 1].xaxis.set_label_position('top')
#         axs[i, 0].tick_params(axis='x', which='major', rotation=90,
#                         labelbottom=False, bottom=False, top=True, labeltop=True)
#         axs[i, 1].tick_params(axis='x', which='major', rotation=90,
#                         labelbottom=False, bottom=False, top=True, labeltop=True)
#     else:
#         axs[i, 0].set_xticks([])
#         axs[i, 1].set_xticks([])
# for ax in axs.flatten():
#     for _, spine in ax.spines.items():
#         spine.set_visible(True)
# f.tight_layout()
# f.subplots_adjust(hspace=0.1)
# f.savefig(out_path + 'selected_heatmap.pdf')
# f_cbar.tight_layout()
# f_cbar.savefig(out_path + 'selected_heatmap_cbar.pdf')
# plt.close('all')

# # ============================================================================
# # GENES SHARED ACROSS GENE SETS

# def gene_set_intersections(this_cell, these_ids):
#     # Get gene sets in this cell state's DEGs
#     hs_homolog_sets = {}
#     cr_gene_sets = {}
#     for this_id in these_ids:
#         hs_homolog_list = []
#         cr_gene_list = []
#         for g_cr in pathway_dict_ciona[this_id]:
#             if g_cr in degs_cr[this_cell]:
#                 cr_gene_list.append(g_cr)
#                 these_hs_homologs = hf.ciona2human[g_cr]
#                 for g in these_hs_homologs:
#                     if g in pathway_dict[this_id]:
#                         if g not in hs_homolog_list:
#                             hs_homolog_list.append(g)
#         hs_homolog_sets[this_id] = set(hs_homolog_list)
#         cr_gene_sets[this_id] = set(cr_gene_list)

#     # Get intersection of all pathways
#     intersection_all_hs = hs_homolog_sets[these_ids[0]]
#     intersection_all_cr = cr_gene_sets[these_ids[0]]
#     if len(these_ids) > 1:
#         # Human genes
#         for this_id in these_ids[1:]:
#             intersection_all_hs = intersection_all_hs.intersection(
#                 hs_homolog_sets[this_id])
#         # Ciona genes
#         for this_id in these_ids[1:]:
#             intersection_all_cr = intersection_all_cr.intersection(
#                 cr_gene_sets[this_id])
#     homologs_of_intersection_cr = [hf.ciona2human[g] for g in intersection_all_cr]

#     print('\nKEGG PATHWAYs:\n  '
#           + '\n  '.join([pathway_names[x] for x in these_ids]) + '\n')
#     print('Human genes intersection:\n  ' + '\n  '.join(intersection_all_hs))
#     print('Ciona genes intersection:\n  ' + '\n  '.join(intersection_all_cr))
#     print('Human homologs of Ciona intersection:\n  '
#         + '\n  '.join([' / '.join(x) for x in homologs_of_intersection_cr]))

#     return hs_homolog_sets, cr_gene_sets

# # ------------------------------------

# this_cell = 'URG-3 (small)'
# these_ids = [x for x in pathway_names if (
#     (pathway_names[x] == 'Pentose phosphate pathway')
# )]
# hs_homolog_sets, cr_gene_sets = gene_set_intersections(this_cell, these_ids)

# this_cell = 'RSC'
# these_ids = [x for x in pathway_names if (
#     (pathway_names[x] == 'Mannose type O-glycan biosynthesis')
# )]
# hs_homolog_sets, cr_gene_sets = gene_set_intersections(this_cell, these_ids)
