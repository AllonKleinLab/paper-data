import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import gseapy as gp
import sys
import os
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import fdrcorrection
from tqdm import tqdm

# Change this path to point to folder containing helper functions scripts
path_to_repo_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(os.path.join(path_to_repo_dir, 'helper_functions'))
import scrna_helper_functions as hf

# ============================================================================


def reorder_matrix(centr, row_order='none',
                  col_order='none',
                  method_row='single', metric_row='cosine',
                  method_col='single', metric_col='cosine',
                  count_sort_row='descending',
                  count_sort_col='ascending',
                  save_as=''):
    # For hierarchical clustering, use a matrix with no nan or inf
    centr2 = centr.copy()
    centr2.replace(np.inf, 300, inplace=True)
    centr2.replace(np.nan, 300, inplace=True)

    # 1. REORDER CELL TYPES
    # a. hierarchical clustering
    if row_order == 'hierarchical_clustering':
        with plt.style.context('tal_paper'):
            f = plt.figure(figsize=(3, 3))
            ax = plt.subplot(1, 1, 1)
            dendr = hf.hierarchical_clustering(
                centr2, ax, optimal_ordering=True,
                labels=centr.index, count_sort=count_sort_row,
                method=method_row, metric=metric_row,
            )
            plt.tight_layout()
            plt.savefig(out_path + f'{save_as}_cell_dendrogram.pdf')
            plt.close()
            
        cluster_order = dendr['ivl'][::-1]
        centr = centr.loc[cluster_order]

    # b. manual order
    elif isinstance(row_order, list):
        centr = centr.loc[row_order]

    # c. descending row mean
    elif row_order == 'descending_row_mean':
        row_order_idx = np.argsort(centr.mean(axis=1))[::-1]
        row_order = centr.index[row_order_idx]
        centr = centr.loc[row_order]

    # 1. REORDER GENES
    # a. hierarchical clustering
    if col_order == 'hierarchical_clustering':
        with plt.style.context('tal_paper'):
            f = plt.figure(figsize=(3, 3))
            ax = plt.subplot(1, 1, 1)
            dendr = hf.hierarchical_clustering(
                centr2.T, ax, optimal_ordering=True,
                labels=centr.columns, count_sort=count_sort_col,
                method=method_col, metric=metric_col,
            )
            plt.tight_layout()
            plt.savefig(out_path + f'{save_as}_gene_dendrogram.pdf')
            plt.close()
            
        # For gene ordering, put not expressed genes at the end
        col_order = dendr['ivl'][::-1]
        centr = centr[col_order]
    
    # b. manual order
    elif isinstance(col_order, list):
        centr = centr.loc[col_order]

    # c. descending row mean
    elif col_order == 'descending_row_mean':
        col_order_idx = np.argsort(centr.mean(axis=1))[::-1]
        col_order = centr.index[col_order_idx]
        centr = centr.loc[col_order]

    return centr


# ============================================================================
# IMPORT DATA

# Import counts matrix
adata = sc.read_h5ad(hf.path.Crob_adata_file)

# Set up output folder
out_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        'go_enrichment_analysis_output')
if not os.path.exists(out_path): os.mkdir(out_path)

# ============================================================================
# GO term gene sets with Ciona genes

for go_gene_sets in ['GO_Molecular_Function_2023', 'GO_Biological_Process_2023',
                     'GO_Cellular_Component_2023']:
    print(go_gene_sets)

    go_human = gp.get_library(name=go_gene_sets, organism='human')

    # GO terms to look at
    go_list = []

    go_ciona = {}
    for term in tqdm(go_human):
        go_ciona[term] = []
        this_gene_set = go_human[term]
        
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
                    if ciona_g not in go_ciona[term]:
                        go_ciona[term].append(ciona_g)
        
        # If >=10 human genes have orthologs
        if count_human_genes_with_orthologs >= 10:
            # AND if >=10 ciona genes in final list
            if len(go_ciona[term]) >= 10:
                go_list.append(term)


    # ============================================================================
    # Differentially expressed genes (Wilcoxon rank-sum test)

    # Wilcoxon rank-sum test
    sc.tl.rank_genes_groups(adata, groupby='cell_type')

    # ============================================================================
    # HYPERGEOMETRIC TEST / FISHER'S EXACT TEST

    try:
        fisher_results = pd.read_csv(os.path.join(out_path, 'fisher_results.tsv'),
                                     sep='\t')
        fisher_results['Cell type'] = fisher_results['Cell type'].astype('category')
        fisher_results['GO term'] = fisher_results['GO term'].astype('category')

    except:

        # Cell types to look at
        cell_type_list = adata.obs['cell_type'].cat.categories

        # Define the world (all genes which could be in a set)
        gene_bool = np.array(np.sum(adata.raw.X != 0, axis=0) > 20).flatten()
        # gene_bool = adata.var['highly_variable']
        full_set = adata.var.index[gene_bool]

        # Initialize pandas dataframe to save info
        fisher_results = {
            'Cell type': [],
            'GO term': [],
            'Enrichment': [],
            'p-value': [],
            'num_genes_overlap': [],
            'num_genes_in_GO_set': [],
        }

        # Loop through every pair of clusters / cell types
        count = 0
        for i, cell_type in enumerate(cell_type_list):
            print(f'{cell_type} ({i}/{len(cell_type_list)})')
            for go_term in tqdm(go_list):

                # Get lists of genes
                genes_go = set(go_ciona[go_term])

                # Get list of genes for this cell type
                wilcoxon_df = sc.get.rank_genes_groups_df(adata, group=cell_type)
                genes_c = set(wilcoxon_df.loc[
                    (
                        (wilcoxon_df['logfoldchanges'] > 1)
                        * (wilcoxon_df['pvals_adj']< 0.01)
                    ),
                    'names'
                ].values)

                # Filter so gene sets not in "world" are ignored
                genes_go = {g for g in genes_go if g in full_set}
                genes_c = {g for g in genes_c if g in full_set}

                # Create contingency table
                # Yes Ciona cell type and yes GO term
                a = len(genes_c.intersection(genes_go))
                # Yes Ciona cell type and not GO term
                b = len(genes_c.difference(genes_go))
                # Not Ciona cell type and yes GO term
                c = len(genes_go.difference(genes_c))
                # Not Ciona cell type and not GO term
                # a + b + c + d = len(full_set)
                d = len(full_set) - a - b - c

                contingency_table = np.array([[a, b], [c, d]])

                # Perform Fisher's exact test
                [enrichment, pval] = fisher_exact(contingency_table,
                                                alternative='greater')
                # fisher_results.loc[count, 'Cell type'] = cell_type
                # fisher_results.loc[count, 'GO term'] = go_term
                # fisher_results.loc[count, 'Enrichment'] = enrichment
                # fisher_results.loc[count, 'p-value'] = pval
                fisher_results['Cell type'].append(cell_type)
                fisher_results['GO term'].append(go_term)
                fisher_results['Enrichment'].append(enrichment)
                fisher_results['p-value'].append(pval)
                fisher_results['num_genes_overlap'].append(a)
                fisher_results['num_genes_in_GO_set'].append(a + c)

                count += 1

        fisher_results = pd.DataFrame(fisher_results)
        fisher_results['Cell type'] = fisher_results['Cell type'].astype('category')
        fisher_results['GO term'] = fisher_results['GO term'].astype('category')

        # --------------------------------
        # MULTIPLE HYPOTHESIS CORRECTION

        # Run FDR Benjamini-Hochberg correction
        [rejected, pval_fdr] = fdrcorrection(fisher_results['p-value'],
                                            alpha=0.05)

        fisher_results['null rejected'] = rejected
        fisher_results['FDR'] = pval_fdr

        # Save results
        fisher_results[fisher_results['null rejected']].to_csv(
            os.path.join(out_path, f'{go_gene_sets}_fisher_results.tsv'), sep='\t'
        )
