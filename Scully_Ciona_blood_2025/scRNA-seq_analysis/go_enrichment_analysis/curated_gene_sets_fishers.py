import os, sys
import pickle
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import fdrcorrection
from gseapy import prerank

# Change this path to point to folder containing helper functions scripts
path_to_repo_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(os.path.join(path_to_repo_dir, 'helper_functions'))
import scrna_helper_functions as hf

# ============================================================================
# OUTPUT PATH

out_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        'curated_gene_sets_fishers_output')
if not os.path.exists(out_path): os.mkdir(out_path)

# ============================================================================

# Ciona dataset
adata_cr = sc.read_h5ad(hf.path.Crob_adata_file)

# ============================================================================
# DEGs for each species

try:
    with open(os.path.join(out_path, 'degs_dict.p'), 'rb') as f:
        degs_cr = pickle.load(f)

except:
    print(f'DEGs for C. robusta')
    sc.tl.rank_genes_groups(adata_cr, groupby='cell_type', method='wilcoxon')
    degs_cr = {}
    for cell_type in tqdm(adata_cr.obs['cell_type'].cat.categories):
        df = sc.get.rank_genes_groups_df(adata_cr, cell_type)
        these_genes = set(df['names'][(df['logfoldchanges'] > 2)
                                    * (df['pvals_adj'] < 0.05)])
        degs_cr[cell_type] = these_genes

    with open(os.path.join(out_path, 'degs_dict.p'), 'wb') as f:
        pickle.dump(degs_cr, f)

# ============================================================================
# IMPORT CURATED GENE SETS

gene_set_df = pd.read_csv(os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    'curated_gene_lists.csv'
))
gene_set_dict = {}
gene_set_dict_human = {}
for gene_set in gene_set_df['Gene set name'].unique():
    these_rows = np.where(gene_set_df['Gene set name'] == gene_set)[0]
    gene_set_dict[gene_set] = [
        x for x in gene_set_df.loc[these_rows, 'Ciona gene(s)'] if x!='none'
    ]
    gene_set_dict_human[gene_set] = [
        gene_set_df.loc[i, 'Human homolog(s) or gene description'] for i in these_rows
        if gene_set_df.loc[i, 'Ciona gene(s)']!='none'
    ]
    if len(gene_set_dict[gene_set]) == 0:
        del gene_set_dict[gene_set]
        del gene_set_dict_human[gene_set]

# ============================================================================
# HYPERGEOMETRIC TEST / FISHER'S EXACT TEST

species_data_list = [
    ('cr', adata_cr, degs_cr, gene_set_dict),
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
            'Gene set': [],
            'Enrichment': [],
            'p-value': [],
            'num_genes_overlap': [],
            'num_genes_in_gene_set': [],
        }

        # Loop through every pair of clusters / cell types
        count = 0
        for i, cell_type in tqdm(enumerate(cell_type_list)):
            # print(f'{cell_type} ({i}/{len(cell_type_list)})')
            for gene_set_name in sp_pathway_dict:
                # Get lists of genes
                genes_this_set = set(sp_pathway_dict[gene_set_name])

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
                # fisher_results.loc[count, 'Cell type'] = cell_type
                # fisher_results.loc[count, 'KEGG pathway'] = kegg_term
                # fisher_results.loc[count, 'Enrichment'] = enrichment
                # fisher_results.loc[count, 'p-value'] = pval
                fisher_results['Cell type'].append(cell_type)
                fisher_results['Gene set'].append(gene_set_name)
                fisher_results['Enrichment'].append(enrichment)
                fisher_results['p-value'].append(pval)
                fisher_results['num_genes_overlap'].append(a)
                fisher_results['num_genes_in_gene_set'].append(a + c)

                count += 1

        fisher_results = pd.DataFrame(fisher_results)
        fisher_results['Cell type'] = fisher_results['Cell type'].astype('category')
        fisher_results['Gene set'] = fisher_results['Gene set'].astype('category')

        # --------------------------------
        # MULTIPLE HYPOTHESIS CORRECTION

        # Run FDR Benjamini-Hochberg correction
        [rejected, pval_fdr] = fdrcorrection(fisher_results['p-value'],
                                            alpha=0.05)

        fisher_results['null rejected'] = rejected
        fisher_results['FDR'] = pval_fdr

        # Save results
        fisher_results[fisher_results['null rejected']].to_csv(
            os.path.join(out_path2, f'curated_gene_sets_fisher_results.tsv'), sep='\t'
        )
