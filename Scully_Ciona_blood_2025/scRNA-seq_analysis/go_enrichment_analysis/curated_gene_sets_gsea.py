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
                        'curated_gene_sets_gsea_output')
if not os.path.exists(out_path): os.mkdir(out_path)

# ============================================================================

# Ciona dataset
adata_cr = sc.read_h5ad(hf.path.Crob_adata_file)

# ============================================================================
# DEGs for each species

# except:
print(f'DEGs for C. robusta')
sc.tl.rank_genes_groups(adata_cr, groupby='cell_type', method='wilcoxon')

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
# GSEA

gsea_results = []
for cell_state in adata_cr.obs['cell_type'].cat.categories:
    # Rank genes based on Wilcoxon rank-sum test
    result = adata_cr.uns['rank_genes_groups']
    genes = pd.DataFrame(result['names'])[cell_state]
    log_fc = pd.DataFrame(result['scores'])[cell_state]
    rnk = pd.DataFrame({'gene': genes, 'score': log_fc})
    rnk = rnk.sort_values(by='score', ascending=False)

    # Run GSEA
    pre_res = prerank(
        rnk=rnk,
        gene_sets=gene_set_dict,
        permutation_num=1000,
        min_size=5
    )

    # Plot
    term2 = pre_res.res2d.Term
    axes = pre_res.plot(terms=term2[:5])
    # plt.subplots_adjust(top=0.9, bottom=0.1, right=0.75, left=0.25)
    plt.savefig(os.path.join(out_path, f'gseapy_{cell_state.replace("/", "-")}.pdf'),
                bbox_inches='tight')
    plt.close()

    # Save data
    for i in range(pre_res.res2d.shape[0]):
        row = pre_res.res2d.iloc[i, :]
        if (row['FDR q-val'] < 0.05) or (row['FWER p-val'] < 0.05):
            gsea_results.append([
                cell_state,
                row['Term'],
                row['NES'],
                row['FDR q-val'],
                row['FWER p-val'],
                row['Lead_genes'],
            ])
gsea_results = pd.DataFrame(gsea_results, columns=['Cell state', 'Gene set',
                                                   'NES', 'FDR', 'FWER', 'Lead_genes'])

# Save results
gsea_results[gsea_results['FWER'] < 0.05].to_csv(
    os.path.join(out_path, f'curated_gene_sets_gsea_results.tsv'), sep='\t'
)
