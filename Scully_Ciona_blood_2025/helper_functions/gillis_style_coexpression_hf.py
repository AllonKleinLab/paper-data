import os, sys
import time
import itertools
import numpy as np
import matplotlib.pyplot as plt
import scanpy as sc
from tqdm import tqdm
from scipy.stats import rankdata

# Change this path to point to folder containing gene_hf.py
# This imports dictionaries and functions for easily converting gene ids
path_to_dropbox = os.environ['PATH_TO_DROPBOX']
sys.path.append(path_to_dropbox + 'klein_lab/resources/helper_functions')
import scrna_helper_functions as hf

# ============================================================================


def my_spearmanr(matrix, axis=0):
    '''
    It seems that scipy.stats.spearmanr returns a nan instead of any results
    if too many columns are all 0s. This function will not, returning nans for
    columns of all 0s alongside actual results for non-zero columns.
    '''
    # Rank the data along each column
    ranked_matrix = np.apply_along_axis(rankdata, axis, matrix)
    
    # Take correlation of ranked matrix
    if axis == 0:
        corr_matrix = np.corrcoef(ranked_matrix.T)
    elif axis == 1:
        corr_matrix = np.corrcoef(ranked_matrix)

    return corr_matrix


def pseudobulk_data(adata, joint_gene_list, res=20):
    sc.tl.leiden(adata, resolution=res, key_added=f'leiden_res{res}')
    pseudobulk = hf.cluster_centroids(adata, groupby=f'leiden_res{res}')
    pseudobulk_arr = np.array(pseudobulk[joint_gene_list])
    return pseudobulk_arr


def coexpr_network_from_pseudobulk(pseudobulk):
    network_spearman = my_spearmanr(pseudobulk, axis=0)

    # Replace NaNs with median value on a per-row (per-gene) basis
    for i, row in enumerate(network_spearman):
        median_value = np.nanmedian(row)
        if np.isnan(median_value):  # happens if whole row is nans
            median_value = 0
        network_spearman[i, :] = np.nan_to_num(row, nan=median_value)
    # Remove connection between gene and itself
    for i in range(network_spearman.shape[0]):
        network_spearman[i, i] = 0

    # Rank
    # ranked_network = np.apply_along_axis(rankdata, 1, network_spearman)
    ranked_network = rankdata(network_spearman).reshape(network_spearman.shape)

    return ranked_network


def coexpr_network(adata, gene_list, res=20):
    pseudobulk_arr = pseudobulk_data(adata, gene_list, res=res)
    ranked_network = coexpr_network_from_pseudobulk(pseudobulk_arr)
    return ranked_network


def get_top_n(ranked_network, n=10):
    # row, col = np.where(ranked_network > (np.nanmax(ranked_network) - n))
    # top_n_genes = np.zeros(ranked_network.shape, dtype='bool')
    # top_n_genes[row, col] = True
    
    # Boolean matrix to store positions of top 10 values
    boolean_matrix = np.zeros(ranked_network.shape, dtype=bool)

    # Get the indices of the top 10 values for each row
    top_n_indices = np.argsort(ranked_network, axis=1)[:, -n:]

    # Set True for the positions of top 10 values
    for i in range(ranked_network.shape[0]):
        boolean_matrix[i, top_n_indices[i]] = True
    
    return boolean_matrix


def coexpression_auroc_for_one_gene_pair(i, j, top_n_genes_list,
                                         order_by_sp_list, joint_gene_list,
                                         example_gene_to_plot=[],
                                         save_as=''):
    """
    i = gene index for species 1 gene
    j = gene index for species 2 gene
    """
    plot_this_gene_bool = ((i==j) and
                           (joint_gene_list[i] in example_gene_to_plot))

    if plot_this_gene_bool:
        f = plt.figure(figsize=(5, 3))

    # ------------------------------------
    # Species 1 >> 2

    top_n_genes1 = top_n_genes_list[0]
    order_by_sp2 = order_by_sp_list[1]

    # Order by top genes in second species
    sort_top_genes1_by_sp2 = top_n_genes1[i, order_by_sp2[j, :]]

    # Calculate AUROC
    # true_pos_count = sort_top_genes1_by_sp2
    # false_pos_count = 1 - sort_top_genes1_by_sp2
    tpr = (np.cumsum(sort_top_genes1_by_sp2)
        / np.sum(sort_top_genes1_by_sp2))
    fpr = (np.cumsum(1 - sort_top_genes1_by_sp2)
        / np.sum(1 - sort_top_genes1_by_sp2))
    # integrate by trapezoidal sum
    auroc1 = np.trapz(y=tpr, x=fpr)

    if plot_this_gene_bool:
        ax = plt.subplot(1, 2, 1)
        ax.plot([0, 1], [0, 1], '--', color='#cccccc')
        ax.plot(fpr, tpr, 'k')
        ax.set_xlabel('FPR')
        ax.set_ylabel('TPR')
        ax.set_xlim((-0.01, 1.01)); ax.set_ylim((-0.01, 1.01))
        ax.text(1, 0, f'AUROC={auroc1:.2f}', transform=ax.transAxes,
                ha='right', va='bottom',)

    # ------------------------------------
    # Species 2 >> 1

    top_n_genes2 = top_n_genes_list[1]
    order_by_sp1 = order_by_sp_list[0]

    # Order by top genes in second species
    sort_top_genes2_by_sp1 = top_n_genes2[j, order_by_sp1[i, :]]

    # Calculate AUROC
    # true_pos_count = sort_top_genes2_by_sp1
    # false_pos_count = 1 - sort_top_genes2_by_sp1
    tpr = (np.cumsum(sort_top_genes2_by_sp1)
        / np.sum(sort_top_genes2_by_sp1))
    fpr = (np.cumsum(1 - sort_top_genes2_by_sp1)
        / np.sum(1 - sort_top_genes2_by_sp1))
    # integrate by trapezoidal sum
    auroc2 = np.trapz(y=tpr, x=fpr)

    if plot_this_gene_bool:
        ax = plt.subplot(1, 2, 2)
        ax.plot([0, 1], [0, 1], '--', color='#cccccc')
        ax.plot(fpr, tpr, 'k')
        ax.set_xlabel('FPR')
        ax.set_ylabel('TPR')
        ax.set_xlim((-0.01, 1.01)); ax.set_ylim((-0.01, 1.01))
        ax.text(1, 0, f'AUROC={auroc2:.2f}', transform=ax.transAxes,
                ha='right', va='bottom',)

    # ------------------------------------
    # Average of both directions
    auroc_mean = np.mean([auroc1, auroc2])

    if plot_this_gene_bool:
        plt.suptitle(f'Mean AUROC={auroc_mean:.2f}')
        plt.tight_layout()
        plt.savefig(save_as + f'{joint_gene_list[i]}_auroc.pdf')
        plt.close()

    return auroc_mean


def coexpression_auroc_for_one_gene_pair_complex_orthology(
        i, j, top_n_genes_list, order_by_sp_list, gene_list1, gene_list2,
        orthology_1to2, orthology_2to1, example_gene_to_plot=[], save_as=''):
    """
    i = gene index for species 1 gene
    j = gene index for species 2 gene

    top_n_genes_list : list of numpy arrays, outputs from get_top_n for each
        species
    order_by_sp_list : list of numpy arrays, rank-ordered genes in order of
        decreasing correlation
    gene_list1 : list, gene list for species 1
    gene_list2 : list, gene list for species 2
    orthology_1to2 : gene orthology dict from species 1 to 2, objects in dict
        are lists of orthologous genes
    orthology_2to1 : gene orthology dict from species 2 to 1, objects in dict
        are lists of orthologous genes

    example_gene_to_plot : list of tuples of strings, (sp1 gene, sp2 gene)
    save_as : str, used for filenaming
    """
    plot_this_gene_bool = ((gene_list1[i], gene_list2[j])
                           in example_gene_to_plot)

    if plot_this_gene_bool:
        f = plt.figure(figsize=(5, 3))

    # ------------------------------------
    # Species 1 >> 2

    top_n_genes1 = top_n_genes_list[0]
    order_by_sp2 = order_by_sp_list[1]

    # Order by top genes in second species
    top_n_coexpr_factors1 = gene_list1[np.where(top_n_genes1[i])[0]]
    sp2_orthologs_of_top_n = []
    for g in top_n_coexpr_factors1:
        if g in orthology_1to2:
            sp2_orthologs_of_top_n += orthology_1to2[g]

    ordered_sp2_genes = gene_list2[order_by_sp2[j]]
    ordered_sp2_genes_bool = np.array([(g in sp2_orthologs_of_top_n)
                                       for g in ordered_sp2_genes])

    # Calculate AUROC
    tpr = (np.cumsum(ordered_sp2_genes_bool)
           / np.sum(ordered_sp2_genes_bool))
    fpr = (np.cumsum(1 - ordered_sp2_genes_bool)
           / np.sum(1 - ordered_sp2_genes_bool))
    # integrate by trapezoidal sum
    auroc1 = np.trapz(y=tpr, x=fpr)

    if plot_this_gene_bool:
        ax = plt.subplot(1, 2, 1)
        ax.plot([0, 1], [0, 1], '--', color='#cccccc')
        ax.plot(fpr, tpr, 'k')
        ax.set_xlabel('FPR')
        ax.set_ylabel('TPR')
        ax.set_xlim((-0.01, 1.01)); ax.set_ylim((-0.01, 1.01))
        ax.text(1, 0, f'sp1 to sp2\nAUROC={auroc1:.2f}',
                transform=ax.transAxes, ha='right', va='bottom',)

    # ------------------------------------
    # Species 2 >> 1

    top_n_genes2 = top_n_genes_list[1]
    order_by_sp1 = order_by_sp_list[0]

    # Order by top genes in second species
    top_n_coexpr_factors2 = gene_list2[np.where(top_n_genes2[j])[0]]
    sp1_orthologs_of_top_n = []
    for g in top_n_coexpr_factors2:
        if g in orthology_2to1:
            sp1_orthologs_of_top_n += orthology_2to1[g]

    ordered_sp1_genes = gene_list1[order_by_sp1[i]]
    ordered_sp1_genes_bool = np.array([(g in sp1_orthologs_of_top_n)
                                       for g in ordered_sp1_genes])

    # Calculate AUROC
    tpr = (np.cumsum(ordered_sp1_genes_bool)
           / np.sum(ordered_sp1_genes_bool))
    fpr = (np.cumsum(1 - ordered_sp1_genes_bool)
           / np.sum(1 - ordered_sp1_genes_bool))
    # integrate by trapezoidal sum
    auroc2 = np.trapz(y=tpr, x=fpr)

    if plot_this_gene_bool:
        ax = plt.subplot(1, 2, 2)
        ax.plot([0, 1], [0, 1], '--', color='#cccccc')
        ax.plot(fpr, tpr, 'k')
        ax.set_xlabel('FPR')
        ax.set_ylabel('TPR')
        ax.set_xlim((-0.01, 1.01)); ax.set_ylim((-0.01, 1.01))
        ax.text(1, 0, f'sp2 to sp1\nAUROC={auroc2:.2f}',
                transform=ax.transAxes, ha='right', va='bottom',)

    # ------------------------------------
    # Average of both directions
    auroc_mean = np.mean([auroc1, auroc2])

    if plot_this_gene_bool:
        plt.suptitle(f'Mean AUROC={auroc_mean:.2f}')
        plt.tight_layout()
        plt.savefig(save_as + f'{gene_list1[i]}_{gene_list2[j]}_auroc.pdf')
        plt.close()

    return auroc_mean


def generate_unique_pairs(n_pairs, range1, range2, forbid_same_num=True,
                          forbidden_pairs=None):
    """
    Get a list of n_pairs tuples, where:
    - first number in each tuple is between range1[0] and range1[1]
    - second number in each tuple is between range2[0] and range2[1]
    - numbers in tuples are unique, i.e. no (0, 0) or (1, 1)
    - there are no repeat tuples in the list

    forbid_same_num : True if not allowing i=j for gene pairs, otherwise False

    Function adapted from chat GPT
    """
    result = set()
    while len(result) < n_pairs:
        # Generate random number pairs in bulk
        a = np.random.randint(range1[0], range1[1] + 1, size=n_pairs*2)
        b = np.random.randint(range2[0], range2[1] + 1, size=n_pairs*2)
        
        if forbid_same_num:
            # Filter out pairs where both elements are the same
            mask = a != b
            valid_pairs = zip(a[mask], b[mask])

            # Add to the set and discard duplicates
            result.update(valid_pairs)

        elif forbidden_pairs is not None:
            # Exclude forbidden pairs
            candidate_pairs = set(zip(a, b))
            valid_pairs = candidate_pairs - set(forbidden_pairs)

            # Add to the set
            result.update(valid_pairs)

        else:
            result.update(zip(a, b))

    # Convert the set back to a list and trim to exact size
    result = list(result)[:n_pairs]
    
    return result


def calculate_coexpression_auroc(ranked_network_list, joint_gene_list, n=10,
                                 example_gene_to_plot=[], save_as='',
                                 random_seed=0):
    num_genes = ranked_network_list[0].shape[0]

    # -----
    # 1. For each gene (row), get top n co-expressed genes
    print(f'Getting top coexpression partners for each gene in each species...')
    top_n_genes_list = []
    for i in range(len(ranked_network_list)):
        top_n_genes_list.append(get_top_n(ranked_network_list[i], n=n))

    # -----
    # 2. For each gene (row), get ordered genes from highest to lowest
    # coexpression in the second species
    print(f'Ordering by high to low coexpression in each species...')
    order_by_sp_list = []
    for i in range(len(ranked_network_list)):
        order_by_sp_list.append(np.argsort(ranked_network_list[i],
                                           axis=1)[:, ::-1])
    
    # -----
    # 3. Calculate AUROC

    print('Calculating two-way AUROCs...')

    # auroc = np.zeros((num_genes, num_genes))
    auroc_orthologs = []
    auroc_non_orthologs = []

    # Get random set of non-orthologous gene pairs for background
    num_random_pairs = 100*num_genes
    np.random.seed(random_seed)
    random_non_orthologous_pairs = generate_unique_pairs(
        num_random_pairs, [0, num_genes-1], [0, num_genes-1])

    print('  orthologous gene pairs')
    for i in tqdm(range(num_genes)):
        auroc_orthologs.append(coexpression_auroc_for_one_gene_pair(
            i=i,
            j=i,    # set j = i since we only look at orthologous pairs
            top_n_genes_list=top_n_genes_list,
            order_by_sp_list=order_by_sp_list,
            joint_gene_list=joint_gene_list,
            example_gene_to_plot=example_gene_to_plot,
            save_as=save_as,
        ))
    
    print(f'  {num_random_pairs} random non-orthologous gene pairs')
    for (i, j) in tqdm(random_non_orthologous_pairs):
        auroc_non_orthologs.append(coexpression_auroc_for_one_gene_pair(
            i=i,
            j=j,
            top_n_genes_list=top_n_genes_list,
            order_by_sp_list=order_by_sp_list,
            joint_gene_list=joint_gene_list,
            example_gene_to_plot=example_gene_to_plot,
            save_as=save_as,
        ))
    
    return auroc_orthologs, auroc_non_orthologs


def calculate_coexpression_auroc_complex_orthology(
        ranked_network_list, gene_list1, gene_list2, orthology_1to2,
        orthology_2to1, n=10, example_gene_to_plot=[], save_as='',
        random_seed=0, num_random_pairs=None,
        orthology_dict=hf.ciona2human):

    # -----
    # 1. For each gene (row), get top n co-expressed genes
    print(f'Getting top coexpression partners for each gene in each species...')
    top_n_genes_list = []
    for i in range(len(ranked_network_list)):
        top_n_genes_list.append(get_top_n(ranked_network_list[i], n=n))

    # -----
    # 2. For each gene (row), get ordered genes from highest to lowest
    # coexpression in the second species
    print(f'Ordering by high to low coexpression in each species...')
    order_by_sp_list = []
    for i in range(len(ranked_network_list)):
        order_by_sp_list.append(np.argsort(ranked_network_list[i],
                                           axis=1)[:, ::-1])
    
    # -----
    # 3. Get orthologous and non-orthologous gene pairs

    print('Getting orthologous gene lists...')
    full_orthologous_gene_list = []
    for i, cr_g in enumerate(gene_list1):
        hs_g_list = orthology_dict[cr_g]
        for hs_g in hs_g_list:
            if hs_g in gene_list2:
                j = np.where(hs_g == gene_list2)[0][0]
                if (i, j) not in full_orthologous_gene_list:
                    full_orthologous_gene_list.append((i, j))

    print('Getting random non-orthologous gene lists...')
    # Get random set of non-orthologous gene pairs for background
    if num_random_pairs is None:
        num_random_pairs = 100*len(full_orthologous_gene_list)
    np.random.seed(random_seed)
    random_non_orthologous_pairs = generate_unique_pairs(
        num_random_pairs, [0, len(gene_list1)-1], [0, len(gene_list2)-1],
        forbid_same_num=False, forbidden_pairs=full_orthologous_gene_list)

    # -----
    # 4. Calculate AUROC

    print('Calculating two-way AUROCs...')

    # auroc = np.zeros((num_genes, num_genes))
    auroc_orthologs = []
    auroc_non_orthologs = []

    print('  orthologous gene pairs')
    for i, j in tqdm(full_orthologous_gene_list):
        auroc_orthologs.append(
            coexpression_auroc_for_one_gene_pair_complex_orthology(
                i=i,
                j=j,
                top_n_genes_list=top_n_genes_list,
                order_by_sp_list=order_by_sp_list,
                gene_list1=gene_list1,
                gene_list2=gene_list2,
                orthology_1to2=orthology_1to2,
                orthology_2to1=orthology_2to1,
                example_gene_to_plot=example_gene_to_plot,
                save_as=save_as,
            )
        )
    
    print(f'  {num_random_pairs} random non-orthologous gene pairs')
    for (i, j) in tqdm(random_non_orthologous_pairs):
        auroc_non_orthologs.append(
            coexpression_auroc_for_one_gene_pair_complex_orthology(
                i=i,
                j=j,
                top_n_genes_list=top_n_genes_list,
                order_by_sp_list=order_by_sp_list,
                gene_list1=gene_list1,
                gene_list2=gene_list2,
                orthology_1to2=orthology_1to2,
                orthology_2to1=orthology_2to1,
                example_gene_to_plot=example_gene_to_plot,
                save_as=save_as,
            )
        )
    
    return auroc_orthologs, auroc_non_orthologs, full_orthologous_gene_list


def coexpression_conservation_from_auroc(auroc_orthologs, auroc_non_orthologs,
                                         gene_list,
                                         example_gene_to_plot=[], save_as=''):
    print('Calculating coexpression conservation score...')

    coexpr_conservation = []
    # Combine auroc lists, first N are the N genes
    auroc_all = np.array(auroc_orthologs + auroc_non_orthologs)
    auroc_no_nan = np.nan_to_num(auroc_all, nan=np.nanmedian(auroc_all))
    auroc_ranked = rankdata(auroc_no_nan, axis=None, method='dense')
    for iGene in tqdm(range(len(auroc_orthologs))):
        # This AUROC rank vs. all other gene pairs
        this_gene_auroc_rank = auroc_ranked[iGene]
        coexpr_conservation.append((this_gene_auroc_rank
                                    / np.max(auroc_ranked)))
        
        if gene_list[iGene] in example_gene_to_plot:
            this_gene_auroc = auroc_orthologs[iGene]

            plt.figure(figsize=(3, 3))
            ax = plt.subplot(1, 1, 1)
            freq, bins = np.histogram(auroc_no_nan, bins=50)
            plt.plot(bins[:-1], freq, 'k')
            # Red point where this gene is
            histogram_idx = np.where(this_gene_auroc <= bins)[0][0]
            plt.plot(this_gene_auroc, freq[histogram_idx-1], 'd',
                    color='#ff00aa', markeredgecolor='k')
            plt.xlabel('AUROC')
            plt.ylabel('Frequency')
            ax.set_ylim((0, ax.get_ylim()[1]))
            plt.title('Coexpression\nconservation='
                    + f'{coexpr_conservation[iGene]:.3f}')
            plt.tight_layout()
            plt.savefig(save_as + f'{gene_list[iGene]}_auroc_rank.pdf')
            plt.close()
    
    return coexpr_conservation


def calculate_coexpression_conservation(ranked_network_list, joint_gene_list,
                                        example_gene_to_plot=[], save_as=''):
    auroc_orthologs, auroc_non_orthologs = calculate_coexpression_auroc(
        ranked_network_list, joint_gene_list,
        example_gene_to_plot=example_gene_to_plot, save_as=save_as)
    coexpr_conservation = coexpression_conservation_from_auroc(
        auroc_orthologs, auroc_non_orthologs, joint_gene_list,
        example_gene_to_plot=example_gene_to_plot, save_as=save_as)

    return coexpr_conservation


def calculate_coexpression_conservation_complex_orthology(
        ranked_network_list, gene_list1, gene_list2, orthology_1to2,
        orthology_2to1, n=20, num_random_pairs=None,
        example_gene_to_plot=[], save_as='', orthology_dict=hf.ciona2human):
    auroc_orthologs, auroc_non_orthologs, gene_pairs = \
        calculate_coexpression_auroc_complex_orthology(
            ranked_network_list, gene_list1, gene_list2, orthology_1to2,
            orthology_2to1, n=n, num_random_pairs=num_random_pairs,
            example_gene_to_plot=example_gene_to_plot, save_as=save_as,
            orthology_dict=orthology_dict,
        )
    coexpr_conservation = coexpression_conservation_from_auroc(
        auroc_orthologs, auroc_non_orthologs, gene_pairs,
        example_gene_to_plot=example_gene_to_plot, save_as=save_as)
    
    outputs = {
        'coexpr_conservation': coexpr_conservation,
        'auroc_orthologs': auroc_orthologs,
        'auroc_non_orthologs': auroc_non_orthologs,
        'gene_pairs': gene_pairs,
    }

    return outputs
