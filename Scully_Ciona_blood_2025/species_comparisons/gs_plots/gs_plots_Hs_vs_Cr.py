import os, sys
import pickle
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm
from itertools import combinations

# Change this path to point to folder containing helper functions scripts
path_to_repo_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(os.path.join(path_to_repo_dir, 'helper_functions'))
import scrna_helper_functions as hf
from gs_plot_helper_functions import *

# Plot style
plt.style.use('tal_paper')

# ============================================================================
# OUTPUT PATH

out_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        'gs_plots_Hs_vs_Cr_output')
if not os.path.exists(out_path): os.mkdir(out_path)

# ============================================================================
# THESE SPECIES

species1 = hf.Species('Ciona', 'robusta')
species2 = hf.Species('Homo', 'sapiens')

sp1_color = '#008c67'
sp2_color = '#8918d3'

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

print(f'DEGs for {species1.name}')
sc.tl.rank_genes_groups(adata_cr, groupby='cell_type', method='wilcoxon')
degs_cr = {}
for cell_type in adata_cr.obs['cell_type'].cat.categories:
    df = sc.get.rank_genes_groups_df(adata_cr, cell_type)
    these_genes = set(df['names'][(df['logfoldchanges'] > 2)
                                  * (df['pvals_adj'] < 0.05)])
    degs_cr[cell_type] = these_genes

print(f'DEGs for {species2.name}')
sc.tl.rank_genes_groups(adata_hs, groupby='cell_type_coarse',
                        method='wilcoxon')
degs_hs = {}
for cell_type in adata_hs.obs['cell_type_coarse'].cat.categories:
    df = sc.get.rank_genes_groups_df(adata_hs, cell_type)
    these_genes = set(df['names'][(df['logfoldchanges'] > 2)
                                  * (df['pvals_adj'] < 0.05)])
    degs_hs[cell_type] = these_genes

# ============================================================================
# GENE HOMOLOGY FOR HUMAN VS. CIONA

cr_gene_list = adata_cr.var.index
hs_gene_list = adata_hs.var.index

print('\nGetting OrthoFinder orthologs for all gene pairs...')
orthofinder_all = []
for cr_gene in hf.ciona2human:
    for hs_gene in hf.ciona2human[cr_gene]:
        if (cr_gene in cr_gene_list) and (hs_gene in hs_gene_list):
            orthofinder_all.append([cr_gene, hs_gene, 1])

# ============================================================================
# Functions for getting numbers needed for modified Sankey plots


def get_gene_subsets(gene_sets):
    """
    Given a dictionary of gene sets per cell state, returns a dictionary 
    giving the set of genes in each combination. For example, for 3 cell
    states, gives the set of genes which are:
    - only in state1
    - only in state2
    - only in state3
    - in state1 and state2, but not state3
    - in state2 and state2, but not state1
    - in state1 and state3, but not state2
    - in all 3 of state1, state2, and state3
    This is used for determining how much ribbons in GS plots should overlap.

    Args:   
        gene_sets (dict) : {"ciona1": set1, "ciona2": set2, "ciona3": set3, ...}
    
    Returns:
        subsets (dict) : Keys are tuples of cell states, e.g. ('state1',) for
            genes only in state1, or ('state2', 'state3') for genes in state2
            and state3 but not state1. Values are sets, listing the genes
            corresponding to each state combination.
    """
    subsets = {}

    # Generate all non-empty combinations of set names
    set_names = list(gene_sets.keys())
    for r in range(1, len(set_names) + 1):
        for combo in combinations(set_names, r):
            # Get the intersection of the selected sets
            intersection_set = set.intersection(*(gene_sets[name]
                                                  for name in combo))
            
            # Get the union of sets not in this combination (to exclude genes
            # appearing elsewhere)
            other_sets = [gene_sets[name] for name in set_names
                          if name not in combo]
            exclusion_set = set.union(*other_sets) if other_sets else set()

            # Subtract genes that appear in any other set
            subsets[combo] = intersection_set - exclusion_set

    return subsets


def get_gene_dict_sets(sp1_cell:str, sp2_cell_list:list, sp1='human',
                       gene_orthology='orthofinder'):
    """
    Given the names of cell states to compare, get dictionaries containing
    shared genes (i.e. genes which have a homolog in the DEGs of the other
    species's cell types' DEGs)
    
    Args:
        sp1_cell (str) : The name of the single species1 cell state, i.e.
            the top rectangle in the GS plot.
        sp2_cell_list (str) : The names of the species 2 cell states, i.e.
            the bottom rectangles in teh GS plot.
        sp1 (str) : The name of the species whose cell state is on top,
            either 'human' or something else.
        gene_orthology (str) : The type of gene homology to use. We only
            use OrthoFinder (gene_orthology='orthofinder') in this script.
    
    Returns:
        hs_genes_dict (dict) : Keys are human cell state names, values are
            dictionaries with the following keys:
            - 'all' contains the set of enriched DEGs in this human cell state
            - 'shared_w_{xx}' contains the set of enriched DEGs which have a 
              homolog in the DEGs of the other species's cell state {xx}
        cr_genes_dict (dict) : Same set-up as hs_genes_dict, but with the
            species swapped.
    """
    # Initialize dictionaries
    if sp1 == 'human':
        hs_genes_dict = {sp1_cell: {}}
        cr_genes_dict = {cr_cell: {} for cr_cell in sp2_cell_list}
    else:
        cr_genes_dict = {sp1_cell: {}}
        hs_genes_dict = {hs_cell: {} for hs_cell in sp2_cell_list}

    # Determine orthology type
    if gene_orthology == 'orthofinder':
        orthology_pairs = orthofinder_all

    # Orthology dicts
    cr2hs = {}; hs2cr = {}
    for dr_g, hs_g, _ in orthology_pairs:
        # Dr to Hs conversion
        if dr_g not in cr2hs: cr2hs[dr_g] = [hs_g]
        else: cr2hs[dr_g].append(hs_g)
        # Hs to Dr conversion
        if hs_g not in hs2cr: hs2cr[hs_g] = [dr_g]
        else: hs2cr[hs_g].append(dr_g)

    for cr_cell in cr_genes_dict:
        for hs_cell in hs_genes_dict:
            # x = centr_z_dr.loc[cr_cell, [x[0] for x in orthology_pairs]]
            # y = centr_z_hs.loc[hs_cell, [x[1] for x in orthology_pairs]]
            # w = np.array([x[2] for x in orthology_pairs])

            # Counting number above threshold in both species
            # count_shared = np.sum((x >= th).values * (y >= th).values
            #                       * (w > th_orthology))
            # enriched_genes_dr = x[(x >= th).values * (w >= th_orthology)].index.unique()
            # enriched_genes_hs = y[(y >= th).values * (w >= th_orthology)].index.unique()
            # genes_dr_with_enriched_ortholog = x[
            #     (x >= th).values * (y >= th).values * (w > th_orthology)].index.unique()
            # genes_hs_with_enriched_ortholog = y[
            #     (x >= th).values * (y >= th).values * (w > th_orthology)].index.unique()
            enriched_genes_cr = set([g for g in degs_cr[cr_cell] if g in cr2hs])
            enriched_genes_hs = set([g for g in degs_hs[hs_cell] if g in hs2cr])
            genes_dr_w_enriched_ortholog = [
                g for g in enriched_genes_cr
                if len(set(cr2hs[g]).intersection(enriched_genes_hs)) > 0]
            genes_hs_w_enriched_ortholog = [
                g for g in enriched_genes_hs
                if len(set(hs2cr[g]).intersection(enriched_genes_cr)) > 0]

            cr_genes_dict[cr_cell]['all'] = enriched_genes_cr
            hs_genes_dict[hs_cell]['all'] = enriched_genes_hs
            cr_genes_dict[cr_cell][f'shared_w_{hs_cell}'] = set(
                genes_dr_w_enriched_ortholog)
            hs_genes_dict[hs_cell][f'shared_w_{cr_cell}'] = set(
                genes_hs_w_enriched_ortholog)

    return hs_genes_dict, cr_genes_dict


def get_sankey_plot_numbers(sp1_genes_dict, sp2_genes_dict):
    """
    Takes the gene dicts outputted by get_gene_dict_sets() and outputs the
    numbers in the format needed to make a GS plot.

    Args:
        sp1_genes_dict (dict) : Output of get_gene_dict_sets for whichever
            species has a single cell state (will become the top rectangle)
        sp2_genes_dict (dict) : Output of get_gene_dict_sets for whichever
            species has multiple cell states (will become the bottom
            rectangles)
    
    Returns:
        sp1_rect_dicts (dict) : Dictionary indicating the number of DEGs per
            species1 (top rectangle) cell state
        sp2_rect_dicts (dict) : Dictionary indicating the number of DEGs per
            species2 (bottom rectangles) cell states
        ribbon_widths (list) : Contains information for plotting GS plot flows
    """

    # Rectangle widths
    sp1_rects_dict = {x: len(sp1_genes_dict[x]['all']) for x in sp1_genes_dict}
    sp2_rects_dict = {x: len(sp2_genes_dict[x]['all']) for x in sp2_genes_dict}

    # Ribbon widths
    sp1_cell = [x for x in sp1_rects_dict][0]
    set_overlap_sp2 = get_gene_subsets({x: sp2_genes_dict[x][f'shared_w_{sp1_cell}']
                                       for x in sp2_genes_dict})
    set_overlap_sp1 = get_gene_subsets({x: sp1_genes_dict[sp1_cell][f'shared_w_{x}']
                                       for x in sp2_genes_dict})
    # [top_rect, bot_rect, top_width, bot_width, overlap_with_next]
    ribbon_widths = []
    for sp2_cell_set in set_overlap_sp1:
        gene_list_sp1 = set_overlap_sp1[sp2_cell_set]
        gene_list_sp2 = set_overlap_sp2[sp2_cell_set]

        count = 1
        for sp2_cell in sp2_cell_set:
            # Overlap with next if they are in this same sp2_cell_set
            if count >= len(sp2_cell_set):
                overlap_with_next = 0
            else:
                overlap_with_next = len(gene_list_sp1)
            
            this_ribbon = [
                sp1_cell,    # top rect label
                sp2_cell,    # bottom rect label
                len(gene_list_sp1), # top ribbon width
                len(gene_list_sp2), # bottom ribbon width
                overlap_with_next
            ]
            if len(gene_list_sp1)!=0 or len(gene_list_sp2)!=0:
                ribbon_widths.append(this_ribbon)
            count += 1
    
    return sp1_rects_dict, sp2_rects_dict, ribbon_widths


# ============================================================================
# Calculate Jaccard similarity score based on proportion of shared genes

jaccard_hs = pd.DataFrame(
    index=[f'Hs_{x}' for x in adata_hs.obs['cell_type_coarse'].cat.categories],
    columns=[f'Cr_{x}' for x in adata_cr.obs['cell_type'].cat.categories]
)
jaccard_cr = pd.DataFrame(
    index=[f'Hs_{x}' for x in adata_hs.obs['cell_type_coarse'].cat.categories],
    columns=[f'Cr_{x}' for x in adata_cr.obs['cell_type'].cat.categories]
)

for hs_cell in jaccard_hs.index:
    for cr_cell in jaccard_hs.columns:
        hs_cell2 = hs_cell.split('_')[1]
        cr_cell2 = cr_cell.split('_')[1]
        hs_genes, cr_genes = get_gene_dict_sets(
            hs_cell2, [cr_cell2], gene_orthology='orthofinder')
        
        count_hs = len(hs_genes[hs_cell2][f'all'])
        count_cr = len(cr_genes[cr_cell2][f'all'])
        count_shared_hs = len(hs_genes[hs_cell2][f'shared_w_{cr_cell2}'])
        count_shared_cr = len(cr_genes[cr_cell2][f'shared_w_{hs_cell2}'])
        jaccard_hs.loc[hs_cell, cr_cell] = (
            count_shared_hs / (count_hs + count_cr - count_shared_hs))
        jaccard_cr.loc[hs_cell, cr_cell] = (
            count_shared_cr / (count_hs + count_cr - count_shared_cr))

# Overall score is the average of Jaccard index in both directions
jaccard = (jaccard_cr + jaccard_hs) / 2
jaccard = jaccard.astype(float)

# ============================================================================
# Plot for human cell

fig_height = 1.25
n = 5

def plot_for_hs_cell(hs_cell, fig_width=2.5, ortholog='orthofinder',
                     out_path=out_path):
    """
    Function to make and save a GS plot for a given human cell, selecting the
    other species's cell states which are most similar (based on Jaccard
    similarity score) to plot against it.
    """

    # --------------------------------
    # Get Ciona cell types to plot
    cr_cell_list = set(jaccard.columns)
    cr_cell_list = [x.replace('Cr_', '') for x in cr_cell_list]

    # Sort Ciona cell types by all-gene similarity
    order = np.argsort(jaccard.loc['Hs_' + hs_cell,
                                       ['Cr_'+x for x in cr_cell_list]])[::-1]
    cr_cell_list = list(np.array(cr_cell_list)[order.values])

    # Keep the top n
    cr_cell_list = cr_cell_list[:n]

    if len(cr_cell_list) == 0:
        print(f'{hs_cell}: No similar cell types')
        return None

    # --------------------------------
    # Get numbers for plot
    hs_genes, cr_genes = get_gene_dict_sets(
        hs_cell, cr_cell_list, gene_orthology=ortholog)
    hs_rects_dict, cr_rects_dict, ribbon_widths = get_sankey_plot_numbers(
        hs_genes, cr_genes)
    hs_rects_unmatched = {x: len(degs_hs[x]) - hs_rects_dict[x]
                          for x in hs_rects_dict}
    cr_rects_unmatched = {x: len(degs_cr[x]) - cr_rects_dict[x]
                          for x in cr_rects_dict}
    
    # Filter to only include cell types with >N genes
    cr_rects_dict = {x: cr_rects_dict[x] for x in cr_rects_dict
                     if (cr_rects_dict[x] > 10) and (x in [y[1] for y in ribbon_widths])}
    ribbon_widths = [x for x in ribbon_widths if x[1] in cr_rects_dict]

    # Rename cell types for plotting
    hs_rects_dict = {x.replace(' ', '\n'): hs_rects_dict[x] for x in hs_rects_dict}
    cr_rects_dict = {x.replace(' ', '\n'): cr_rects_dict[x] for x in cr_rects_dict}
    hs_rects_unmatched = {x.replace(' ', '\n'): hs_rects_unmatched[x] for x in hs_rects_unmatched}
    cr_rects_unmatched = {x.replace(' ', '\n'): cr_rects_unmatched[x] for x in cr_rects_unmatched}
    ribbon_widths = [[x[0].replace(' ', '\n')] + [x[1].replace(' ', '\n')] + x[2:]
                     for x in ribbon_widths]

    # Plot
    plot_modified_sankey(hs_rects_dict, cr_rects_dict, ribbon_widths,
                         top_rect_unmatched=hs_rects_unmatched,
                         bot_rect_unmatched=cr_rects_unmatched,
                         rect_height=0.15, vertical_rect_spacing=1,
                         scalebar_size=100, figsize=(fig_width, fig_height),
                         color_top=sp2_color, color_bot=sp1_color)
    plt.savefig(os.path.join(out_path, f'{hs_cell}.pdf'))
    plt.close()


# ------------------------------------
# Make and save plots
out_path2 = os.path.join(out_path, 'plots_human_top')
if not os.path.exists(out_path2): os.mkdir(out_path2)
plot_for_hs_cell('hematopoietic progenitor', fig_width=3, ortholog='orthofinder', out_path=out_path2)
plot_for_hs_cell('neutrophil', fig_width=3, ortholog='orthofinder', out_path=out_path2)
plot_for_hs_cell('monocyte', fig_width=3, ortholog='orthofinder', out_path=out_path2)
plot_for_hs_cell('macrophage', fig_width=3, ortholog='orthofinder', out_path=out_path2)
plot_for_hs_cell('NK cell', fig_width=3, ortholog='orthofinder', out_path=out_path2)
plot_for_hs_cell('erythrocyte', fig_width=3, ortholog='orthofinder', out_path=out_path2)
plot_for_hs_cell('T cell', fig_width=3, ortholog='orthofinder', out_path=out_path2)


# ============================================================================
# Plot for ciona cell

fig_height = 1.25
n = 5

def plot_for_cr_cell(cr_cell, fig_width=2.5, ortholog='orthofinder',
                     out_path=out_path):
    """
    Function to make and save a GS plot for a given cell state in the non-
    human species, selecting the human cell states which are most similar
    (based on Jaccard similarity score) to plot against it.
    """

    # Flip similarity scores
    this_jaccard = jaccard.T

    # Get human cell types to plot
    hs_cell_list = set(this_jaccard.columns)
    hs_cell_list = [x.replace('Hs_', '') for x in hs_cell_list]

    # Sort human cell types by all-gene similarity
    order = np.argsort(this_jaccard.loc['Cr_' + cr_cell,
                                       ['Hs_'+x for x in hs_cell_list]])[::-1]
    hs_cell_list = list(np.array(hs_cell_list)[order.values])

    # Keep top n
    hs_cell_list = hs_cell_list[:n]

    if len(hs_cell_list) == 0:
        print(f'{hs_cell}: No similar cell types')
        return None

    # --------------------------------
    # Get numbers for plot
    hs_genes, cr_genes = get_gene_dict_sets(
        cr_cell, hs_cell_list, sp1='ciona',
        gene_orthology=ortholog)
    cr_rects_dict, hs_rects_dict, ribbon_widths = get_sankey_plot_numbers(
        cr_genes, hs_genes)
    hs_rects_unmatched = {x: len(degs_hs[x]) - hs_rects_dict[x]
                          for x in hs_rects_dict}
    cr_rects_unmatched = {x: len(degs_cr[x]) - cr_rects_dict[x]
                          for x in cr_rects_dict}
    
    # Filter to only include cell types with >N genes
    hs_rects_dict = {x: hs_rects_dict[x] for x in hs_rects_dict
                     if (hs_rects_dict[x] > 10) and (x in [y[1] for y in ribbon_widths])}
    ribbon_widths = [x for x in ribbon_widths if x[1] in hs_rects_dict]

    # Rename cell types for plotting
    hs_rects_dict = {x.replace(' ', '\n'): hs_rects_dict[x] for x in hs_rects_dict}
    # cr_rects_dict = {x.replace(' ', '\n'): cr_rects_dict[x] for x in cr_rects_dict}
    hs_rects_unmatched = {x.replace(' ', '\n'): hs_rects_unmatched[x] for x in hs_rects_unmatched}
    # cr_rects_unmatched = {x.replace(' ', '\n'): cr_rects_unmatched[x] for x in cr_rects_unmatched}
    ribbon_widths = [[x[0]] + [x[1].replace(' ', '\n')] + x[2:]
                     for x in ribbon_widths]

    # Plot
    plot_modified_sankey(cr_rects_dict, hs_rects_dict, ribbon_widths,
                         top_rect_unmatched=cr_rects_unmatched,
                         bot_rect_unmatched=hs_rects_unmatched,
                         rect_height=0.15, vertical_rect_spacing=1,
                         scalebar_size=100, figsize=(fig_width, fig_height),
                         color_top=sp1_color, color_bot=sp2_color)
    plt.savefig(os.path.join(out_path, f'{cr_cell.replace("/","-")}.pdf'))
    plt.close()


# ------------------------------------
# Make and save plots
out_path2 = os.path.join(out_path, 'plots_ciona_top')
if not os.path.exists(out_path2): os.mkdir(out_path2)
plot_for_cr_cell('GA', fig_width=3, ortholog='orthofinder', out_path=out_path2)
plot_for_cr_cell('HA-1 (phag.)', fig_width=3, ortholog='orthofinder', out_path=out_path2)
plot_for_cr_cell('HA-2 (phag.)', fig_width=3, ortholog='orthofinder', out_path=out_path2)
plot_for_cr_cell('HA-3', fig_width=3, ortholog='orthofinder', out_path=out_path2)
plot_for_cr_cell('HA-4', fig_width=3, ortholog='orthofinder', out_path=out_path2)
plot_for_cr_cell('HA-5', fig_width=3, ortholog='orthofinder', out_path=out_path2)
plot_for_cr_cell('BLC', fig_width=3, ortholog='orthofinder', out_path=out_path2)
plot_for_cr_cell('URG-1 (large)', fig_width=3, ortholog='orthofinder', out_path=out_path2)
plot_for_cr_cell('LGH/MC', fig_width=3, ortholog='orthofinder', out_path=out_path2)
plot_for_cr_cell('cMPP', fig_width=3, ortholog='orthofinder', out_path=out_path2)


# ============================================================================
# PLOT JACCARD SIMILARITY IN HEATMAP

jaccard.to_csv(os.path.join(out_path, f'jaccard_similarity_matrix.csv'))
rows_all_zeros = jaccard.index[jaccard.sum(axis=1) == 0]
cols_all_zeros = jaccard.columns[jaccard.sum(axis=0) == 0]

# Plot as ordered from high to low similarity score for Ciona
marginals = jaccard.max(axis=0)
sp1_cell_type_order = jaccard.columns[np.argsort(marginals)[::-1]]
sp2_max_match_to_sp1 = jaccard.idxmax(axis=1)
sp2_cell_type_order = []
for sp1_cell in sp1_cell_type_order:
    these_sp2_cells = sp2_max_match_to_sp1.index[
        np.where(sp2_max_match_to_sp1 == sp1_cell)[0]]
    marginals = jaccard.loc[these_sp2_cells].max(axis=1)
    these_sp2_cells_ordered = these_sp2_cells[np.argsort(marginals)[::-1]]
    sp2_cell_type_order += list(these_sp2_cells_ordered)
f = plt.figure(figsize=(2.7, 4.7))
ax = plt.subplot(1, 1, 1)
sns.heatmap(
    jaccard.loc[sp2_cell_type_order, sp1_cell_type_order].T,
    cmap='rocket_r',
    xticklabels=True,
    yticklabels=True,
    vmin=0, vmax=0.12,
)
ax.xaxis.set_ticks_position('top')  # Move ticks to top
ax.xaxis.set_label_position('top')  # Move label position to top
ax.set_xticklabels([x.split('_')[-1] for x in sp2_cell_type_order], rotation=90)
ax.set_yticklabels([x.split('_')[-1] for x in sp1_cell_type_order])
plt.tight_layout()
plt.savefig(os.path.join(out_path, f'jaccard_similarity.pdf'))
plt.close()

# Plot transposed matrix
marginals = jaccard.max(axis=1)
sp1_cell_type_order = jaccard.index[np.argsort(marginals)[::-1]]
sp2_max_match_to_sp1 = jaccard.idxmax(axis=0)
sp2_cell_type_order = []
for sp1_cell in sp1_cell_type_order:
    these_sp2_cells = sp2_max_match_to_sp1.index[
        np.where(sp2_max_match_to_sp1 == sp1_cell)[0]]
    marginals = jaccard[these_sp2_cells].max(axis=0)
    these_sp2_cells_ordered = these_sp2_cells[np.argsort(marginals)[::-1]]
    sp2_cell_type_order += list(these_sp2_cells_ordered)
f = plt.figure(figsize=(5, 2.3))
ax = plt.subplot(1, 1, 1)
sns.heatmap(
    jaccard.loc[sp1_cell_type_order, sp2_cell_type_order],
    cmap='rocket_r',
    xticklabels=True,
    yticklabels=True,
    vmin=0, vmax=0.12,
)
ax.xaxis.set_ticks_position('top')  # Move ticks to top
ax.xaxis.set_label_position('top')  # Move label position to top
ax.set_yticklabels([x.split('_')[-1] for x in sp1_cell_type_order])
ax.set_xticklabels([x.split('_')[-1] for x in sp2_cell_type_order], rotation=90)
plt.tight_layout()
plt.savefig(os.path.join(out_path, f'jaccard_similarity.T.pdf'))
plt.close()
