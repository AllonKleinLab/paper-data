import sys
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
from tqdm import tqdm

from plotting_helper_functions import *
from general_helper_functions import *
from highly_variable_genes import *

# ============================================================================
# Helper functions for UMAP plots

def pl_umap(adata, color=None, title=None, groups=None, gene_symbols=None,
            legend_loc='right margin', legend_fontoutline=None, layer=None,
            vmin=None, vmax=None, vcenter=None, categorical_rand_order=None,
            cmap=None, palette=None, s=None, alpha=None, na_color='lightgray',
            ncol=4, width_scale=5, height_scale=5, spines=True):
    """
    Calls sc.pl.umap, but outputs saved files with neater formatting.
    """
    # Initialize plot
    if not isinstance(color, list): color = [color]
    nrow, ncol = get_nrow_ncol(len(color), ncol=ncol)
    f, ax = plt.subplots(nrow, ncol, figsize=(width_scale * ncol,
                                              height_scale * nrow))

    # Make each subplot in color list
    for i, c in enumerate(color):
        # Format inputs
        iTitle = title[i] if isinstance(title, list) else title
        # iGroups = groups[i] if isinstance(groups, list) else groups
        iGene_symbols = gene_symbols[i] if isinstance(gene_symbols, list) \
            else gene_symbols
        iVmin = vmin[i] if isinstance(vmin, list) else vmin
        iVmax = vmax[i] if isinstance(vmax, list) else vmax
        iVcenter = vcenter[i] if isinstance(vcenter, list) else vcenter
        iS = s[i] if isinstance(s, list) else s
        iAlpha = alpha[i] if isinstance(alpha, list) else alpha
        iLegend_loc = legend_loc[i] if isinstance(legend_loc, list) else \
            legend_loc
        iLegend_fontoutline = legend_fontoutline[i] if \
            isinstance(legend_fontoutline, list) else legend_fontoutline
        iCmap = cmap[i] if isinstance(cmap, list) else cmap
        iCategorical_rand_order = (
            categorical_rand_order[i]
            if isinstance(categorical_rand_order, list)
            else categorical_rand_order)

        # Check if palette is a list of palettes (list of lists)
        if isinstance(palette, list):
            iPalette = (palette[i] if any(isinstance(l, list) for l in palette)
                        else palette)
        else: iPalette = palette
        iNa_color = na_color[i] if isinstance(na_color, list) else na_color
        if nrow == 1 and ncol == 1: iAx = ax
        elif nrow == 1 and ncol > 1: iAx = ax[i]
        else: iAx = ax.flatten()[i]
        
        # Make subplot
        if not iCategorical_rand_order == True:
            sc.pl.umap(adata, color=c, show=False, ax=iAx,
                groups=groups, title=iTitle, legend_loc=iLegend_loc,
                legend_fontoutline=iLegend_fontoutline,
                vmin=iVmin, vmax=iVmax, vcenter=iVcenter,
                s=iS, alpha=iAlpha, cmap=iCmap, palette=iPalette,
                na_color=iNa_color, gene_symbols=iGene_symbols,
                layer=layer)
        
        # For randomizing z-order of categorical labels, use plt.scatter
        elif iCategorical_rand_order == True:
            umap_x = adata.obsm['X_umap'][:, 0].copy()
            umap_y = adata.obsm['X_umap'][:, 1].copy()

            color_values = adata.obs[c].copy()

            # Set color palette
            if iPalette is None:
                iPalette = adata.uns[f'{c}_colors'].copy()
            color_values = color_values.cat.rename_categories(
                np.arange(len(color_values.cat.categories)))
            color_for_scatter = np.array(
                [iPalette[x] for x in color_values.values])
            
            # Set point size
            if iS is None:
                # Use same default as sc.pl.umap
                iS = 120000 / len(color_values)
            iS = iS / 4

            # Randomize plotting order
            index_order = np.arange(len(color_for_scatter))
            np.random.shuffle(index_order)

            # Make UMAP
            plt.scatter(umap_x[index_order], umap_y[index_order],
                        c=color_for_scatter[index_order],
                        s=iS, linewidths=0, alpha=iAlpha)
            
            # Label axes etc.
            plt.xlabel('UMAP1')
            plt.ylabel('UMAP2')
            if iTitle is None: plt.title(c)
            else: plt.title(iTitle)
            plt.tick_params(left=False, right=False, labelleft=False, 
                            labelbottom=False, bottom=False)
        
        if not spines: iAx.axis('off')
    
    # Return
    plt.tight_layout()
    return f


def pl_umap_separate(adata, color:str, title=None, palette=None, s=None,
            alpha=None, ncol=4, width_scale=5, height_scale=5):
    obs_cats = np.unique(adata.obs[color])
    ncol = min(4, len(obs_cats))
    nrow = int(np.ceil(len(obs_cats) / ncol))
    f, ax = plt.subplots(nrow, ncol,
                         figsize=(width_scale*ncol, height_scale*nrow))

    for i, c in enumerate(obs_cats):
        # Format inputs
        iTitle = title[i] if isinstance(title, list) else c
        iS = s[i] if isinstance(s, list) else s
        iAlpha = alpha[i] if isinstance(alpha, list) else alpha
        if nrow == 1 and ncol == 1: iAx = ax
        elif nrow == 1 and ncol > 1: iAx = ax[i]
        else: iAx = ax.flatten()[i]
        iPalette = palette

        sc.pl.umap(adata, color=color, groups=c, s=iS, title=iTitle,
                   ax=iAx, show=False, palette=iPalette, alpha=iAlpha)
        iAx.legend().remove()
    
    f.tight_layout()
    return f


# ============================================================================
# Other helper functions for scanpy


def get_vmax(gene_list, adata, q=95):
    """
    !! Most of the time for sc.pl.umap you should use vmax='pXX', e.g. 'p99'
    for the 99th percentile - see scanpy docs.

    For a single gene or list of genes, return the 90th percentile of the
    non-zero values. Used for setting vmax in sc.pl.umap()
    """
    if isinstance(gene_list, str):
        gene_list = [gene_list]

    # Initialize output
    output_list = []

    for gene in gene_list:
        this_row = adata.raw[:, adata.raw.var_names.isin([gene])].X.todense()
        if np.all(this_row==0):
            output_list.append(0)
        else:
            output_list.append(np.percentile(this_row[this_row!=0], q))
    
    return output_list


def cluster_centroids(adata, groupby, clusters=None, use_raw=True,
    gene_list=None, layer=None):
    """
    Return the mean expression across clusters

    Inputs
        adata : AnnData object, single cell dataset
        groupby : The name of the adata.obs variable to group cells by before
            taking the mean expression level.
        clusters : list, subset of clusters to return data for
    
    Outputs
        A pandas DataFrame containing mean expression of each gene
            across each group of the groupby variable
    """
    # If no list of clusters given, get from adata.obs[groupby]
    if clusters is None:
        clusters = adata.obs[groupby].unique().categories
    num_clusters = len(clusters)
    if gene_list is None: num_genes = len(adata.var.index)
    else: num_genes = len(gene_list)

    # Get adata.X and number of genes
    if gene_list is None: gene_list = adata.var.index.tolist()
    if use_raw:
        gene_bool = [(i in gene_list) for i in adata.raw.var.index]
        # gene_list_ordered = [i for i in adata.var.index if i in gene_list]
        X = adata.raw.X[:, gene_bool]
        gene_list_ordered = list(adata.raw.var.index[gene_bool])
    elif layer is not None:
        gene_bool = [(i in gene_list) for i in adata.raw.var.index]
        X = adata.layers[layer]
        gene_list_ordered = list(adata.raw.var.index[gene_bool])
    else:
        gene_bool = [(i in gene_list) for i in adata.var.index]
        # gene_list_ordered = [i for i in adata.var.index if i in gene_list]
        X = adata.X[:, gene_bool]
        gene_list_ordered = list(adata.var.index[gene_bool])
        
        # Error if gene_list not all in X
        X_gene_list = np.sort(adata[:, gene_bool].var.index.tolist())
        actual_gene_list = np.sort(list(gene_list))
        try: list_match = np.all(X_gene_list == actual_gene_list)
        except: list_match = False
        if list_match == False:
            print('gene_list does not match genes in adata.X')
            return

    # Initiate array to save means
    exp_mean = np.zeros((num_clusters, num_genes))

    for i, cl in enumerate(clusters):
        cluster_X = X[adata.obs[groupby] == cl]

        # Save mean of each gene
        exp_mean[i, :] = cluster_X.mean(axis=0)

    return pd.DataFrame(exp_mean, columns=gene_list_ordered, index=clusters)
