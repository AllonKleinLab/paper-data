import sys
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
from tqdm import tqdm
from matplotlib.colors import ListedColormap

# ============================================================================
# Colors & colormaps


def hex2rgb(color):
    """
    Convert hexadecimal string to RGB list representation of a color.
    """
    # Remove '#' if at the start of color
    color = color.split('#')[-1]

    if len(color) == 6:
        # Parse colors into R, G, B values from 0 to 256
        rgb_hex = [color[:2], color[2:4], color[4:6]]
        rgb = [int(i, 16) for i in rgb_hex]
        return rgb
    else:
        print('Wrong formatting for hexadecimal color.')


def cmap(color_list):
    """
    Create a custom cmap, gradient between two given colors. Colors are
    entered as hexadecimal strings.
    """
    # Catch incorrect inputs
    if not isinstance(color_list, list):
        print('Expected list input.')
        return
    elif len(color_list) < 2:
        print('Expected list of length 2 or greater.')
        return

    # Parse colors into R, G, B values from 0 to 255
    rgb_list = [hex2rgb(c) for c in color_list]

    # Build color map between list of colors
    N = 256
    n = int(np.ceil(N / (len(rgb_list) - 1)))
    vals = np.ones((N, 4))
    for i in range(len(rgb_list) - 1):
        # Indices to set in vals
        low_i = i * n
        high_i = (i + 1) * n
        # If last color and there is a remainder when dividing N
        if high_i > N:
            n = n - (high_i - N)
            high_i = N
            
        # Linspace of colors
        rgb1 = rgb_list[i]
        rgb2 = rgb_list[i+1]
        vals[low_i : high_i, 0] = np.linspace(rgb1[0]/256, rgb2[0]/256, n)
        vals[low_i : high_i, 1] = np.linspace(rgb1[1]/256, rgb2[1]/256, n)
        vals[low_i : high_i, 2] = np.linspace(rgb1[2]/256, rgb2[2]/256, n)

    return ListedColormap(vals)


pink2blue = ['#ab003f', '#ec6697', '#f2f2f2', '#637ed6', '#334faa']
cmap_pink2blue = cmap(pink2blue)
cmap_pink2blue_r = cmap(pink2blue[::-1])

# ============================================================================
# UMAP plotting


def pl_umap_separate(adata, color:str, title=None, palette=None, s=None,
            alpha=None, ncol=4, width_scale=5, height_scale=5):
    """
    Plots categorical information on separate UMAPs rather than together on
    a single UMAP. Makes it easier to see where there is overlap.
    """
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
# Density plot


def neighborhood_composition(adata, obs, num_neighbors=200, use_rep='X_pca',
                             use_existing_neighbor_graph=False):
    """
    This function returns a pandas DataFrame with the number of neighbors in
    a given obs category.
    
    INPUTS:
    adata           scanpy AnnData object, the scRNA-seq data
    obs             string, name of column in adata.obs with categorical
                    values
    num_neighbors   int, the size of the neighborhood for counting how many
                    neighbors are from each obs value. Default = 200.
    cmap            matplotlib cmap or string cmap name, the cmap used for
                    plotting the calculated density.

    OUTPUTS:
    f               matplotlib figure object, containing the plotted densities
    Note - this function also adds obs to adata with the calculated densities:
    adata.obs['{val1}_over_{val2}_density']
    """

    # ------------------------------------
    # 1. Get N nearest neighbors
    if not use_existing_neighbor_graph:
        N = num_neighbors
        print(f'Getting {N} neighbors...')
        sc.pp.neighbors(adata, n_neighbors=N, use_rep=use_rep)
    # cells, neighbors = np.where(adata.uns['neighbors']['distances'].toarray()
    #                             > 0)
    cells, neighbors = np.where(adata.obsp['distances'].toarray() > 0)
    # cell[i] and neighbors[i] are the i^th pair of neighbors

    # ------------------------------------
    # 2. For each cell, get number of neighbors in each sample/condition
    print(f'Counting # of cells from each sample in neighborhood...')

    # Initialize
    neighborhood_df = pd.DataFrame(
        np.zeros([adata.shape[0], len(np.unique(adata.obs[obs]))]),
        index=adata.obs.index,
        columns=np.unique(adata.obs[obs])
    )
    neighborhood_df = neighborhood_df.astype(int)
    
    # Loop through cells, count # neighbors for each condition
    for iCell in tqdm(range(adata.shape[0])):
        # Get neighbors + their sample for this cell
        these_neighbors = neighbors[np.where(cells == iCell)[0]]
        neighbor_vals = adata.obs.loc[adata.obs.index[these_neighbors], obs]
        
        # Save number of neighbors for each sample
        for cond in neighborhood_df.columns:
            neighborhood_df.loc[adata.obs.index[iCell], cond] = \
                np.sum(neighbor_vals == cond)

    return neighborhood_df


def get_nrow_ncol(nplots, ncol=4):
    """
    For a given number of plots in a subplot figure, returns the number of
    rows and columns required to fit that many plots. Can adjust the number of
    columns.
    """
    ncol = min(ncol, nplots)
    nrow = int(np.ceil(nplots / ncol))

    return (nrow, ncol)


def sample_density(adata, obs, obs_vals, num_neighbors=200, use_rep='X_pca',
                   cmap=cmap_pink2blue, sort_order=True,
                   log_offset=10**-10):
    """
    This function calculates the relative density of two values of a
    categorical observation in scRNA-seq data, then plots that density on a
    UMAP. For a given cell, the relative density is:
        log2( (# cells in neighborhood with value1) /
              (# cells in neighborhood with value2) )
    
    If more than one pair of values is given, a subplot showing all
    pairs' relative density is returned.

    INPUTS:
    adata           scanpy AnnData object, the scRNA-seq data
    obs             string, name of column in adata.obs with categorical
                    values
    obs_vals        list of tuples, names of pairs of values in adata.obs to
                    compare densities. E.g. [(val1, val2)] will give a density
                    plot of:
                    (# cells in a neighborhood with adata.obs[obs] == val1) /
                    (# cells in a neighborhood with adata.obs[obs] == val2)
    num_neighbors   int, the size of the neighborhood for counting how many
                    neighbors are from each obs value. Default = 200.
    cmap            matplotlib cmap or string cmap name, the cmap used for
                    plotting the calculated density.

    OUTPUTS:
    f               matplotlib figure object, containing the plotted densities
    Note - this function also adds obs to adata with the calculated densities:
    adata.obs['{val1}_over_{val2}_density']
    """

    # Values listed in obs_vals must be in adata.obs[obs]
    for (val1, val2) in obs_vals:
        if (val1 not in np.unique(adata.obs[obs])):
            print(f'{val1} is not a value in adata.obs[{obs}]')
            return
        elif (val2 not in np.unique(adata.obs[obs])):
            print(f'{val2} is not a value in adata.obs[{obs}]')
            return

    # ------------------------------------
    # 1. For each cell, get number of neighbors in each sample/condition

    neighborhood_df = neighborhood_composition(adata, obs, use_rep=use_rep,
                                               num_neighbors=num_neighbors)

    # ------------------------------------
    # 2. Normalize

    # Add sample/condition information for cell itself
    neighborhood_df += adata.obs[obs].str.get_dummies()

    # Normalize so that each row adds to 1
    # neighborhood_df = neighborhood_df / num_neighbors

    # ------------------------------------
    # 3. Calculate & plot fold change between each pair of values in obs_vals

    nrow, ncol = get_nrow_ncol(len(obs_vals))
    f, ax = plt.subplots(nrow, ncol, figsize=(ncol*5, nrow*5.1))

    for i, (val1, val2) in enumerate(obs_vals):

        # Get fold change between conditions
        adata.obs[f'{val1}_over_{val2}_density'] = np.log2(
            (neighborhood_df.loc[:, val1] + log_offset)
            / (neighborhood_df.loc[:, val2] + log_offset)
        )
        
        # Get max absolute value of 99th and 1st percentile
        vmax = max(
            abs(np.percentile(adata.obs[f'{val1}_over_{val2}_density'], 99)),
            abs(np.percentile(adata.obs[f'{val1}_over_{val2}_density'], 1))
        )
        
        # Plot
        if nrow == 1 and ncol == 1: iAx = ax
        elif nrow == 1 and ncol > 1: iAx = ax[i]
        else: iAx = ax.flatten()[i]
        sc.pl.umap(adata, color=f'{val1}_over_{val2}_density',
                   s=30, show=False, ax=iAx, cmap=cmap, vmax=vmax, vmin=-vmax,
                   title=('$\log_2$({val1} / {val2})\nin neighborhood '
                          + '(N={num_neighbors})'),
                   sort_order=sort_order)
        plt.tight_layout()
    
    return f
