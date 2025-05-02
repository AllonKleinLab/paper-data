import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
from matplotlib.colors import ListedColormap
from scipy.cluster.hierarchy import dendrogram, linkage

# ============================================================================
# Plotting helper functions

def get_nrow_ncol(nplots, ncol=4):
    """
    For a given number of plots in a subplot figure, returns the number of
    rows and columns required to fit that many plots. Can adjust the number of
    columns.
    """
    ncol = min(ncol, nplots)
    nrow = int(np.ceil(nplots / ncol))

    return (nrow, ncol)


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


def rgb2hex(color, max_value=255):
    """
    Convert hexadecimal string to RGB list representation of a color.
    """
    if len(color) == 3:
        # Parse colors into R, G, B values from 0 to 256
        rgb_hex = '#'
        for c in color:
            # If on scale from 0 to 1, scale to 0 to 255 as int
            if max_value != 255:
                c = c * (255 / max_value)
            if not isinstance(c, int): c = int(c)

            # Add two digit hexadecimal
            c_hex = hex(c).split('0x')[-1]
            if len(c_hex) == 1: c_hex = '0' + c_hex
            rgb_hex += c_hex
        return rgb_hex
    else:
        print('Wrong formatting for RGB color.')


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


def add_to_cmap_end(existing_cmap, color, replace='vmin'):
    """
    Add color to the end of existing_cmap such that the vmin value of the
    original cmap is replaced by color.

    To instead add the color such that vmax is replaced, set replace='vmax'
    """
    # Parse colors into R, G, B values from 0 to 255
    color_rgba = np.array(hex2rgb(color) + [255]) / 255
    color_rgba = color_rgba.reshape(1, 4)

    # Parse cmap if string input given
    try:
        cmap_obj = sns.color_palette(existing_cmap, as_cmap=True)
        cmap_obj_list = [cmap_obj(x) for x in range(cmap_obj.N)]
        existing_cmap_arr = np.array(cmap_obj_list)
    except:
        if existing_cmap == 'Reds':
            existing_cmap_arr = plt.cm.Reds(np.linspace(0, 1, 256))
        else:
            existing_cmap_arr = existing_cmap

    # Add color to cmap
    new_cmap_arr = existing_cmap_arr.copy()
    if replace == 'vmin':
        new_cmap_arr[0, :] = color_rgba
    elif replace == 'vmax':
        new_cmap_arr[-1, :] = color_rgba
    
    return ListedColormap(new_cmap_arr)


# ============================================================================
# Colors

# Colors and palettes
purple = '#A67CC8'
purple_text = '#9A72B9'
palette1 = [
    '#3e5dc4',
    '#daa53b',
    '#C72652',
    '#53C0A7',
    '#B778D5',
]
# Copied from scanpy pallete on github:
# https://github.com/scverse/scanpy/blob/main/src/scanpy/plotting/palettes.py
palette_102 = [
    # "#000000",  # remove the black, as often, we have black colored annotation
    "#FFFF00",
    "#1CE6FF",
    "#FF34FF",
    "#FF4A46",
    "#008941",
    "#006FA6",
    "#A30059",
    "#FFDBE5",
    "#7A4900",
    "#0000A6",
    "#63FFAC",
    "#B79762",
    "#004D43",
    "#8FB0FF",
    "#997D87",
    "#5A0007",
    "#809693",
    "#6A3A4C",
    "#1B4400",
    "#4FC601",
    "#3B5DFF",
    "#4A3B53",
    "#FF2F80",
    "#61615A",
    "#BA0900",
    "#6B7900",
    "#00C2A0",
    "#FFAA92",
    "#FF90C9",
    "#B903AA",
    "#D16100",
    "#DDEFFF",
    "#000035",
    "#7B4F4B",
    "#A1C299",
    "#300018",
    "#0AA6D8",
    "#013349",
    "#00846F",
    "#372101",
    "#FFB500",
    "#C2FFED",
    "#A079BF",
    "#CC0744",
    "#C0B9B2",
    "#C2FF99",
    "#001E09",
    "#00489C",
    "#6F0062",
    "#0CBD66",
    "#EEC3FF",
    "#456D75",
    "#B77B68",
    "#7A87A1",
    "#788D66",
    "#885578",
    "#FAD09F",
    "#FF8A9A",
    "#D157A0",
    "#BEC459",
    "#456648",
    "#0086ED",
    "#886F4C",
    "#34362D",
    "#B4A8BD",
    "#00A6AA",
    "#452C2C",
    "#636375",
    "#A3C8C9",
    "#FF913F",
    "#938A81",
    "#575329",
    "#00FECF",
    "#B05B6F",
    "#8CD0FF",
    "#3B9700",
    "#04F757",
    "#C8A1A1",
    "#1E6E00",
    "#7900D7",
    "#A77500",
    "#6367A9",
    "#A05837",
    "#6B002C",
    "#772600",
    "#D790FF",
    "#9B9700",
    "#549E79",
    "#FFF69F",
    "#201625",
    "#72418F",
    "#BC23FF",
    "#99ADC0",
    "#3A2465",
    "#922329",
    "#5B4534",
    "#FDE8DC",
    "#404E55",
    "#0089A3",
    "#CB7E98",
    "#A4E804",
    "#324E72",
]

# ============================================================================
# MISC PLOTTING

def hierarchical_clustering(matrix, ax=None, no_plot=False,
                            method='single', metric='cosine',
                            labels=None, orientation='left',
                            optimal_ordering=False, count_sort='descending',
                            leaf_font_size=None, color='#888888'):
    linkage_data = linkage(matrix, method=method, metric=metric,
                           optimal_ordering=optimal_ordering)
    dendr = dendrogram(linkage_data, labels=labels, no_plot=no_plot,
                       count_sort=count_sort, ax=ax,
                       orientation=orientation, leaf_font_size=leaf_font_size,
                       link_color_func=lambda _: color)
    if ax is not None:
        ax.spines[['left', 'bottom', 'top', 'right']].set_visible(False)
        ax.set_xticks([])
    return dendr
