import os
import numpy as np
from tqdm import tqdm
from skimage import io
from matplotlib_scalebar.scalebar import ScaleBar

from general_helper_functions import *
from plotting_helper_functions import *

# ============================================================================
# Plotting helper functions

# Color constants
# standards: https://www.webucator.com/article/python-color-constants-module/
color_constants = {
    'red': '#ff0000',
    'green': '#00ff00',
    'blue': '#0000ff',
    'magenta': '#ff00ff',
    'yellow': '#ffff00',
    'cyan': '#00ffff',
    'white': '#ffffff',
    'black': '#000000',
    'gray': '#808080',
}

# ============================================================================
# Opening images in napari

def import_image_large_scan(path_to_images:str, x_range=None, y_range=None,
                            z_range=None):
    """
    INPUT
    path_to_images : string, path to folder containing tiled images
    x_range : list of length 2 (default None), range (inclusive) of
        x-positions to include
    y_range : list of length 2 (default None), range (inclusive) of
        y-positions to include
    y_range : list of length 2 (default None), range (inclusive) of
        z-slices to include

    OUTPUT
    full_img : mutidimensional numpy array, output image. Dimensions are
        [xy_positions] x [z_slices] x [channels] x [img_dim1] x [img_dim2]
    xy_pos : numpy array, xy positions corresponding to full_img xy_positions
        dimension. E.g. xy_pos[0, :] gives the x and y position corresponding
        to full_img[0, ...]
    """
    # Check input format is correct
    for i, range in enumerate([x_range, y_range, z_range]):
        if isinstance(range, list):
            if len(range) != 2:
                raise Exception(f'{["x","y","z"][i]}_range must be length 2')
        # If not a list and also not None
        elif range is not None:
            raise Exception(f'{["x","y","z"][i]}_range must be a list')

    print(f'Importing from {path_to_images}')
    file_list = os.listdir(path_to_images)
    tmp_img_dict = {}
    for file in tqdm(np.sort(file_list)):
        if (file[:5] == 'tile_') and (file[-4:] == '.tif'):
            # Parse information from file name
            file_info = file.split('.')[0].split('_')
            if file_info[1][0] == 'x': x = int(file_info[1][1:])
            else: raise Exception('Incorrect file name formatting - no x')
            if file_info[2][0] == 'y': y = int(file_info[2][1:])
            else: raise Exception('Incorrect file name formatting - no y')
            if '_z' in file:
                has_z = True
                if file_info[3][0] == 'z': z = int(file_info[3][1:])
                else: raise Exception('Incorrect file name formatting - no z')
            else: has_z = False

            if has_z:
                # Check if in specificied ranges
                in_range_bool = (
                    ((x_range is None) or (x>=x_range[0] and x<=x_range[-1]))
                    and ((y_range is None) or (y>=y_range[0] and y<=y_range[-1]))
                    and ((z_range is None) or (z>=z_range[0] and z<=z_range[-1]))
                )
                # Load this image
                if in_range_bool:
                    img = io.imread(path_to_images + file)
                    tmp_img_dict[(x, y, z)] = img
            
            elif not has_z:
                # Check if in specificied ranges
                in_range_bool = (
                    ((x_range is None) or (x>=x_range[0] and x<=x_range[-1]))
                    and ((y_range is None) or (y>=y_range[0] and y<=y_range[-1]))
                )
                # Load this image
                if in_range_bool:
                    img = io.imread(path_to_images + file)
                    tmp_img_dict[(x, y)] = img

    # Format into numpy array
    if has_z:
        xyz_pos = np.array([xyz for xyz in tmp_img_dict])
        xy_pos = np.unique(xyz_pos[:, :2], axis=0)
        full_img = np.zeros(
            (len(xy_pos), np.max(xyz_pos[:, -1]))   # xy pos, z stacks
            + img.shape     # width, height, channels
        )
        for i, xy in enumerate(xy_pos):
            x, y = xy
            these_z = np.sort(xyz_pos[np.all(xyz_pos[:, :2] == xy, axis=1), -1])
            full_img[i] = np.stack([tmp_img_dict[(x, y, z)] for z in these_z])
    else:
        xy_pos = np.array([xy for xy in tmp_img_dict])
        full_img = np.zeros(
            (len(xy_pos), 1)   # xy pos, z stacks
            + img.shape     # width, height, channels
        )
        for i, xy in enumerate(xy_pos):
            x, y = xy
            full_img[i] = tmp_img_dict[(x, y)]

    return full_img, xy_pos


# ============================================================================
# Image analysis

def make_bbox(bbox_extents, padding=0):
    """Get the coordinates of the corners of a
    bounding box from the extents

    Parameters
    ----------
    bbox_extents : list (4xN)
        List of the extents of the bounding boxes for each of the N regions.
        Should be ordered: [min_row, min_column, max_row, max_column]

    Returns
    -------
    bbox_rect : np.ndarray
        The corners of the bounding box. Can be input directly into a
        napari Shapes layer.
    """
    minr = bbox_extents[0] - padding
    minc = bbox_extents[1] - padding
    maxr = bbox_extents[2] + padding
    maxc = bbox_extents[3] + padding

    bbox_rect = np.array(
        [[minr, minc], [maxr, minc], [maxr, maxc], [minr, maxc]]
    )
    bbox_rect = np.moveaxis(bbox_rect, 2, 0)

    return bbox_rect


def normalize_clims(image, contrast_limits):
    if image.dtype == 'bool':
        im_norm = image.astype('int16')
    else:
        # normalize data by clims
        im_norm = (
            (image - contrast_limits[0])
            / (contrast_limits[1] - contrast_limits[0])
        )
    im_norm[im_norm > 1] = 1
    return np.array(im_norm)


def normalize_clims_from_viewer(im_layer):
    return normalize_clims(im_layer.data, im_layer.contrast_limits)


def blended_img(image_list, cmap_maps_list, clims_list):
    """
    Get RGB image of given layers
    Adapted from https://forum.image.sc/t/saving-image-from-napari/50379/7

    INPUTS:
    image_list : list of array-like objects, each array-like is one image
        slice for one channel.
    cmap_maps_list : mappings for the given colormap into RGB. Can be found
        in napari viewers with viewer.layers['layer_name'].colormap.map
    clims_list : list length 2, contrast limits for image
    """
    if not len(image_list) == len(cmap_maps_list) == len(clims_list):
        raise ValueError('Lengths of image_list, cmap_maps_list, and '
                         + 'clims_list must match')
    blended = np.zeros(image_list[0].shape + (4,))
    for i in range(len(image_list)):
        normalized_data = normalize_clims(image_list[i], clims_list[i])
        colormapped_data = cmap_maps_list[i](normalized_data.flatten())
        colormapped_data = colormapped_data.reshape(normalized_data.shape
                                                    + (4,))
        blended = blended + colormapped_data
    
    blended[..., 3] = 1 # set alpha channel to 1
    blended[blended > 1] = 1
    return blended


def blended_img_from_viewer(viewer, layers=None, alpha=None):
    """
    Get RGB image of given layers
    Taken from https://forum.image.sc/t/saving-image-from-napari/50379/7
    """
    if layers is None:
        layers_list = viewer.layers
    else: layers_list = layers

    blended = np.zeros(viewer.layers[layers_list[0]].data.shape + (4,))
    for layer_name in layers_list:
        layer = viewer.layers[layer_name]
        normalized_data = normalize_clims(layer)
        colormapped_data = layer.colormap.map(normalized_data.flatten())
        colormapped_data = colormapped_data.reshape(normalized_data.shape
                                                    + (4,))

        if alpha is not None:
            if layer_name in alpha: channel_alpha = alpha[layer_name]
        else: channel_alpha = 1
        blended = blended + (channel_alpha * colormapped_data)
    
    blended[..., 3] = 1 # set alpha channel to 1

    return blended


# ============================================================================
# Saving & displaying images

def plot_img_set_title(ax, label, color='white'):
    ax.text(.01, .99, label, ha='left', va='top', transform=ax.transAxes,
            color=color)
    return ax
    

def plot_img_add_scalebar(ax, dx, units='um', length=None, labeled=False):
    if labeled == False:
        scalebar = ScaleBar(dx, units, fixed_value=length, color='white',
                            location='lower right', frameon=False,
                            label_loc='none', scale_loc='none')
                            # label_formatter = lambda x, y:'')
    else:
        scalebar = ScaleBar(dx, units, fixed_value=length, color='white',
                            location='lower right', frameon=False)
    ax.add_artist(scalebar)
    return ax


def plot_img(ax, image, label, um_per_px, blended=False, color='white',
             scalebar_len=None, scalebar_label=False):
    if not blended:
        ax.imshow(image, cmap='gray', vmin=0, vmax=1)
        if color == 'gray': color = 'white'
        plot_img_set_title(ax, label, color=color)
    elif blended:
        ax.imshow(image)
        plot_img_set_title(ax, label, color='gray')
    
    ax.axis('off')
    plot_img_add_scalebar(ax, um_per_px, length=scalebar_len,
                          labeled=scalebar_label)
    return ax
