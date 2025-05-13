import os
import sys
import time
# import nd2
import napari
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from nd2reader import ND2Reader
from magicgui import magicgui
from tqdm import tqdm

# Change this path to point to folder containing helper functions scripts
path_to_repo_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(os.path.join(path_to_repo_dir, 'helper_functions'))
import imaging_helper_functions as hf

plt.style.use('tal_light_spine')

# ============================================================================
# SAMPLE INFO

def parse_selected_sample():
    # Parse arguments from command line, if any
    args = ['sample', 'image_folder', 'z_range', 'fov_range']
    arg_dict = hf.get_args([], args)

    # Parse tags
    # --sample, the sample name
    if 'sample' not in arg_dict:
        print('')
        print('PLEASE SPECIFY:')
        print('  --sample : the sample name for analysis')
        print('  --image_folder : path to raw widefield & confocal images')
        print('  --z_range : two integers separated by commas, the range of '
              + 'z slices to use')
        print('  --fov_range : two integers separated by commas, the range of '
              + 'XY positions to use')
        print('')
        raise Exception('Must specify sample with --sample tag')
    else: this_sample = arg_dict['sample']

    # --image_folder, the path to raw images
    if 'image_folder' in arg_dict:
        image_folder = arg_dict['image_folder']
    else:
        image_folder = os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            'raw_images'
        )

    # --z_range, the range of z-slices
    if 'z_range' in arg_dict:
        if ',' in arg_dict['z_range']:
            z_range = [int(x) for x in arg_dict['z_range'].split(',')]
        else:
            raise ValueError('--z_range tag must be followed by two '
                            + 'integers separated by commas.')
    else: z_range = None

    # --fov_range, the range of XY positions
    if 'fov_range' in arg_dict:
        if ',' in arg_dict['fov_range']:
            fov_range = [int(x) for x in arg_dict['fov_range'].split(',')]
        else:
            raise ValueError('--fov_range tag must be followed by two '
                            + 'integers separated by commas.')
    else: fov_range = None

    return [this_sample, image_folder, z_range, fov_range]


# ------------------------------------
# Read additional sample info from SAMPLE_INFO.csv

def import_from_sample_info_file(image_folder, this_sample):
    tmp_df = pd.read_csv(os.path.join(image_folder, 'SAMPLE_INFO.csv'),
                         header=2, index_col='Sample')
    [p1, p2, p3] = tmp_df.loc[this_sample, ['488', '561', '640']].values
    # cluster = tmp_df.loc[this_sample, 'Cluster(s)']
    
    return {p1:'488', p2:'561', p3:'640'}


# ------------------------------------
# Get folder names for wf-live, wf-fixed, and confocal images

def find_image_files(image_folder, this_sample):
    # Widefield live images
    live_path = os.path.join(image_folder, 'widefield')
    file_list = [x for x in os.listdir(live_path) if f'{this_sample}_live' in x]
    if len(file_list) == 1: live_file = os.path.join(live_path, file_list[0])
    elif len(file_list) < 1: raise Exception('No wf-live imaging file found')
    else: raise Exception('More than one wf-live imaging file found')

    # Widefield fixed images
    fixed_path = os.path.join(image_folder, 'widefield')
    file_list = [x for x in os.listdir(fixed_path) if f'{this_sample}_fixed' in x]
    if len(file_list) == 1: fixed_file = os.path.join(fixed_path, file_list[0])
    elif len(file_list) < 1: raise Exception('No wf-fixed imaging file found')
    else: raise Exception('More than one wf-fixed imaging file found')

    # Confocal fluorescence images
    cn_path = os.path.join(image_folder, 'confocal')
    file_list = [x for x in os.listdir(cn_path) if f'{this_sample}' in x]
    if len(file_list) == 1: cn_file = os.path.join(cn_path, file_list[0])
    elif len(file_list) < 1: raise Exception('No confocal imaging file found')
    else: raise Exception('More than one confocal imaging file found')

    return [live_file, fixed_file, cn_file]


# ============================================================================
# Get gene names and expression

def get_gene_names(these_probes):
    path_to_probe_ids = os.path.dirname(os.path.abspath(__file__))
    probes_df = pd.read_excel(os.path.join(path_to_probe_ids,
                                           'HCR_PROBES_MASTER_LIST.xlsx'))
    probes2genes = {probes_df.iloc[i, 0]: probes_df.iloc[i, 1]
                for i in probes_df.index}
    
    gene_list = []
    for p in these_probes:
        # Get gene names
        gene_list.append('.'.join(probes2genes[p].split('.')[:3]))

        # Double check fluorescence channel assignments
        b = probes2genes[p].split('_')[-1]
        b_dict = {'B1':'561', 'B2':'488', 'B3':'640'}
        if b_dict[b] != these_probes[p]:
            Exception(f'{p} ({these_probes[p]}) does not match '
                      + f'{b} probes ({b_dict[b]}).')

    return gene_list


# ============================================================================
# Select subset of FOVs
# (XY positions follow a snake path, starting in the upper left corner)

def select_subset_FOV(live, fixed, cn, fov_range):
    # If fov_range is None, look at all fields of view
    if fov_range is None:
        fov_range = (0, live.sizes['v'])

    # Get number of images on each side of the grids
    if live.sizes['v'] != fixed.sizes['v']:
        raise Exception('Wf live and fixed images have different fields of view')
    cn_fields_side = np.sqrt(cn.sizes['v'])
    wf_fields_side = np.sqrt(live.sizes['v'])
    if not (cn_fields_side.is_integer() and wf_fields_side.is_integer()):
        raise Exception('Image does not have a square number of fields of view')
    cn_fields_side = int(cn_fields_side)
    wf_fields_side = int(wf_fields_side)

    # Get matrix showing order of xy positions
    cn_xy = np.arange(cn_fields_side**2).reshape(cn_fields_side, cn_fields_side)
    for i in range(cn_fields_side):
        if not (i/2).is_integer(): cn_xy[i, :] = cn_xy[i, :][::-1]
    wf_xy = np.arange(wf_fields_side**2).reshape(wf_fields_side, wf_fields_side)
    for i in range(wf_fields_side):
        if not (i/2).is_integer(): wf_xy[i, :] = wf_xy[i, :][::-1]

    # Rotate xy positions for widefield
    # wf_xy_rot = np.rot90(wf_xy, k=-1, axes=(0,1))
    # Not necessary after station 10 upgrade!
    wf_xy_rot = wf_xy

    # Get subset of xy positions
    if cn_fields_side == wf_fields_side:
        wf_xy_cropped = wf_xy_rot
    else:
        diff = wf_fields_side - cn_fields_side
        if diff < 0:
            raise Exception('Too many fields of view in confocal images')
        elif not (diff/2).is_integer():
            raise Exception('Fields of view error: confocal cannot be centered on'
                            + ' widefield fields of view')
        low_idx = int(diff/2)
        high_idx = wf_fields_side - int(diff/2)
        wf_xy_cropped = wf_xy_rot[low_idx:high_idx, low_idx:high_idx]

    cn_xy_idx = cn_xy.flatten()
    wf_xy_idx = wf_xy_cropped.flatten()

    # Use only specified field of view range
    fov_min = max(0, fov_range[0])
    fov_max = min(fov_range[1], len(cn_xy_idx))
    if fov_min is not None and fov_max is not None:
        cn_xy_idx = cn_xy_idx[fov_min:fov_max]
        wf_xy_idx = wf_xy_idx[fov_min:fov_max]

    return wf_xy_idx, cn_xy_idx


# ============================================================================
# Convert from PIM to numpy array

def pim_to_array(img, channel_order, z_range=None, xy_list=None, maxip=True):
    channels = [ch.lower() for ch in img.metadata['channels']]
    full_image = []

    if xy_list is None: xy_list = range(img.sizes['v'])
    if 'z' in img.iter_axes:
        if (z_range is None): z_range = range(img.sizes['z'])
        else: z_range = range(z_range[0], z_range[1])

    # Loop through xy ('v') positions)
    for xy in tqdm(xy_list):
        # Dict ensures channels are saved in expected order
        this_image = {ch_name: [] for ch_name in channel_order}

        # Loop through channels
        for i, ch in enumerate(channels):

            # For widefield: no separate z slices
            if 'z' not in img.iter_axes:
                this_image[ch] = img[img.sizes['c'] * xy + i]   # correct xy

            # For confocal: has z slices, option to take MaxIP
            elif maxip:
                # Take bottom slice for dic
                if 'dic' in ch:
                    this_image[ch] = img[
                        img.sizes['z'] * img.sizes['c'] * xy    # correct xy
                        + i * img.sizes['z']    # correct channel
                        + 1]                    # correct z slice
                # MaxIP for fluorescence channels (all others)
                else:
                    channel_maxip = np.zeros(img[0].shape)
                    for z in z_range:
                        z_slice = img[
                            img.sizes['z'] * img.sizes['c'] * xy  # correct xy
                            + i * img.sizes['z']    # correct channel
                            + z]                    # correct z slice
                        channel_maxip = np.maximum(channel_maxip, z_slice)
                    this_image[ch] = channel_maxip

            # For confocal: has z slices, option to keep z stacks
            else:
                channel_stack = np.zeros((img.sizes['z'], img[0].shape[0],
                                          img[0].shape[1]))
                # Loop through and save all z slices
                for z in z_range:
                    z_slice = img[
                        img.sizes['z'] * img.sizes['c'] * xy  # correct xy
                        + i * img.sizes['z']    # correct channel
                        + z]                    # correct z slice
                    channel_stack[z, :, :] = z_slice
                
                this_image[ch] = channel_stack
        
        # Add image with channels (as multi-dim array) to full_image list
        # print(this_image)
        full_image.append(np.stack([this_image[ch] for ch in this_image]))

    # Save as multidimensional array
    return np.stack(full_image)


# ============================================================================
# Load images from PIM, then to numpy array, then display in napari viewer

def load_images_into_viewer(live_file, fixed_file, cn_file, z_range,
                            fov_range, wf_channel_order, cn_channel_order,
                            channel_names_dict):
    # --------------------------------
    # Load images
    print('Loading images...')
    t = time.time()

    live = ND2Reader(live_file)
    voxel_live = live.metadata['pixel_microns']
    fixed = ND2Reader(fixed_file)
    voxel_fixed = fixed.metadata['pixel_microns']
    cn = ND2Reader(cn_file)
    voxel_cn = cn.metadata['pixel_microns']

    # Set up multidimensional slicing
    # See http://soft-matter.github.io/pims/v0.6.1/multidimensional.html
    for img in [live, fixed, cn]:
        if img == cn: img.iter_axes = ['v', 'c', 'z']
        else: img.iter_axes = ['v', 'c']
        img.bundle_axes = ['y', 'x']

    # Get subset of images (xy positions) to import
    [wf_xy_idx, cn_xy_idx] = select_subset_FOV(live, fixed, cn, fov_range)

    # Convert to numpy arrays
    print('Live widefield images')
    live_img = pim_to_array(live, wf_channel_order, xy_list=wf_xy_idx)
    print('Fixed widefield images')
    fixed_img = pim_to_array(fixed, wf_channel_order, xy_list=wf_xy_idx)
    print('Confocal images (with MaxIP)')
    cn_img = pim_to_array(cn, cn_channel_order, z_range, cn_xy_idx)

    # Rotate images from widefield, since camera is rotated
    # EDIT: Not necessary after station 10 upgrade!
    # live_img_rot = np.rot90(live_img, k=-1, axes=(2, 3))
    live_img_rot = live_img
    # fixed_img_rot = np.rot90(fixed_img, k=-1, axes=(2, 3))
    fixed_img_rot = fixed_img

    # --------------------------------
    # Display in viewer
    viewer = napari.Viewer()

    # Widefield images
    # wf_channel_names = ['DIC', 'DIC_Red', 'DIC_Green', 'DIC_Blue', 'DAPI']
    wf_channel_names = [channel_names_dict[c] for c in wf_channel_order]
    viewer.add_image(live_img_rot, channel_axis=1,
                    scale=[voxel_live, voxel_live],
                    colormap=['gray', 'red', 'green', 'blue', 'blue'],
                    name=['live_pre-FISH_'+x for x in wf_channel_names])
    viewer.add_image(fixed_img_rot, channel_axis=1,
                    scale=[voxel_fixed, voxel_fixed],
                    colormap=['gray', 'red', 'green', 'blue', 'blue'],
                    name=['fixed_pre-FISH_'+x for x in wf_channel_names])

    # Confocal images
    cn_channel_names = [channel_names_dict[c] for c in cn_channel_order]
    viewer.add_image(cn_img, channel_axis=1, scale=[voxel_cn, voxel_cn],
                     colormap='gray',
                     name=(['post-FISH_' + cn_channel_names[0]]
                           + cn_channel_names[1:]))
    
    # Add scale bar
    viewer.scale_bar.visible = True
    viewer.scale_bar.unit = "um"

    print('Total time to load: ', time.time() - t)
    return viewer


# ============================================================================
# Adjust contrast

@magicgui(call_button='(Re)load contrast limits')
def load_clims(viewer, out_path):
    try:
        contr_limits = pd.read_csv(os.path.join(out_path, 'contrast_limits.txt'),
                                   index_col=0)
        # Change channel colors
        for ch in viewer.layers:
            if 'points_' not in ch.name:
                viewer.layers[ch.name].contrast_limits = [
                    contr_limits.loc[ch.name, 'contrast_limits1'],
                    contr_limits.loc[ch.name, 'contrast_limits2']
                ]
        # viewer.status = '(Re)loaded contrast limits'
        print('(Re)loaded contrast limits')
    except:
        # viewer.status = 'No contrast limits found - please set new ones'
        print('No contrast limits found - please set new ones')


@magicgui(call_button='Save current contrast limits')
def save_current_clims(viewer, out_path):
    with open(os.path.join(out_path, 'contrast_limits.txt'), 'w') as f:
        f.write('channel,contrast_limits1,contrast_limits2\n')
        for layer in viewer.layers:
            if 'points_' not in layer.name:
                contr_limits_str = [str(x) for x in layer.contrast_limits]
                f.write(f'{layer.name},{",".join(contr_limits_str)}\n')

    # viewer.status = 'Contrast limits saved'
    print('Contrast limits saved')


@magicgui(call_button="Fluorescent channels gray")
def channels_to_gray(viewer, channel_list):
    for c in channel_list:
        viewer.layers[c].colormap = 'gray'

@magicgui(call_button="Fluorescent channels color")
def channels_to_color(viewer, channel2color_dict):
    for c in channel2color_dict:
        viewer.layers[c].colormap = channel2color_dict[c]

# ============================================================================
# Adjust layer positioning

@magicgui(call_button='Shift live images',
          x={"min": -float('inf'), "max": float('inf')},
          y={"min": -float('inf'), "max": float('inf')})
def shift_live_images(viewer, x:float=0, y:float=0):
    live_image_layers = [i for i in viewer.layers if ('live' in i.name) and
                                                     ('point' not in i.name)]
    for layer in live_image_layers:
        layer.translate = [0, y, x]

# ============================================================================
# Set up viewer for point placement

def make_point_layers(viewer, point_layers):
    for layer in point_layers:
        viewer.add_points(name=layer, scale=viewer.layers[0].scale, ndim=3,
                          size=30, face_color='#ffb300', edge_color='#ffb300',
                          symbol=point_layers[layer])


def set_marker_genes_visible(viewer, gene_list):
    # Set layer visibility for point placement
    for layer in viewer.layers:
        lname = layer.name
        if (lname in [gene_list[1], gene_list[2]]) or ('points_' in lname):
            layer.visible = True
        else:
            layer.visible = False


def set_layer_visibility(viewer, layer_list):
    # Set layer visibility for point placement
    for layer in viewer.layers:
        lname = layer.name
        if (lname in layer_list) or ('points_' in lname):
            layer.visible = True
        else:
            layer.visible = False


def reorder_layers(viewer, gene_list):
    # Reorder layers for point placement
    layer_order_old = [l.name for l in viewer.layers]
    layer_order_new = (
        [l for l in layer_order_old if 'live_pre-FISH_DIC_' in l]
        + ['live_pre-FISH_DAPI']
        + [l for l in layer_order_old if 'fixed_pre-FISH_DIC_' in l]
        + ['fixed_pre-FISH_DAPI']
        + ['live_pre-FISH_DIC']
        + ['fixed_pre-FISH_DIC']
        + ['post-FISH_DIC']
        + ['DAPI']
        + gene_list
        + [l for l in layer_order_old if 'points_' in l]
    )
    viewer.layers[:] = [viewer.layers[name] for name in layer_order_new]

    set_marker_genes_visible(viewer, gene_list)

    # Go to first XY position (0)
    current_step = list(viewer.dims.current_step)
    current_step[0] = 0
    viewer.dims.current_step = tuple(current_step)


# ============================================================================
# Cropped images of each cell

# ------------------------------------
# GUI function to save cropped cell figures

@magicgui(call_button='Clear points layers')
def clear_points_layers(viewer):
    for layer_name in point_layers:
        viewer.layers[layer_name].data = []


@magicgui(call_button='Save cropped cell images')
def save_cropped_cell_images(viewer, gene_list, out_path, cell_id: int = 0,
                             fluor_point=None, fixed_point=None,
                             live_point=None, marker_pos_cell=True,
                             box_size=15, dpi=200):
    # Set up folders for saved images
    if marker_pos_cell:
        out_path2 = hf.add_output_path(os.path.join(out_path, 'aligned_cells_marker_pos'))
        out_path_live = hf.add_output_path(os.path.join(out_path2, 'with_preFISH_live_fixed'))
        out_path_fixed = hf.add_output_path(os.path.join(out_path2, 'with_preFISH_fixed'))
        out_path_fluor = hf.add_output_path(os.path.join(out_path2, 'no_preFISH'))
    else:
        out_path2 = hf.add_output_path(os.path.join(out_path, 'aligned_cells_marker_neg'))
        out_path_live = hf.add_output_path(os.path.join(out_path2, 'with_live_only'))
        out_path_fixed = hf.add_output_path(os.path.join(out_path2, 'with_live_and_fixed'))
        out_path_fluor = hf.add_output_path(os.path.join(out_path2, 'with_FISH'))
    out_path_info = hf.add_output_path(os.path.join(out_path2, 'cell_info'))

    # Re-set contrast limits
    # load_clims()

    # Get points
    if (fluor_point is None) and (fixed_point is None) and (live_point is None):
        fluor_point = (viewer.layers['points_fluor(Control)'].data
                       - (viewer.layers['post-FISH_DIC'].translate
                          /viewer.layers['post-FISH_DIC'].scale))
                    #    + viewer.layers['points_fluor(Control)'].translate)
        fixed_point = (viewer.layers['points_fixed(Alt)'].data
                       - (viewer.layers['fixed_pre-FISH_DIC'].translate
                          /viewer.layers['fixed_pre-FISH_DIC'].scale))
                    #    + viewer.layers['points_fixed(Alt)'].translate)
        live_point = (viewer.layers['points_live(Shift)'].data
                      - (viewer.layers['live_pre-FISH_DIC'].translate
                         /viewer.layers['live_pre-FISH_DIC'].scale))

    # Define box size for rows and columns
    # for box_size in box_size_list:
    box_r_um = box_c_um = box_size
    box_r_px = np.ceil(box_r_um / viewer.layers[0].scale[1])
    box_c_px = np.ceil(box_c_um / viewer.layers[0].scale[2])

    # Define channels (& sets of channels) to plot
    if False:#this_sample == 'sx':
        """
        On a per-experiment basis, use this section if a sample has two marker
        genes for two different clusters. This way overlays can be saved
        separately for different clusters, rather than including an overlay of
        genes which label different cell states.

        For this example, not relevant, so it will always be skipped
        """
        if marker_pos_cell:
            gene_list_for_overlay1 = [gene_list[1]]
            gene_list_for_overlay2 = [gene_list[2]]
        else:
            gene_list_for_overlay1 = [gene_list[0], gene_list[1]]
            gene_list_for_overlay2 = [gene_list[0], gene_list[2]]
        viewer_layer_list = [
            ['live_pre-FISH_DIC'],
            [f'live_pre-FISH_DIC_{x}' for x in ['Red','Green','Blue']],
            ['live_pre-FISH_DAPI'],
            ['fixed_pre-FISH_DIC'],
            [f'fixed_pre-FISH_DIC_{x}' for x in ['Red','Green','Blue']],
            ['fixed_pre-FISH_DAPI'],
            ['post-FISH_DIC'],
            ['DAPI'] + gene_list_for_overlay1,
            ['DAPI'] + gene_list_for_overlay2,
            ['DAPI'],
            [gene_list[0]],
            [gene_list[1]],
            [gene_list[2]],
        ]
    else:
        if marker_pos_cell:
            gene_list_for_overlay = gene_list[1:]
        else:
            gene_list_for_overlay = gene_list
        viewer_layer_list = [
            ['live_pre-FISH_DIC'],
            [f'live_pre-FISH_DIC_{x}' for x in ['Red','Green','Blue']],
            ['live_pre-FISH_DAPI'],
            ['fixed_pre-FISH_DIC'],
            [f'fixed_pre-FISH_DIC_{x}' for x in ['Red','Green','Blue']],
            ['fixed_pre-FISH_DAPI'],
            ['post-FISH_DIC'],
            ['DAPI'] + gene_list_for_overlay,
            ['DAPI'],
            [gene_list[0]],
            [gene_list[1]],
            [gene_list[2]],
        ]

    # Initialize file for cell positioning info
    with open(os.path.join(out_path_info, f'cell{cell_id:02d}_info.csv'), 'w') as f:
        f.write(f'Channels,img,center_row,center_col\n')

    # Save this cell as a figure
    # nrow, ncol = hf.get_nrow_ncol(len(viewer_layer_list), nrow=1)
    nrow = 1; ncol = len(viewer_layer_list)
    f = plt.figure(figsize=(2 * ncol, 2 * nrow))
    count_no_image = 0
    for i, layer_list in enumerate(viewer_layer_list):
        ax = plt.subplot(nrow, ncol, i+1)

        image_list = []; clims_list = []; cmap_maps_list = []
        for image_layer in layer_list:
            if 'live_pre-FISH_' in image_layer:
                this_cell = live_point
            elif 'fixed_pre-FISH_' in image_layer:
                this_cell = fixed_point
            else:
                this_cell = fluor_point

            if len(this_cell) != 0:
                this_cell = this_cell[0]

                # Get image with this cell
                img_idx = int(this_cell[0])
                this_img = viewer.layers[image_layer].data[img_idx, ...]
                
                # Calculate edges of box
                center_r, center_c = int(this_cell[1]), int(this_cell[2])
                max_r = min(center_r + int(box_r_px/2), this_img.shape[0])
                min_r = max(0, center_r - int(box_r_px/2))
                max_c = min(center_c + int(box_c_px/2), this_img.shape[1])
                min_c = max(0, center_c - int(box_c_px/2))
                this_cell_img = this_img[min_r:max_r, min_c:max_c]

                image_list.append(this_cell_img)
                clims_list.append(viewer.layers[image_layer].contrast_limits)
                cmap_maps_list.append(viewer.layers[image_layer].colormap.map)

        # Plot this image
        # Multichannel overlay
        if len(image_list) > 1:
            blended = hf.blended_img(image_list, cmap_maps_list, clims_list)
            hf.plot_img(ax, blended, 'Overlay', blended=True,
                    um_per_px=viewer.layers[0].scale[1], scalebar_len=5)
        # Single channel
        elif len(image_list) == 1:
            image_layer = layer_list[0]
            this_img_norm = hf.normalize_clims(image_list[0], clims_list[0])
            hf.plot_img(ax, this_img_norm, ' '.join(image_layer.split('_')),
                    um_per_px=viewer.layers[0].scale[1], scalebar_len=5,
                    color=viewer.layers[image_layer].colormap.name)
            # Save cell positioning info
            with open(os.path.join(out_path_info, f'cell{cell_id:02d}_info.csv'), 'a') as f:
                f.write(f'{image_layer},{img_idx},{center_r},{center_c}\n')
        # No image to plot
        elif len(image_list) == 0:
            ax.axis('off')
            # Save cell positioning info
            with open(os.path.join(out_path_info, f'cell{cell_id:02d}_info.csv'), 'a') as f:
                f.write(f'{image_layer},CELL MISSING OR UNIDENTIFIABLE,,\n')
            # Track which output folder to save image to
            count_no_image += 1
            
    # Set output folder to save for marker_pos_cell
    if marker_pos_cell:
        if (len(live_point) == 0) and (len(fixed_point) == 0):
            this_output_folder = out_path_fluor
        elif (len(live_point) == 0) and (len(fixed_point) > 0):
            this_output_folder = out_path_fixed
        elif len(live_point) > 0:
            this_output_folder = out_path_live
    # Set output folder to save for marker_neg_cell
    else:
        if (len(fluor_point) > 0):
            this_output_folder = out_path_fluor
        elif (len(fixed_point) > 0):
            this_output_folder = out_path_fixed
        else:
            this_output_folder = out_path_live

    plt.tight_layout()
    plt.savefig(os.path.join(this_output_folder, f'cell{cell_id:02d}_box{box_size}um.pdf'),
                dpi=dpi)
    plt.close()

    # Reset/initalize for next cell
    save_cropped_cell_images.cell_id.value += 1
    clear_points_layers()
    if marker_pos_cell: set_marker_genes_visible(viewer, gene_list)
    else: set_layer_visibility(viewer, ['live_pre-FISH_DIC'])


# ============================================================================
# Regenerate cell files

def get_coords_from_cell_info(cell_info_filename, path_to_cell_info):    
    # Extract info from csv file for relevant channels
    with open(os.path.join(path_to_cell_info, cell_info_filename), 'r') as f:
        lines = f.readlines()
    channel_list = np.array([ch.split(',')[0] for ch in lines])
    ch_idx_live = np.where(channel_list == 'live_pre-FISH_DIC')[0][0]
    ch_idx_fixed = np.where(channel_list == 'fixed_pre-FISH_DIC')[0][0]
    ch_idx_cn = np.where(channel_list == 'post-FISH_DIC')[0][0]
    info_live = lines[ch_idx_live].strip().split(',')
    info_fixed = lines[ch_idx_fixed].strip().split(',')
    info_cn = lines[ch_idx_cn].strip().split(',')

    # Format as points arrays compatible with viewer points layers
    empty_points = np.array([]).reshape((0, 3))
    xyz_live, xyz_fixed, xyz_cn = empty_points, empty_points, empty_points
    if not 'CELL MISSING OR UNIDENTIFIABLE' in info_live:
        xyz_live = np.array([float(info_live[1]), float(info_live[2]), 
                                float(info_live[3])]).reshape((1, 3))
    if not 'CELL MISSING OR UNIDENTIFIABLE' in info_fixed:
        xyz_fixed = np.array([float(info_fixed[1]), float(info_fixed[2]),
                                float(info_fixed[3])]).reshape((1, 3))
    if not 'CELL MISSING OR UNIDENTIFIABLE' in info_cn:
        xyz_cn = np.array([float(info_cn[1]), float(info_cn[2]),
                            float(info_cn[3])]).reshape((1, 3))
    
    return (xyz_live, xyz_fixed, xyz_cn)


@magicgui(call_button='Regenerate images from\nexisting cell coordinates')
def regenerate_cell_images(viewer, gene_list, out_path, marker_pos_cell=True,
                           this_cell_id:int=999):
    # Get path to csv files
    # cell_info_path = out_path + 'aligned_cells/cell_info/'
    # Set up folders for saved images
    if marker_pos_cell:
        out_path2 = hf.add_output_path(os.path.join(out_path, 'aligned_cells_marker_pos'))
    else:
        out_path2 = hf.add_output_path(os.path.join(out_path, 'aligned_cells_marker_neg'))
    out_path_info = hf.add_output_path(os.path.join(out_path2, 'cell_info'))
    csv_files = np.sort([f for f in os.listdir(out_path_info)
                         if f[:4] == 'cell'])
    
    to_print = 'marker+' if marker_pos_cell else 'marker-'
    print(f'Regenerating cropped cell images ({to_print})')
    
    if this_cell_id == 999:
        for c in tqdm(csv_files):
            cell_id = c.split('_')[0].split('cell')[-1]
            cell_id = int(cell_id)
            xyz_live, xyz_fixed, xyz_cn = get_coords_from_cell_info(c,
                                                                    out_path_info)

            save_cropped_cell_images(viewer, gene_list, cell_id=cell_id,
                                    fluor_point=xyz_cn, fixed_point=xyz_fixed,
                                    live_point=xyz_live, out_path=out_path,
                                    marker_pos_cell=marker_pos_cell)
    
    else:
        cell_id = this_cell_id
        csv_file = f'cell{cell_id:02d}_info.csv'
        xyz_live, xyz_fixed, xyz_cn = get_coords_from_cell_info(csv_file,
                                                                out_path_info)

        save_cropped_cell_images(viewer, gene_list, cell_id=cell_id,
                                fluor_point=xyz_cn, fixed_point=xyz_fixed,
                                live_point=xyz_live, out_path=out_path,
                                marker_pos_cell=marker_pos_cell)
    
    print('Done')


# ============================================================================

# def main():
#     ...

# ============================================================================

if __name__ == "__main__":
    # main()

    # Create output folder
    try: out_path = hf.add_output_path(
        os.path.abspath(__file__).split('.')[0] + '_output')
    except: out_path = hf.add_output_path('cropped_cell_images_output')

    # Prompt to get info about sample
    [this_sample, image_folder, z_range, fov_range] = parse_selected_sample()
    [live_file, fixed_file, cn_file] = find_image_files(image_folder,
                                                        this_sample)
    these_probes = import_from_sample_info_file(image_folder, this_sample)
    gene_list = get_gene_names(these_probes)

    out_path = hf.add_output_path(os.path.join(out_path, this_sample))

    # Set desired channel order
    wf_channel_order = ['dic n2', 'dic n2 red', 'dic n2 green', 'dic n2 blue',
                        'widefield blue']
    cn_channel_order = ['dic n2', 'confocal 405', 'confocal 488',
                        'confocal 561', 'confocal 640']
    channel_names_dict = {'dic n2': 'DIC', 'dic n2 red': 'DIC_Red',
                          'dic n2 green': 'DIC_Green',
                          'dic n2 blue': 'DIC_Blue',
                          'widefield blue': 'DAPI',
                        #   'dic': 'post-FISH_DIC',
                          'confocal 405': 'DAPI',
                          'confocal 488': gene_list[0],
                          'confocal 561': gene_list[1],
                          'confocal 640': gene_list[2]}

    # Load images
    viewer = load_images_into_viewer(live_file, fixed_file, cn_file, z_range,
                                     fov_range, wf_channel_order,
                                     cn_channel_order, channel_names_dict)

    cn_channel2color_dict = {
        'post-FISH_DIC': 'gray',
        'DAPI': 'blue',
        gene_list[0]: 'yellow', # ubiquitous control
        gene_list[1]: 'green',  # marker 2
        gene_list[2]: 'magenta',# marker 1
    }

    # Set up for point placements
    point_layers = {
        'points_live(Shift)': {'shape': 'disc', 'key':'Shift'},
        'points_fixed(Alt)': {'shape': 'diamond', 'key':'Alt'},
        'points_fluor(Control)': {'shape': 'cross', 'key':'Control'},
    }
    make_point_layers(viewer,
                      {p:point_layers[p]['shape'] for p in point_layers})
    reorder_layers(viewer, gene_list)

    # Custom mouse function for placing points
    # See https://napari.org/stable/gallery/custom_mouse_functions.html
    @viewer.mouse_drag_callbacks.append
    def add_cell_point(viewer, event):
        new_point = np.divide(event.position, viewer.layers[0].scale)
        # print('Click at', event.position)
        # print('New point at', new_point)
        # print(event.modifiers)
        for l in point_layers:
            if point_layers[l]['key'] in event.modifiers:
                # viewer.layers[l].add(new_point)
                # print('mouse click', new_point)
                # print('translate', viewer.layers[l].translate)
                viewer.layers[l].data = new_point

    # Add widgets (setting variables as needed)
    save_cropped_cell_images.viewer.value = viewer
    save_cropped_cell_images.gene_list.value = gene_list
    save_cropped_cell_images.out_path.value = out_path
    viewer.window.add_dock_widget(save_cropped_cell_images)

    shift_live_images.viewer.value = viewer
    viewer.window.add_dock_widget(shift_live_images)

    clear_points_layers.viewer.value = viewer
    viewer.window.add_dock_widget(clear_points_layers)

    load_clims.viewer.value = viewer
    load_clims.out_path.value = out_path
    viewer.window.add_dock_widget(load_clims)

    save_current_clims.viewer.value = viewer
    save_current_clims.out_path.value = out_path
    viewer.window.add_dock_widget(save_current_clims)

    channels_to_color.viewer.value = viewer
    channels_to_color.channel2color_dict.value = cn_channel2color_dict
    viewer.window.add_dock_widget(channels_to_color)

    channels_to_gray.viewer.value = viewer
    channels_to_gray.channel_list.value = cn_channel2color_dict
    viewer.window.add_dock_widget(channels_to_gray)

    regenerate_cell_images.viewer.value = viewer
    regenerate_cell_images.gene_list.value = gene_list
    regenerate_cell_images.out_path.value = out_path
    viewer.window.add_dock_widget(regenerate_cell_images)

    napari.run()
