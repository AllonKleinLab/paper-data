import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
import h5py
import os

# Scanpy settings
# verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.verbosity = 2

# Set random seed
np.random.seed(seed = 0)

# ============================================================================
# PATHS TO DATA

# Filtered & preprocessed adata
# v1: UMI/bc thresholds used in preprint
#version = 'v1'
# v2: UMI/bc thresholds used in Figure S1
#version = 'v2'
# v3: using 10X Cell Ranger's automatic barcode filtering
version = 'v3'

adata_path = (f'filtering_and_preprocessing_output/{version}/')

# Unfiltered counts matrix, molecule info
data_path = 'raw_unfiltered_data/'

# Create output path
out_path = 'data_quality_metrics/'
if not os.path.exists(out_path): os.mkdir(out_path)
out_path += version + '/'
if not os.path.exists(out_path): os.mkdir(out_path)

# Samples/libraries
libs = [
    'GSM8869530_Cr_blood_asw',
    'GSM8869531_Cr_blood_mannitol'
]

# ============================================================================
# LOOP THROUGH LIBRARIES

for l in libs:

    print(l)

    # --------------------------------
    # 1. Get reads and UMIs per barcode
    # Read HDF5 file: https://docs.h5py.org/en/latest/quick.html
    # Molinfo contents: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/molecule_info
    with h5py.File(data_path + l + '_molecule_info.h5', "r") as molinfo:

        # Barcode corresponding to each barcode index
        bc_dict = dict([*enumerate(np.array(molinfo['barcodes']))])

        # Molecules associated to each barcode
        _, umis_per_barcode = np.unique(np.array(molinfo['barcode_idx']),
                                        return_counts=True)

        # Genes associated to each barcode
        genes_per_barcode = pd.DataFrame(
                {
                    'bc_idx': np.array(molinfo['barcode_idx']),
                    'feature_idx': np.array(molinfo['feature_idx'])
                },
            ).groupby(['bc_idx', 'feature_idx']).size().groupby(['bc_idx']).size()

        # Reads per barcode
        per_bc_df = pd.DataFrame(
                {
                    'bc_index': np.array(molinfo['barcode_idx']),
                    'reads': np.array(molinfo['count'])
                },
            ).groupby('bc_index').sum()
        per_bc_df['UMIs'] = umis_per_barcode
        per_bc_df['genes'] = genes_per_barcode
        per_bc_df.index = per_bc_df.index.map(lambda x : bc_dict[x].decode("utf-8"))

    # --------------------------------
    # 2. Add info from adata

    # Import adata
    adata = sc.read_h5ad(adata_path + 'adata.h5ad')
    adata = adata[adata.obs['sample'] == l]

    # Remove "-1" from the end of adata barcode indices
    bc_adata = adata.obs.index.map(lambda x: x.split('-')[0])

    # Label in per_bc_df which cells passed the filter
    per_bc_df['passed_filter'] = per_bc_df.index.isin(bc_adata)

    # Add reads per barcode to adata
    adata.obs['reads'] = per_bc_df.loc[bc_adata, :]['reads'].tolist()

    # --------------------------------
    # 3. Save csv summary of data quality statistics

    with open(out_path + f'{l}_quality_stats.csv', 'w') as f:

        # Data about read mapping
        f.write('MAPPING')

        f.write('\nValid reads mapped confidently to transcriptome,')
        f.write(repr(np.sum(per_bc_df['reads'])))

        # Data about barcodes which are cells
        f.write('\n\nCELLS')

        f.write('\nEstimated number of cells,')
        f.write(repr(np.sum(per_bc_df['passed_filter'])))

        f.write('\nReads confidentally mapped to transcriptome in droplets with cells,')
        num = np.sum(per_bc_df.loc[per_bc_df['passed_filter'], 'reads'])
        denom = np.sum(per_bc_df['reads'])
        f.write(f'{num} ({np.round(100*num/denom, 1)}%)')

        f.write('\nMean reads per cell,')
        f.write(repr(per_bc_df.loc[per_bc_df['passed_filter'], 'reads'].mean()))
        f.write('\nMedian reads per cell,')
        f.write(repr(per_bc_df.loc[per_bc_df['passed_filter'], 'reads'].median()))

        f.write('\nMean UMIs per cell,')
        f.write(repr(per_bc_df.loc[per_bc_df['passed_filter'], 'UMIs'].mean()))
        f.write('\nMedian UMIs per cell,')
        f.write(repr(per_bc_df.loc[per_bc_df['passed_filter'], 'UMIs'].median()))

        f.write('\nMean genes per cell,')
        f.write(repr(per_bc_df.loc[per_bc_df['passed_filter'], 'genes'].mean()))
        f.write('\nMedian genes per cell,')
        f.write(repr(per_bc_df.loc[per_bc_df['passed_filter'], 'genes'].median()))

        f.write('\nTotal genes detected,')
        f.write(repr(np.sum(adata.var[f'total_counts-{l}'] != 0)))

        f.write('\n')
