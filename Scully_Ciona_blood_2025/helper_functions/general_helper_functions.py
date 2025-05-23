import os
import sys
import getopt
import numpy as np
import pickle
from zipfile import ZipFile
from types import SimpleNamespace
from scipy.cluster.hierarchy import dendrogram, linkage

# ============================================================================
# Paths to data and useful folders

path = SimpleNamespace()

path.repo_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# !! FIRST MAKE SURE TO DOWNLOAD C. ROBUSTA DATA INTO THE RIGHT SUBFOLDER
# See data/README.md
# Path to C. robusta blood dataset, downloadable from GEO:
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE296253
# GSE296253_combined_processed_data.h5ad
path.Crob_adata_file = os.path.join(
    path.repo_dir,
    'data',
    'GSE296253_combined_processed_data.h5ad'
)

# External datasets
path.external_datasets = os.path.join(
    path.repo_dir,
    'data',
)

# ============================================================================
# Gene homology and annotations

path_to_dict = os.path.join(path.repo_dir, 'helper_functions', 'gene_dictionaries')

# Use gene homology determined by OrthoFinder
with open(os.path.join(path_to_dict, 'crob2hsap.pickle'), 'rb') as fp:
    ciona2human = pickle.load(fp)
    human2ciona = pickle.load(fp)
# with open(path_to_dict + 'hsap2drer.pickle', 'rb') as fp:
#     ciona2human = pickle.load(fp)
#     human2ciona = pickle.load(fp)

# TF gene names
with open(os.path.join(path_to_dict, 'tf_dict.pickle'),'rb') as fp:
    tf_dict = pickle.load(fp)

# ============================================================================
# Species class for easy name variations

class Species:
    def __init__(self, genus, species):
        genus = genus.title()
        species = species.lower()
        self.genus = genus
        self.species = species
        self.name_full = f'{genus} {species}'
        self.name = f'{genus[0]}. {species}'
        self.Xxyz = f'{genus[0]}{species[:3]}'
        self.Xx = f'{genus[0]}{species[0]}'

    def __str__(self):
        return self.name_Xx


# ============================================================================
# Creating output paths for scripts

def add_output_path(dir_name=''):
    if not os.path.exists(dir_name): os.mkdir(dir_name)
    return dir_name


# ============================================================================
# Miscellaneous

def get_args(short_list, long_list=None):
    """
    Parse and return command line arguments in a dictionary
    """
    # Format input variables
    short_str = ':'.join(short_list) + ':'
    arg_dict = {f'-{i}': i for i in short_list}
    if long_list is not None:
        long_str = [i + ' =' for i in long_list]
        arg_dict.update({f'--{i} ': i for i in long_list})

    try:
        opts, _ = getopt.getopt(sys.argv[1:], short_str, long_str)
        
    except:
        short_error = ['-'+i for i in short_list].join(', ')
        long_error = ['--'+i for i in long_list].join(', ')
        print(f"Input error. Expect {short_error} or {long_error}")

    return {arg_dict[i[0]]:i[1] for i in opts}


def hierarchical_clustering_leaf_order(matrix, method='single',
    metric='cosine', labels=None, optimal_ordering=False):
    linkage_data = linkage(matrix, method=method, metric=metric,
                           optimal_ordering=optimal_ordering)
    dendr = dendrogram(linkage_data, labels=labels, count_sort='descending')
    return dendr['ivl'][::-1]
