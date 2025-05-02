# ============================================================================
# GOALS OF THIS SCRIPT
# Evaluate the genotype demultiplexing output from Vireo. If needed adjust the
# labels to better reflect annotations we are confident in.
# ============================================================================

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import os, sys
import vireoSNP
# from vireoSNP.utils import io_utils as vio
# from scipy import sparse
from scipy.io import mmread

# Change this path to point to folder containing tal_helper_functions.py
path_to_dropbox = os.environ['PATH_TO_DROPBOX']
sys.path.append(path_to_dropbox + 'klein_lab/resources/helper_functions')
import scrna_helper_functions as hf

# Scanpy settings
sc.settings.verbosity = 2

# Plotting settings
plt.style.use('tal_light_spine')

# ============================================================================
# Function written and shared by Nico

def genotype_distance_matrix(library1, library2, plot=None, verbose=False):
    from vireoSNP.plot import heat_matrix
    from vireoSNP.utils.vcf_utils import match_VCF_samples
    res = match_VCF_samples(f'from_o2/demultiplex_{library1}/GT_donors.vireo.vcf.gz',
                            f'from_o2/demultiplex_{library2}/GT_donors.vireo.vcf.gz',
                            GT_tag1 = 'GT', GT_tag2='GT', )
    if plot == 'matched':
        fig = plt.figure(dpi=150, figsize=(4,4))
        heat_matrix(res['matched_GPb_diff'],
                                  res['matched_donors1'] ,
                                  res['matched_donors2'] , 
                                  interpolation='none')
        plt.ylabel(library1)
        plt.xlabel(library2)
        plt.title(f"{library1} & {library2}\nOverlap: {res['matched_n_var']} SNPs")
        plt.tight_layout()
        # plt.show()
    elif plot == 'full':
        fig = plt.figure(dpi=150, figsize=(4,4))
        heat_matrix(res['full_GPb_diff'],
                                  res['full_donors1'] ,
                                  res['full_donors2'] , 
                                  interpolation='none')
        plt.ylabel(library1)
        plt.xlabel(library2)
        plt.title(f"{library1} & {library2}")
        plt.tight_layout()
        # plt.show()
    title_string = f"{library1} & {library2}"
    res.update({'title': title_string})
    return res


# ============================================================================

version = '4.1'
out_path = f'donor_matching_between_libraries_output/'
if not os.path.exists(out_path): os.mkdir(out_path)
out_path += f'{version}/'
if not os.path.exists(out_path): os.mkdir(out_path)

# genotype_distance_matrix('Cr_blood_1', 'Cr_blood_2', plot='matched')
# plt.savefig(out_path + 'donor_matched.pdf')
# plt.savefig(out_path + 'donor_matched.png', dpi=200)
# plt.close()

# genotype_distance_matrix('Cr_blood_1', 'Cr_blood_2', plot='full')
# plt.savefig(out_path + 'donor_full.pdf')
# plt.savefig(out_path + 'donor_full.png', dpi=200)
# plt.close()

genotype_distance_matrix('C_rob_230518_1', 'C_rob_230518_2', plot='matched')
plt.savefig(out_path + 'donor_matched.pdf')
plt.savefig(out_path + 'donor_matched.png', dpi=200)
plt.close()

genotype_distance_matrix('C_rob_230518_1', 'C_rob_230518_2', plot='full')
plt.savefig(out_path + 'donor_full.pdf')
plt.savefig(out_path + 'donor_full.png', dpi=200)
plt.close()

# --------------

library1 = 'C_rob_230518_1'
library2 = 'C_rob_230518_2'
from vireoSNP.utils.vcf_utils import match_VCF_samples
res = match_VCF_samples(f'from_o2/demultiplex_{library1}/GT_donors.vireo.vcf.gz',
                        f'from_o2/demultiplex_{library2}/GT_donors.vireo.vcf.gz',
                        GT_tag1 = 'GT', GT_tag2='GT', )
df = pd.DataFrame(res['matched_GPb_diff'], index=res['matched_donors1'],
                  columns=res['matched_donors2'])
with plt.style.context('tal_paper_spine'):
    f = plt.figure(figsize=(2.2, 2))
    sns.heatmap(df.iloc[::-1, :],
                xticklabels=True, yticklabels=True, cmap='mako', annot=True)
    plt.xlabel('Library ' + library2[-8:])
    plt.ylabel('Library ' + library1[-8:])
    plt.yticks(rotation=0)
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.savefig(out_path + 'for_paper.pdf')
    plt.close()
