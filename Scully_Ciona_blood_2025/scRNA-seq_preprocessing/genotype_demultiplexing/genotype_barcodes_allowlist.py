import scanpy as sc
import subprocess
import os

out_path = 'genotype_barcodes_allowlist_output/'
if not os.path.exists(out_path): os.mkdir(out_path)

# ============================================================================
# GENEROUSLY FILTER DATA

# Run filtering.py with very generous thresholds to get all barcodes we might
# want to have genotyping information on
bashCommand = ('python3 filtering.py --version version_for_genotype '
               + '--max_mito_pct 20,20 --min_num_UMI 500,800 '
               + '--min_num_genes 0,0 --run_scrublet False')
bashCommand_list = bashCommand.split(' ')
process = subprocess.run(bashCommand, shell=True)

# ============================================================================
# IMPORT FILTERED ADATA

print('Getting allowlisted barcodes')

# Get paths to folders, create folder for saving figures
filt_path = 'filtering_output/version_for_genotype/'
sample = 'Cr_blood_mannitol'    # library name

# Import expression data
adata = sc.read_h5ad(filt_path + f'filtered_{sample}.h5ad')

# ============================================================================
# FOR EACH LIBRARY, EXPORT LIST OF CELL BARCODES
#
# To check the right barcode formatting I ran the following in an interactive
# O2 session:
#
# module load gcc/9.2.0
# module load samtools/1.14
# samtools view count_Cr_blood_mannitol/outs/possorted_genome_bam.bam | head -n 1
# 
# This returned the first line of that BAM file, which included
# "CB:Z:ACCCTCAAGACAACAT-1"
# We want the barcode to match this formatting.

libs = [str(x) for x in adata.obs.library.unique()]

for lib in libs:
    bc_list = adata[adata.obs.library == lib].obs.index.unique()
    bc_fmt = ['-'.join(x.split('-')[:2]) for x in bc_list]

    with open(out_path + f'{lib}_bc_allowlist.tsv', 'w') as f:
        for bc in bc_fmt: f.write(bc + '\n')
