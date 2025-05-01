#SBATCH -p short                           # Partition to run in
#SBATCH --mem=100G                         # Memory total in MiB (for all cores)
#SBATCH -o logs/genotyping_%a.o            # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e logs/genotyping_%a.e            # File to which STDERR will be written, including job ID (%j)
#SBATCH --mail-type=END                    # Email notification type
#SBATCH --mail-user=tds395@g.harvard.edu   # Email to which notifications will be sent
#SBATCH --array=0-1                        # Split job into this many workers

# List libraries
# ** Make sure to adjust the --array argument of SBATCH to match the number of libraries
lib=("Cr_blood_1" "C_rob_230518_2")

# cellsnp-lite
# See documentation here: https://cellsnp-lite.readthedocs.io/en/latest/manual.html, https://github.com/single-cell-genetics/cellsnp-lite

# Make output folder
out_folder="genotyping"

mkdir logs -p

module purge
module load gcc

source /home/ng136/miniconda3/etc/profile.d/conda.sh
conda activate base


# Inputs
BAM_PATH="/n/groups/klein/tal/230518_ciona_robusta_10X/10x_data/HT_KY21_with_Ensembl_mito/count_${lib[${SLURM_ARRAY_TASK_ID}]}/outs/possorted_genome_bam.bam"
VCF_PATH="/n/groups/klein/tal/230518_ciona_robusta_10X/genotype_demultiplexing/snp_calling_output/cellSNP.base.vcf.gz"
WHITELIST="bc_whitelists/${lib[${SLURM_ARRAY_TASK_ID}]}_bc_whitelist.tsv"
OUTFOLDER="${out_folder}_${lib[${SLURM_ARRAY_TASK_ID}]}"

mkdir $OUTFOLDER -p

echo ${lib[${SLURM_ARRAY_TASK_ID}]}
N_CELLS=$(cat ${WHITELIST} | wc -l)
echo ${N_CELLS}

cellsnp-lite \
    -s $BAM_PATH \
    -R $VCF_PATH \
    -b $WHITELIST \
    -O $OUTFOLDER \
    -p 22 \
    --minMAF 0.1 \
    --minCOUNT 100 \
    --chrom Chr1,Chr2,Chr3,Chr4,Chr5,Chr6,Chr7,Chr8,Chr9,Chr10,Chr11,Chr12,Chr13,Chr14,UAContig1,UAContig2,UAContig3,UAContig4,UAContig5,UAContig6,UAContig7,UAContig8,UAContig9,UAContig10,UAContig11,UAContig12,UAContig13,UAContig14,UAContig15,UAContig16,UAContig17,UAContig18,UAContig19,UAContig20,UAContig21,UAContig22,UAContig23,UAContig24,UAContig25,UAContig26,UAContig27,UAContig28,UAContig29,UAContig30,UAContig31,UAContig32,UAContig33,UAContig34,UAContig35,UAContig36,UAContig37,UAContig38,UAContig39,UAContig40,UAContig41,UAContig42,UAContig43,UAContig44,UAContig45,UAContig46,UAContig47,UAContig48,UAContig49,UAContig50,UAContig51,UAContig52,UAContig53 \
    --gzip \
    --genotype