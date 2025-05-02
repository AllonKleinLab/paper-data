# cellsnp-lite
# See documentation here: https://cellsnp-lite.readthedocs.io/en/latest/manual.html, https://github.com/single-cell-genetics/cellsnp-lite

# Make output folder
out_folder="snp_matrices"
mkdir $out_folder -p

mkdir logs -p

module purge
module load gcc

source /home/ng136/miniconda3/etc/profile.d/conda.sh
conda activate base


for lib in Cr_blood_mannitol Cr_blood_asw
do
    # Inputs
    BAM_PATH="/n/groups/klein/tal/220428_ciona_blood_10x_mannitol/10x_raw_data/count_${lib}/outs/possorted_genome_bam.bam"
    #VCF_PATH="None"
    WHITELIST="bc_whitelists/${lib}_bc_whitelist.tsv"
    OUTFOLDER="${out_folder}_${lib}"

    echo ${lib}
    N_CELLS=$(cat ${WHITELIST} | wc -l)
    echo ${N_CELLS}


#    sbatch -p medium --job-name generate_snp_${lib} --mem 128000 --cores 8 -t 48:00:00 --mail-type=TIME_LIMIT_90,FAIL,END -o logs/generate_snp_${lib}.o -e logs/generate_snp_${lib}.e --wrap """cellsnp-lite -s $BAM_PATH -b $WHITELIST -O $OUTFOLDER -R $VCF_PATH -p 20 --minMAF 0.1 --minCOUNT 20 --gzip --genotype"""

    sbatch -p medium --job-name generate_snp_${lib} --mem 128000 --cores 8 -t 05-00:00:00 --mail-type=TIME_LIMIT_90,FAIL,END -o logs/generate_snp_${lib}.o -e logs/generate_snp_${lib}.e --wrap \
        """cellsnp-lite \
            -s $BAM_PATH \
            -b $WHITELIST \
            -O $OUTFOLDER \
            -p 22 \
            --minMAF 0.1 \
            --minCOUNT 100 \
            --chrom Chr1,Chr2,Chr3,Chr4,Chr5,Chr6,Chr7,Chr8,Chr9,Chr10,Chr11,Chr12,Chr13,Chr14,UAContig1,UAContig2,UAContig3,UAContig4,UAContig5,UAContig6,UAContig7,UAContig8,UAContig9,UAContig10,UAContig11,UAContig12,UAContig13,UAContig14,UAContig15,UAContig16,UAContig17,UAContig18,UAContig19,UAContig20,UAContig21,UAContig22,UAContig23,UAContig24,UAContig25,UAContig26,UAContig27,UAContig28,UAContig29,UAContig30,UAContig31,UAContig32,UAContig33,UAContig34,UAContig35,UAContig36,UAContig37,UAContig38,UAContig39,UAContig40,UAContig41,UAContig42,UAContig43,UAContig44,UAContig45,UAContig46,UAContig47,UAContig48,UAContig49,UAContig50,UAContig51,UAContig52,UAContig53 \
            --gzip \
            --genotype"""


done
