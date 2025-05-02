
# VIREO
# See documentation here: https://vireosnp.readthedocs.io/en/latest/manual.html

# Make output folders
out_folder="demultiplex"
mkdir logs -p

# Load modules & environment
module purge
module load gcc

source /home/ng136/miniconda3/etc/profile.d/conda.sh
conda activate base


for lib in Cr_blood_mannitol
do
    # Folder with CellSNP output / SNP matrices
    CELLSNP_OUTPUT="snp_matrices_${lib}"

    # Output folder for demultiplexing
    VIREO_OUTPUT="${out_folder}_${lib}"

    VIREO_EXECUTABLE="/home/ng136/miniconda3/bin/vireo"

    # Input the number of individuals
    N_DONORS=5

    echo Demultiplexing ${lib} \(${N_DONORS} donors\)

    sbatch -p short --job-name snp_demux_SCC_${i} --mem 128000 --cores 4 -t 12:00:00 --mail-type=TIME_LIMIT_90,FAIL,END -o logs/demultiplex_${lib}.o -e logs/demultiplex_${lib}.e --wrap """$VIREO_EXECUTABLE -c $CELLSNP_OUTPUT -N $N_DONORS -o $VIREO_OUTPUT"""


done
