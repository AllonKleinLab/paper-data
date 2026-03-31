#!/bin/bash
#SBATCH -c 3                               # Request three cores
#SBATCH -t 2:00:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=100G                          # Memory total in MiB (for all cores)
#SBATCH --mail-type=FAIL                   # ALL email notification type
#SBATCH --mail-user=q_wu@g.harvard.edu     # Email to which notifications will be sent

# Define the output folder as a variable
OUTPUT_DIR="/n/groups/klein/qiu/exp_0070_analysis/220718_Thorough_Analysis/240315_Paper_Writing/Code/IL17_Synergy_Project_Code/Job_Statuses/"

# Ensure the directory exists (create it if necessary)
mkdir -p "$OUTPUT_DIR"

# Capture the hostname
HOSTNAME=$(hostname)

# Set the filenames for STDOUT and STDERR using job name (%x), job ID (%j), and node name (%N)
OUT_FILE="${OUTPUT_DIR}${SLURM_JOB_NAME}_${SLURM_JOB_ID}.out"
ERR_FILE="${OUTPUT_DIR}${SLURM_JOB_NAME}_${SLURM_JOB_ID}.err"

# Redirect STDOUT and STDERR to the files
exec > "$OUT_FILE" 2> "$ERR_FILE"

### JOB CODE ####

module load gcc/6.2.0
source activate pydeseq2_env
cd /n/groups/klein/qiu/exp_0070_analysis/220718_Thorough_Analysis/240315_Paper_Writing/Code/IL17_Synergy_Project_Code/GLM_Implementation/

python 241002_GLM_Update.py

rm *.out