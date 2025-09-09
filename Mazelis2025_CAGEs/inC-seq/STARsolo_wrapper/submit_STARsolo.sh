#!/bin/sh                                                                                                                                                                                                       

logdir='log_STAR/'
mkdir -p ${logdir}
module load gcc
module load python/3.13.1
module load samtools

n_libs=2
N2=$n_libs
ref = /path/to/STAR/reference/

jid1=$( sbatch -p short --array 1-${N2} -c 6 --job-name STAR  --mem 64000 -t 0-2:00 -o "$logdir"STAR%a.out -e "$logdir"STAR%a_2.out --wrap """ python3 STARsolo_wrapper.py my_analysis.yaml 'map' \$((\$SLURM_ARRAY_TASK_ID-1)) $ref """ )
jid1=${jid1##* }
echo $jid1
