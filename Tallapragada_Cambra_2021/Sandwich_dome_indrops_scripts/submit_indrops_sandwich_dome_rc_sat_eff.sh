
logdir="logs_rc/" # path for writing log files
yaml="20200722_sandwich_dome_rc_sat_eff.yaml"
nlibs=5 # number of libraries to process
nlanes=4


######################################################
# In general, no need to change anything below here. #
######################################################

N1=$nlanes
N3=$nlibs # number of workers for "sort" step
N4=200 # number of workers for "quantify" step

# bug recently introduced that messes up multiple
# workers in the aggregate step
# N5=$nlibs # number of workers for "aggregate" step

module load gcc/6.2.0
module load bowtie/1.2.2
module load java/jdk-1.8u112
module load rsem/1.3.0

mkdir -p ${logdir}

pyenv=/n/groups/klein/sam/pyndrops
indrops=/n/groups/klein/sam/pipelines/indrops_common/indrops/indrops.py
myPython=/n/groups/klein/sam/pyndrops/bin/python
source ${pyenv}/bin/activate ${pyenv}

## PART 4 (quantify)
jid4=$( sbatch -p short --array 1-${N4} --job-name Q --mem 3000 -t 06:00:00 -o "$logdir"quant_worker_%a.out -e "$logdir"quant_worker_%a.out --wrap """ ${myPython} $indrops $yaml quantify --min-reads 250 --min-counts 0 --total-workers ${N4} --worker-index \$((\$SLURM_ARRAY_TASK_ID-1)) """)
jid4=${jid4##* }
echo $jid4

## PART 5 (aggregate)
# max number of workers = number of libraries
jid5=$( sbatch --dependency=afterany:$jid4 -p short --job-name A --mem 8000 -t 02:00:00 -o "$logdir"agg.out -e "$logdir"agg.out --wrap """ ${myPython} $indrops $yaml aggregate """ )
jid5=${jid5##* }
echo $jid5
