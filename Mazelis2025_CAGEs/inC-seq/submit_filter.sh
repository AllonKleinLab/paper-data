logdir="logs_out/" # path for writing log files
yaml="my_analysis.yaml"
nlibs=2 # number of libraries to process
nlanes=1 # number of lanes

######################################################
# In general, no need to change anything below here. #
######################################################

N1=$nlanes
N3=$nlibs # number of workers for "sort" step

module load gcc
module load perl/5.40.1
module load java/jdk-23.0.1
module load samtools
mkdir -p ${logdir}

pyenv=/path/to/py27/
indrops=inC-seq_filtering.py
myPython=/path/to/py27//bin/python
source ${pyenv}/bin/activate ${pyenv}


## Run Filtering step
jid1=$( sbatch -p short --array 1-${N1} --job-name F --mem 48000 -t 1:00:00 -o "$logdir"filter_worker_%a.out -e "$logdir"filter_worker_%a.out --wrap """ ${myPython} ${indrops} $yaml filter --total-workers ${N1} --worker-index \$((\$SLURM_ARRAY_TASK_ID-1)) """ )
jid1=${jid1##* }
echo $jid1


