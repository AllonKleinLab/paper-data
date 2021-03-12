run=200714_NB501807_0599_AHVW7HBGXF

sbatch -p short --job-name bcl --mem 8000 -t 6:00:00 -o bcl2fq.o -e bcl2fq.o --wrap """ module load bcl2fastq/2.18.0.12; cd $run; bcl2fastq --use-bases-mask y*,y*,y*,y* --mask-short-adapter-reads 0 --minimum-trimmed-read-length 0 """