run=/n/groups/klein/hailey/organoids_single_cell/FC_06195/n/files/Genetics/BPF-NGS/nextseq/200827_NB501715_0612_AHFHKFBGXG

sbatch -p short --job-name bcl --mem 8000 -t 6:00:00 -o bcl2fq.o -e bcl2fq.o --wrap """ module load bcl2fastq/2.18.0.12; cd $run; bcl2fastq --use-bases-mask y*,y*,y*,y* --mask-short-adapter-reads 0 --minimum-trimmed-read-length 0 """