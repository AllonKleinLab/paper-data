import os
import sys
import pandas as pd
import numpy as np
import gzip
from collections import defaultdict
import ast

#Input     
yaml_file = sys.argv[1]
run_mode= sys.argv[2]
file_to_run= sys.argv[3]
ref= sys.argv[4]


#loading input files from ymal file
def get_files_dir_from_yaml(yaml_file):
    temp1 = []
    temp2 = []
    index = []
    with open(yaml_file, 'r') as file:
        for line in file:
            line.split(':')[0]
            temp1.append(line.split(':')[0])
            temp2.append(line.split(':')[1].strip('\n'))
            if len(line.split(':'))>2:
                index.append(line.split(':')[2].strip('\n')[2:-2])
    
    df = pd.DataFrame([temp1, temp2]).T

    names = [i.split(',')[0][2:-1] for i in (df[df[0] == '       - {library_name'])[1]]
    num_samples = len(names)
    df.loc[(df[df[0] == '       - {library_name']).index.values, 1] = index
    df.loc[(df[df[0] == '       - {library_name']).index.values, 0] = names
    df.index = df[0]
    del df[0]
    del temp1
    del temp2
    
    #--readFilesIn   
    files = []
#    for run in range(0, df.loc['    dir '].shape[0]):
#        files.append(['/'.join(df.loc['    dir '].iloc[run].values.astype(str)[0][2:-1].split('/')[:-5]) + '/' + df.iloc[1].values.astype('str')[0][2:-1] + i + '/filtered_parts/' for iG, i in enumerate(names)])
    files = [os.getcwd() + '/' + df.iloc[1].values.astype('str')[0][2:-1] + i + '/filtered_parts/' for iG, i in enumerate(names)]
    files = np.array(files)
    files = np.unique(files.flatten())
    
    all_files = []
    library = []
    
    for run_loc in files:
        transcript_file = []
        bcd_file = []
        for file in os.listdir(run_loc):
            if file.endswith('.fastq.gz') & (not file.endswith('_barcodes.fastq.gz')):
                all_files.append(run_loc + file)
                transcript_file.append(run_loc + file)
            elif file.endswith('_barcodes.fastq.gz'):
                bcd_file.append(run_loc + file)
        transcript_file.sort()
        bcd_file.sort()
        library.append([transcript_file, bcd_file])
    return [all_files, library, names]

#running STARsolo
def run_STAR(yaml_file, file_to_run, ref'):
    file_to_run = int(file_to_run)
    file_locs, names = get_files_dir_from_yaml(yaml_file)[1:]
    os.system('path/to/STAR-2.7.11b/bin/Linux_x86_64_static/STAR --genomeDir {4} \
--runThreadN 6 \
--readFilesCommand gunzip -c \
--readFilesIn {1} {2} \
--outFileNamePrefix {3}/results_STAR/{0}_ \
--soloCBstart 1 \
--soloCBlen 24 \
--soloUMIstart 25 \
--soloUMIlen 10 \
--soloBarcodeReadLength 0 \
--soloType CB_UMI_Simple \
--soloCBwhitelist None \
--soloFeatures Gene GeneFull \
--genomeSAsparseD 3 \
--outSAMtype BAM SortedByCoordinate \
--outSAMattributes CB CR UR UB UY GN GX \
--outFilterScoreMinOverLread 0.3 \
--outFilterMatchNminOverLread 0.3 \
--soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
--soloUMIfiltering MultiGeneUMI_CR \
--soloUMIdedup 1MM_CR \
--soloMultiMappers EM \
--soloCellFilter EmptyDrops_CR \
--quantMode TranscriptomeSAM GeneCounts '.format(names[file_to_run], ','.join(file_locs[file_to_run][0]), ','.join(file_locs[file_to_run][1]), '/'.join(file_locs[file_to_run][0][0].split('/')[:-3]), ref))

#Mode to run
if run_mode == 'map':
    run_STAR(yaml_file, file_to_run, ref)
