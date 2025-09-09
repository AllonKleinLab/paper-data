# Python code for demultiplexing and filtering raw fastq files from inC-seq protocols.

# inC-seq library layout
Read1 file includes the cDNA sequence (or genome fragment sequence in amplicon and ATAC-seq experiments).  <br />
Read2 file contains a 6bp barcode that is appended during the 3rd barcoding PCR step.  <br />
Read3 file contains sample p5 index (used for demultiplexing), added during the library prep after capsule breakage.  <br />
Read4 contains Ligation-2 barcode [0:6], Ligation-1 barcode [10:16], RT/PCR/Tn5 [21:27] barcode and UMI [27:37] (in the case of scRNAseq). <br />


Adapted and updated from https://github.com/indrops/indrops

# Principle
This code takes raw fastq files, demultiplexes them based on library index sequences in the yaml file, filters the low quality reads with incorect barcodes and outputs two fastq files: 1 for cDNA sequence and 1 for barcode sequence.
The following two files can then be used in STARsolo pipeline 

STARsolo_wrapper contain the script to read .yaml file, infer the outputs of filter step and run STARsolo

### Filtering step
This iterates over sequencing run parts, optionally filtered by a list of sequencing parts, and a list of libraries of interest.

    python inCseq_filtering.py project.yaml filter [--total-workers 1] [--worker-index 0]
          [-r --runs RUNS ] [-l --libraries LIBRARIES ]

    # --runs comma-separated list of runs :               If specified, step will be restricted to run parts coming
    # --libraries comma-separated list of libraries :      If specified, step will be restricted to run parts that 
    #                                                     contain reads from a library in the list
    # 
    # Resulting workload (a list of run parts), will be split among N --total-workers,
    # where worker with --worker-index i will do steps (i, N+i, 2N+i, ...)

This step reads the raw FastQ files as input and filters them:
  - For every raw read, it determines if the read has the expected structure (depending on library version). 
  - For reads with correct structure, it runs Trimmomatic.
  - For reads surviving Trimmomatic, it finds and trims the polyA tail a maximum length of 4, and checks if the reads are still above MIN_LEN.
  - For surviving reads, it determines which fraction of the read is composed of runs of the same base (only considering runs of 5 or more). 
    It rejects reads whose fraction is greater than `low_complexity_filter_arguments:max_low_complexity_fraction`.

As output, for every input run part, this produces a filtered FastQ file for every library contained in that run for the read1 sequence and fully assembled barcode.

A log is created detailing what happened to every input read. An index is created that lists the number of reads found for every barcode. 

## Project YAML file.

A project is composed of a series of sequencing runs, each including one (or several) inC-seq libraries within it. A sequencing run can further be split into several parts (usually by sequencing) to parallelize the analysis. Give example of a chunk
The project yaml file contains the details of all sequencing runs and libraries within a project. 


## Project YAML file

An example YAML file is provided in `my_analysis.yaml`. It should contain the following information:
Markdown
    project_name : "project_name"
    project_dir : "/path/to/project/dir"  #This dir should be user-owned and writable, all output will go into this dir.
    paths : 
      bowtie_index : "/path/to/index" #This index will be built automatically
      # The paths below can be omitted if the relevant directories are already on $PATH
      python_dir : "/path/to/env/bins/"
      java_dir: "/path/to/java/dir/"
      rsem_dir: "/path/to/rsem/dir/"
      samtools_dir: "/path/to/samtools-1.3.1/bin/" #This needs to be version 1.3.1, 1.3 is not good enough!

    sequencing_runs : 
      # A list of sequencing runs which form the project. 
      # Each run should have:
      - name : "MyRun" # The name of the run will be used as a prefix in filenames, so keep it sane.
        version : "vN" # Can be 'v1', 'v2' or 'v3'

      # For a run with a single 'part', and a single library
        dir : "/path/to/run_files/"
        fastq_path : "{read}.fastq.gz" # Read with be replaced by R1, R2, R3, R4 as appropriate.
        library_name : "my_library"

        # This will expect to find the files:
        #    /path/to/run_files/R1.fastq.gz (and R2...)

       # For a run with several parts, but a single library
        dir : "/path/to/run_files/"
        fastq_path : "{split_affix}_{read}.fastq.gz" # Read with be replaced by R1, R2, R3, R4 as appropriate.
        split_affixes : ["L001", "L002"]
        library_name : "my_library"

        # This will expect to find the files:
        #    /path/to/run_files/L001_R1.fastq.gz (and R2...)
        #    /path/to/run_files/L002_R1.fastq.gz (and R2...)

       # For a run with several parts, several libraries, that have already been demultiplexed
        dir : "/path/to/run_files/"
        fastq_path : "{library_prefix}_{split_affix}_{read}.fastq.gz" # Read with be replaced by R1, R2, R3, R4 as appropriate.
        split_affixes : ["L001", "L002"]
        libraries : 
          - {library_name: "test_lib1", library_prefix: "lib1"}
          - {library_name: "test_lib2", library_prefix: "lib2"}

        # This will expect to find the files:
        #    /path/to/run_files/lib1_L001_R1.fastq.gz (and R2...)
        #    /path/to/run_files/lib1_L002_R1.fastq.gz (and R2...)
        #    /path/to/run_files/lib2_L001_R1.fastq.gz (and R2...)
        #    /path/to/run_files/lib2_L002_R1.fastq.gz (and R2...)

       # For a V3 run with several parts, with several libraries that are not already demultiplexed
        dir : "/path/to/run_files/"
        fastq_path : "{library_prefix}_{split_affix}_{read}.fastq.gz" # Read with be replaced by R1, R2, R3, R4 as appropriate.
        split_affixes : ["L001", "L002", "L003", "L004"]
        libraries :  # The library index is what the expected index read sequence (on a NextSeq, this is the reverse complement of the index sequence)
          - {library_name: "test_lib3", library_index: "ATAGAG"}
          - {library_name: "test_lib4", library_index: "AGAGGA"}

        # This will expect to find the files:
        #    /path/to/run_files/lib1_L001_R1.fastq.gz (and R2, R3, R4...)
        #    /path/to/run_files/lib1_L002_R1.fastq.gz (and R2, R3, R4...)
        #    /path/to/run_files/lib1_L003_R1.fastq.gz (and R2, R3, R4...)
        #    /path/to/run_files/lib1_L004_R1.fastq.gz (and R2, R3, R4...)



