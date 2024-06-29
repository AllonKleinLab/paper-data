This directory contains codes for data analyses in Kukreja et. al. 2024 Cell cycle manuscript. ([Link to preprint](https://doi.org/10.1101/2023.07.29.551123))

There are two folders: scRNA-seq_data_analysis and Image_analysis

### **scRNA-seq_data_analysis**:
This folder contains:

- cell_cycle_env.yml: This is the conda environment used to run python jupyter notebooks for scRNA-seq data analysis.
- pip_packages_list.txt: List of packages installed obtainied using the command pip list.
- helper_functions_py3.py: Functions called while running analysis notebooks.
- Notebooks:
	- 1.preprocessing_data.ipynb: used for filtering cells based on total UMI counts and saving data in .h5ad anndata format which is used in subsequent notebooks.
	- 2.classification_14hpf.ipynb: example notebook for classifying cells into cell states.
	- 3.Fig1_find_top_marker_genes.ipynb: used for finding marker genes for each cell state.
	- 4.Fig1_compare_top_markers.ipynb: used for calculating shared markers between different conditions (Fig. 1h).
	- 5.Fig1_number_of_diff_genes_from6.ipynb: used for calculating number of differentially expressed genes starting from 6 hours to later time points (Fig. 1i).
	- 6.Fig1_number_of_leiden_clusters.ipynb: used for calculating number of Leiden clusters at each time point (Fig. 1j).
	- 7.Fig4_cosine_distance_highly_var_genes.ipynb: used for finding cosine distance between different cell states at different time points (Fig. 4b).
	- 8.23-02-21_global_programs.ipynb: used to make plots related to expression of cell division arrest specific global programs (Fig. 4e-h).
	- 9.percent_TTV.ipynb: used for calculating total transcriptomic variance (TTV) and % TTV as defined in Fig. 4c. This notebook generates the plot Fig. 4d.
	- 10.diff_delay_by_classification.ipynb: used for calculating differentiation delay (Fig. 3).
	- 11.abundance_analysis_hua.ipynb: example notebook used for calculating fold change abundance of cell types. This notebook calculates fold change abundance for 24 hpf HUA embryos (Fig. 5a)
	- 12.beta_binomial_test_abundance_hua.ipynb: example notebook for calculating p-values for fold change abundances.
	- 13.R_ibb_script_hua.txt: the notebook above (12.beta_binomial_test_abundance_hua.ipynb) generates this script which needs to be run on R to get p-values.
	- 14.proliferation_score_tree.ipynb: used for calculating proliferation scores and generating Fig. 5g.
	- 15.FinalTricycleAnalysis.R: R script used for calculating cell cycle phase angle.  
- processing_raw_reads: contains scripts for getting counts data from raw fastq reads. the pipeline used can be found here: https://github.com/AllonKleinLab/klunctions/tree/master/sam/inDrops_pipeline

### **Image_analysis folder**:
This folder contains:

- Notebooks:
	- 20231214_z20c_Voltron_GCaMP_Heartbeat_Plotting.ipynb: used for analysis of cardiac calcium imaging data.
	- 20240116_phalloidin_head_volume_measurement.ipynb: used for analysis of cell size vs tissue size for embryo head (Fig. 5).
	- example_plantseg_configuration: configuration files necessary to generate PlantSeg membrane segmentations for cell size vs. tissue size analysis (Fig. 5).
	- 20240120_emi_pH3_run_cellpose.ipynb: used for analysis of phospho-histone3 immunostaining.
	- 20240409_tissue_nuclear_size_comparison_emi_ctrl.ipynb: used for segmenting nuclei in emi1 and control embryos.
	- 20240423_HUA_nuclear_size_run_cellpose.ipynb: used for segmenting nuclei for high-background HUA embryos.
	- CP_20240424_nuclear_high_bg: human-in-the-loop updated Cellpose model for high-background nuclear stain.
