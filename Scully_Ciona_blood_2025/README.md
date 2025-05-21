This folder contains auxillary code for the preprint: [insert citation]

## Data Availability
scRNA-seq data of _C. robusta_ blood can be viewed interactively and downloaded from [https://kleintools.hms.harvard.edu/paper_websites/scully_ciona_robusta_blood/](https://kleintools.hms.harvard.edu/paper_websites/scully_ciona_robusta_blood/). Downloads are also available from the Gene Expression Omnibus (GEO) repository: Series [GSE296253](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE296253) (full processed data and May 2023 libraries), and Sample GSM8869531 from Series [GSE292926](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE292926) (April 2022 library).

Scripts in this directory expect the _C. robusta_ blood dataset, as well as human and zebrafish datasets, to be downloaded and added to the `data/` subfolder. See that subfolder's `README.md` file for more details.

## Code in This Directory

### Subfolder `helper_functions/`
Contains several helper function files used in scripts throughout the folder.

### Subfolder `scRNA-seq_preprocessing/`
- `genotype_demultiplexing/` contains scripts for carrying out SNP demultiplexing, including generation of plots shown in Supplementary Fig. 1.
- `highly_variable_genes.py` is a python script containing a function to identify highly variable genes following methods in Klein et al. 2015.

### Subfolder `HCR_FISH/`
- `napari_image_analysis_pipeline/` contains the custom image analysis pipeline written with napari, used to process image files and to match cells between live and post-HCR FISH rounds of imaging.

### Subfolder `scRNA-seq_analysis/`
- `go_enrichment_analysis/` contains scripts for performing GO enrichment analysis in _C. robusta_.
- `differentiation_hierarchy_analysis/` contains scripts for performing differentiation hierarchy anlaysis on the data.

### Subfolder `cross_species_comparisons/`
- `gs_plots/` contains scripts for making generalized Sankey (GS) plots and calculating Jaccard similarity scores.
- `coexpression_conservation/` contains scripts for quantifying gene coexpression conservation.
