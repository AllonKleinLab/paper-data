This folder contains auxillary code for the preprint: [insert citation]

## Data Availability
scRNA-seq data of _C. robusta_ blood can be viewed interactively and downloaded from [https://kleintools.hms.harvard.edu/paper_websites/scully_ciona_robusta_blood/](https://kleintools.hms.harvard.edu/paper_websites/scully_ciona_robusta_blood/). Downloads are also available from the Gene Expression Omnibus (GEO) repository: Series GSExxxxx (full processed data and May 2023 libraries), and Sample GSM8869531 from Series [GSE292926](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE292926) (April 2023 library).

## Code in This Directory

Mentioned in methods sections:
- GS Plots

### Subfolder `scRNA-seq_preprocessing/`
- `genotype_demultiplexing/` contains scripts for carrying out SNP demultiplexing, including generation of plots shown in Supplementary Fig. 1.
- `highly_variable_genes.py` is a python script containing a function to identify highly variable genes following methods in Klein et al. 2015.

### Subfolder `HCR_FISH/`
- `napari_image_analysis_pipeline/` contains the custom image analysis pipeline written with napari, used to process image files and to match cells between live and post-HCR FISH rounds of imaging.

### Subfolder `scRNA-seq_analysis/`
- `go_enrichment_analysis/` contains scripts for performing GO enrichment analysis in _C. robusta_.
- `differentiation_hierarchy_analysis/` contains scripts for performing differentiation hierarchy anlaysis on the data.

### Subfolder `cross_species_comparisons/`
- `gs_plots/` contains scripts for making generalized Sankey (GS) plots.