

# ScRNAseq-related analyses for Tallapragada/Cambra et al., (in revision) [1]

### Table of contents
- [Data](#data)
- [Methods](#methods)
  * [Contributors to this repo](#contributors-to-this-repo)
  * [From reads to counts](#from-reads-to-counts)
  * [Filtering to annotation](#filtering-to-annotation)
  * [Analyses using annotated data](#analyses-using-annotated-data)
- [References](#references)
- [SPRING](#SPRING-Links)

<small><i><a href='http://ecotrust-canada.github.io/markdown-toc/'>Table of contents generated with markdown-toc</a></i></small>

### Data

Large data files used by the code in this repository are stored on [GEO:GSE164638](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE164638)

| Description | Filename(s) | 
 ---  | --- 
| Raw count data before filtering | <geo_sample_name_><library_name>.abundant_barcodes.tsv.gz, <geo_sample_name_><library_name>.counts.tsv.gz, <geo_sample_name_><libraray_name>.genes.txt.gz [GSE164638](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE164638)|
| Anndata objects combining all libraries from each experiment | GSE164638_adata_mg_annotated_merged_coarse_grain_no_ss_12092020.h5ad.gz (Sandwich and Dome data as in publication); GSE164638_adata_perturbations_merged_annotated_cg_122220.h5ad.gz (perturbation data) [GSE164638](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE164638) |
| Cell annotation files (i.e. adata.obs) | sandwich_dome_metadata.xlsx and perturbations_metadata.xlsx ([here](Preprocessing_to_annotation)*) |
| Quants metrics files (ouput indrops.py) | _sc_data/quants ([here](_sc_data/quants)*) |

### Methods

#### Contributors to this repo

Hailey M. Cambra (HMC) <br>
Allon M. Klein (AMK) <br>

#### From reads to counts

The [indrop.py](https://github.com/indrops) pipeline was used for obtain cells x genes matrices. Yaml files detailing the parameters used and scripts used to run the indrops.py pipeline can be found [here for perturbation data](Perturbation_indrops_scripts) and [here for sandwich and dome culture data](Sandwich_dome_indrops_scripts).


#### Filtering to annotation

Notebooks are organized into directories containing required files to complete analyses. 

Packages not included in the anaconda installation of python, that must be separately installed, to perform analysis:

Scanpy (Single-Cell Analysis in Python): https://scanpy.readthedocs.io/en/stable/installation.html <br>
Scrublet (Single-Cell Remover of Doublets): https://github.com/allonkleinlab/scrublet <br>

All input files used to generate single cell RNA sequencing data analysis are provided
in this directory, together with the GEO files. 

| Methods | Figure panel(s) | Comment | Relevant notebooks | Contributions |
 ---  | --- | --- | --- | ---
| Quality check raw data, preprocessing, visualization, and annotation | 1C, S1C, 6A, S4A | Filter data for background noise (empty droplets), filter cells and genes, exclude cells with high mitochondrial counts, ribosomal counts, and optionally ncRNA; Normalization and Dimensionality Reduction (PCA); Clustering, visualization (SPRING), and annotation | [Part1_Preprocess and Annotate Data](Preprocessing_to_annotation) | HMC



#### Analyses using annotated data

After cells are annotated by cell state, analyze the data for trends in gene expression and state characteristics. 

| Methods | Figure panel(s) | Comment | Relevant notebooks | Contributions |
 ---  | --- | --- | --- | ---
| Calculate and compare abundances of annotated states | 1D, 6B, S4B | Calculate abundances for each cell state; Compare abundances between treatment conditions | [Part2_Abundance Analyses](Abundance_analyses) | HMC
| Marker gene analysis, Cell cycle (G2) score analysis, Gene abundance analysis  | S1D, 6C, 6E, 6F, S4C | Compare marker genes across states, compare gene expression between treatment conditions, calculate and plot G2 score | [Part3_Marker Gene and Cell Cylce Score Analyses](Marker_gene_and_cc_score_analyses) | HMC, AMK
| Gene set enrichment analysis  | 6I | Compare genes enriched in stretch response state to similar cell states reported in literature | [Part4_Gene Set Enrichment Analysis](Gene_set_enrichment_analysis) | HMC, AMK

### SPRING Links

[Inflation perturbation dataset](https://kleintools.hms.harvard.edu/tools/springViewer_1_6_dev.html?datasets/SPRING_private/mouse_intestinal_organoids_21/murine_intestinal_organoid_inflation_perturbation_expt/all_cells)

[Sandwich and Dome culture dataset (no stretch signature annotation)](https://kleintools.hms.harvard.edu/tools/springViewer_1_6_dev.html?datasets/SPRING_private/mouse_intestinal_organoids_21/murine_intestinal_organoid_sandwich_dome_culture_merged_no_stretch_signature/all_cells)

[Sandwich and Dome culture dataset (with stretch signature annotation)](https://kleintools.hms.harvard.edu/tools/springViewer_1_6_dev.html?datasets/SPRING_private/mouse_intestinal_organoids_21/murine_intestinal_organoid_sandwich_dome_culture_merged/all_cells)

### References 

[1] To be updated after publication <br>
[2] Klein, A. M., Mazutis, L., Akartuna, I., Tallapragada, N., Veres, A., Li, V., Peshkin, L., Weitz, D. A., & Kirschner, M. W. (2015). Droplet barcoding for single-cell transcriptomics applied to embryonic stem cells. Cell, 161(5), 1187–1201. [https://doi.org/10.1016/j.cell.2015.04.044](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4441768/) <br>
[3] Weinreb, C., Wolock, S., & Klein, A. M. (2018). SPRING: a kinetic interface for visualizing high dimensional single-cell expression data. Bioinformatics (Oxford, England), 34(7), 1246–1248. [https://doi.org/10.1093/bioinformatics/btx792](https://pubmed.ncbi.nlm.nih.gov/29228172/) <br>
[4] Wolf, F. A., Angerer, P., & Theis, F. J. (2018). SCANPY: large-scale single-cell gene expression data analysis. Genome Biology, 19(1), 15. [https://doi.org/10.1186/s13059-017-1382-0](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-017-1382-0) <br>
[5] Wolock, S. L., Lopez, R., & Klein, A. M. (2019). Scrublet: Computational Identification of Cell Doublets in Single-Cell Transcriptomic Data. Cell systems, 8(4), 281–291.e9. [https://doi.org/10.1016/j.cels.2018.11.005](https://pubmed.ncbi.nlm.nih.gov/30954476/) <br>
