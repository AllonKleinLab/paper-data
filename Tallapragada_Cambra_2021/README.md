

# ScRNAseq-related analyses for Tallapragada/Cambra et al., (in revision) [1]

### Table of contents
- [Data](#data)
- [Methods](#methods)
  * [Contributors to this repo](#contributors-to-this-repo)
  * [From reads to counts](#from-reads-to-counts)
  * [Filtering, doublet removal, visualization and annotation](#filtering--doublet-removal--visualization-and-annotation)
  * [Analyses using annotated data](#analyses-using-annotated-data)
- [References](#references)

<small><i><a href='http://ecotrust-canada.github.io/markdown-toc/'>Table of contents generated with markdown-toc</a></i></small>

### Data
Large data files used by the code in this repository are stored on [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE164638).
| Description | Filename(s) | 
 ---  | --- 
| Raw count data before filtering | <geo_sample_name_><library_name>.abundant_barcodes.tsv.gz, <geo_sample_name_><library_name>.counts.tsv.gz, <libraray_name>.genes.tsv ([GSE164638](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE161771)|
| Anndata objects combining all libraries from each experiment | adata_annotated_sand_dome.h5ad, adata_merged_annotated_perturbations.h5ad |
| Cell annotation files (i.e. adata.obs) | sandwich_dome_metadata.xlsx and perturbations_metadata.xlsx ([here](Preprocessing_to_annotation)*) |
| Quants metrics files (ouput indrops.py) | _sc_data/quants ([here](_sc_data/quants)*) |

### Methods

#### Contributors to this repo
Hailey M. Cambra (HMC)
Allon M. Klein (AMK)

#### From reads to counts
The [indrop.py](https://github.com/indrops) pipeline was used for obtain cells x genes matrices. Yaml files detailing the parameters used can be found [here for perturbation data](Perturbation_indrops_scripts) and [here for sandwich and dome culture data](Sandwich_dome_indrops_scripts).


#### Filtering, doublet removal, visualization and annotation
Notebooks are organized into directories containing required files to complete analyses. 

Packages not included in the anaconda installation of python, that must be separately installed, to perform the entire analysis for the mouse intestinal organoid data (GEO:GSE164368):

Scanpy (Single-Cell Analysis in Python): https://scanpy.readthedocs.io/en/stable/installation.html
Scrublet (Single-Cell Remover of Doublets): https://github.com/allonkleinlab/scrublet

All input files used to generate single cell RNA sequencing data analysis are provided
in this directory, together with the GEO files. 

| Methods | Figure panel(s) | Comment | Relevant notebooks | Contributions |
 ---  | --- | --- | --- | ---
| Quality check raw data, preprocessing, visualization, and annotation | 1C, S1C, 6A, S4A | Filter data for background noise (empty droplets), filter cells and genes, exclude cells with high mitochondrial counts, ribosomal counts, and optionally ncRNA; Counts-normalizing, Logarithmizing, Z-scoring, Dimensionality Reduction (PCA); Clustering, Visualization, and Annotation | [Part1_Preprocess and Annotate Data](Preprocessing_to_annotation) | HMC
| Calculate and compare abundances of annotated states | 1D, 6B, S4B | Calculate abundances for each cell state; Compare abundances between treatment conditions | [Part2_Abundance Analyses](Abundance_analyses) | HMC
| Further Analysis: Marker gene analysis, Cell Cycle (G2) Score analysis, Gene abundance Analysis  | S1D, 6C, 6E, 6F, S4C | Compare marker genes across states, compare gene expression between treatment conditions, Calculate and plot G2 score | [Part3_Marker Gene and Cell Cylce Score Analysese](Marker_gene_and_cc_score_analyses) | HMC



### References   
[1] To be updated after publication  

