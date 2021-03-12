

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
| Anndata objects combining all libraries from each experiment |  |
| Cell annotation files (i.e. adata.obs) | sandwich_dome_metadata.xlsx ([here](Preprocessing_to_annotation/sandwich_dome_metadata.xlsx)*) |

### Methods

#### Contributors to this repo
Hailey M. Cambra (HMC)
Allon M. Klein (AMK)

#### From reads to counts
The [indrop.py](https://github.com/indrops) pipeline was used for obtain cells x genes matrices. Yaml files detailing the parameters used can be found [here for perturbation data](Perturbation_indrops_scripts) and [here for sandwich and dome culture data](Sandwich_dome_indrops_scripts).


#### Filtering, doublet removal, visualization and annotation
Notebooks starting with "part*" focus on the analysis steps up to filtered annotated data.
Data filtering involved repeating cell visualization several times. For example, drawing and clustering
a kNN graph of T cells only revealed two distinct double populations within T cells, the removal of which
required repeating the visualization of all CD45+ cells.

| Methods | Figure panel(s) | Comment | Relevant notebooks | Contributions |
 ---  | --- | --- | --- | ---
| Quality check raw data, preprocessing, visualization, and annotation | 1C, S1C, 6A, S4A | Filter data for background noise (empty droplets), filter cells and genes, exclude cells with high mitochondrial counts, ribosomal counts, and optionally ncRNA; Counts-normalizing, Logarithmizing, Z-scoring, Dimensionality Reduction (PCA); Clustering, Visualization, and Annotation | [Part1_Preprocess and Annotate Data](Preprocessing_to_annotation) | HMC



### References   
[1] To be updated after publication  

