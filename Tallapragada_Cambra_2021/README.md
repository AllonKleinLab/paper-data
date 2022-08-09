

# ScRNAseq-related analyses for Tallapragada, Cambra et al., 2021

### Table of contents
- [Data](#data)
- [Methods](#methods)
  * [Contributors to this repo](#contributors-to-this-repo)
  * [Triple-decker sandwich cultures](#triple-decker-sandwich-cultures)
  * [From reads to counts](#from-reads-to-counts)
  * [Filtering to annotation](#filtering-to-annotation)
  * [Analyses using annotated data](#analyses-using-annotated-data)
- [SPRING](#SPRING-Links)
- [References](#references)

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

#### Triple-decker sandwich cultures

[See our publication on the method used in the paper in Current Protocols [6]](https://pubmed.ncbi.nlm.nih.gov/35030297/). Recently, we've been having a lot of trouble with the sandwich cultures lifting off of the coverslip. This could be for a number of reasons, but we've found that the PolyHEMA is the common denominator in the lifing we have observed. One way to get around this is to use [Aquapel](http://www.aquapel.com/) instead. This is a reagent that is routinely used to coat the microchips used in our inDrops protocol, and the use in that protocol inspired us to try this reagent as a replacement for PolyHEMA. It is a rain-repellant used on windshields, and is quite inexpensive. The trade-off to using this versus the PolyHEMA is that it can be more difficult to evenly spread, which can result in monolayer formation when the organoids eventually touch the coverslip. To achieve full coating of the surface area, you will need to repeat the coating process multiple times. We have grown organoids in sandwich cultures with Aquapel coating and not observed any differences in growth. In the future, we will aim to publish this protocol more widely in one form or another, but for now, we're sharing this for those of you running into the same issue. The protocol for obtaining Aquapel from it's container (be careful!), and it's application to glass coverslips is as follows: <br>

Reagents and Materials: <br>
* Double gloves <br>
* A chemical hood with sash <br>
* A glass syringe or other glass container without excess air for Aquapel <br>
* Glass bottom dish or coverslip to be coated <br>
* A 15 mL falcon tube for breaking open the capsule <br> 
* 0.2 µm PTFE syringe filter to fit onto the syringe <br> 
* 1 mL pipette and pipette tips <br>
* 25 G 5/8" needle <br>
* 0.015" I.D. x 0.043" O.D. 95 Durometer LDPE Micro Medical Tubing <br>
* Sterile hood for applying to glass coverslips <br>
* Vacuum desiccator <br>
* Parafilm <br>
* Aquapel <br>

**Use full PPE during this protocol! Glasses, lab coat, and break the glass capsule inside the chemical fume hood in step 4!**

Protocol: <br>
1. Make your own dish-friendly applicator by cutting a piece of the applicator fabric and stuffing into the back of a 200 µL pipette tip. <br>
2. To get to the glass capsule inside, very carefully cut open the fabric applicator on the outside with one long cut down the middle, to release the capsule. Do not leave the Aquapel in plastic containers (including falcon tubes, etc.) for long periods of time, and do not allow prolonged exposure to the air. <br>
3. Prepare the syringe by adding the 0.2 µm PTE syringe filter onto the opening. <br>
4. Break open the glass container in the 15 mL falcon tube by placing the glass capsule inside and forcing the cap on. The capsule will shatter and drain it's contents into the 15 mL falcon tube. <br>
5. Empty the syringe of air and use a 1 mL pipette to obtain the contents into the pipette tip and use to transfer through the filter into the syringe. You can pull as you push the contents from the pipette tip. <br>
6. Remove any excess air and add the needle to the syringe. You can keep the cap or alternatively, use micro tubing which you can then tie and tape the end of to prevent air from coming into contact with the solution. The latter is how we normally store our Aquapel, but this is because it has multiple uses in our lab. If you need to use Auqapel for this purpose only, feel free to use whatever glass container is useful for you; it's just critical that the Aquapel is not stored in plastic and is kept from exposure to the air. <br>
7. Apply one to two drops of the Aquapel at a time to the coverslip, and use the applicator (made in step 1) to spread it and remove excess Aquapel. <br>
8. Place the coated glass-bottom dish or coverslip into a vacuum desiccator and seal with Parafilm. In about 2-4 hours the coat should have dried. <br>
9. Repeat steps 7-8 for a total of at least three times. <br>

Troubleshooting: <br>
* If you see monolayers forming on patches of your coverslip, then there is insufficient Aquapel coating on that region of the coverslip. Be sure to get the entire surface area (it's easy to miss the very edges, or the very middle if you become too concerned about missing the edges). <br> 
* Check your coverslips before you add any Aquapel -- are there imperfections in the glass? Is it actually entirely clean? This will matter for any coating you add. <br>

Thanks to Alex Ratner and Ignas Mazelis for their guidance in opening up the Aquapel capsules safely and storing the material properly!

#### From reads to counts

The [indrop.py](https://github.com/indrops) pipeline was used for obtain cells x genes matrices. Yaml files detailing the parameters used and scripts used to run the indrops.py pipeline can be found [here for perturbation data](Perturbation_indrops_scripts) and [here for sandwich and dome culture data](Sandwich_dome_indrops_scripts).


#### Filtering to annotation

Notebooks are organized into directories containing required files to complete analyses. 

Packages not included in the anaconda installation of python, that must be separately installed, to perform analysis:

Scanpy (Single-Cell Analysis in Python): https://scanpy.readthedocs.io/en/stable/installation.html <br>
Scrublet (Single-Cell Remover of Doublets): https://github.com/allonkleinlab/scrublet <br>

All input files used to generate single cell RNA sequencing data analysis are provided
in this directory, together with the GEO files. 

| Methods | Figure panel(s) | Comment | Relevant directory | Contributions |
 ---  | --- | --- | --- | ---
| Quality check raw data, preprocessing, visualization, and annotation | 1C, S1C, 6A, S4A | Filter data for background noise (empty droplets), filter cells and genes, exclude cells with high mitochondrial counts, ribosomal counts, and optionally ncRNA; Normalization and Dimensionality Reduction (PCA); Clustering, visualization (SPRING), and annotation | [Preprocess and Annotate Data](Preprocessing_to_annotation) | HMC



#### Analyses using annotated data

After cells are annotated by cell state, analyze the data for trends in gene expression and state characteristics. 

| Methods | Figure panel(s) | Comment | Relevant directory | Contributions |
 ---  | --- | --- | --- | ---
| Calculate and compare abundances of annotated states | 1D, 6B, S4B | Calculate abundances for each cell state; Compare abundances between treatment conditions | [Abundance Analyses](Abundance_analyses) | HMC
| Marker gene analysis, Cell cycle (G2) score analysis, Gene abundance analysis  | S1D, 6C, 6E, 6F, S4C | Compare marker genes across states, compare gene expression between treatment conditions, calculate and plot G2 score | [Marker Gene and Cell Cycle Score Analyses](Marker_gene_and_cc_score_analyses) | HMC, AMK
| Gene set enrichment analysis  | 6I | Compare genes enriched in stretch response state to similar cell states reported in literature | [Gene Set Enrichment Analysis](Gene_set_enrichment_analysis) | HMC, AMK

### SPRING Links

[Inflation perturbation dataset](https://kleintools.hms.harvard.edu/tools/springViewer_1_6_dev.html?datasets/SPRING_private/mouse_intestinal_organoids_21/murine_intestinal_organoid_inflation_perturbation_expt/all_cells)

[Sandwich and Dome culture dataset (no stretch signature annotation)](https://kleintools.hms.harvard.edu/tools/springViewer_1_6_dev.html?datasets/SPRING_private/mouse_intestinal_organoids_21/murine_intestinal_organoid_sandwich_dome_culture_merged_no_stretch_signature/all_cells)

[Sandwich and Dome culture dataset (with stretch signature annotation)](https://kleintools.hms.harvard.edu/tools/springViewer_1_6_dev.html?datasets/SPRING_private/mouse_intestinal_organoids_21/murine_intestinal_organoid_sandwich_dome_culture_merged/all_cells)

### References 

[1] Tallapragada, N. P., Cambra, H. M., Wald, T., Keough Jalbert, S., Abraham, D. M., Klein, O. D., & Klein, A. M. (2021). 
Inflation-collapse dynamics drive patterning and morphogenesis in intestinal organoids. Cell Stem Cell, [https://doi.org/10.1016/J.STEM.2021.04.002](https://pubmed.ncbi.nlm.nih.gov/33915079/) <br>
[2] Klein, A. M., Mazutis, L., Akartuna, I., Tallapragada, N., Veres, A., Li, V., Peshkin, L., Weitz, D. A., & Kirschner, M. W. (2015). Droplet barcoding for single-cell transcriptomics applied to embryonic stem cells. Cell, 161(5), 1187–1201. [https://doi.org/10.1016/j.cell.2015.04.044](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4441768/) <br>
[3] Weinreb, C., Wolock, S., & Klein, A. M. (2018). SPRING: a kinetic interface for visualizing high dimensional single-cell expression data. Bioinformatics (Oxford, England), 34(7), 1246–1248. [https://doi.org/10.1093/bioinformatics/btx792](https://pubmed.ncbi.nlm.nih.gov/29228172/) <br>
[4] Wolf, F. A., Angerer, P., & Theis, F. J. (2018). SCANPY: large-scale single-cell gene expression data analysis. Genome Biology, 19(1), 15. [https://doi.org/10.1186/s13059-017-1382-0](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-017-1382-0) <br>
[5] Wolock, S. L., Lopez, R., & Klein, A. M. (2019). Scrublet: Computational Identification of Cell Doublets in Single-Cell Transcriptomic Data. Cell systems, 8(4), 281–291.e9. [https://doi.org/10.1016/j.cels.2018.11.005](https://pubmed.ncbi.nlm.nih.gov/30954476/) <br>
[6] Cambra, Hailey M., Naren P. Tallapragada, Prabhath Mannam, David T. Breault, and Allon M. Klein. “Triple-Decker Sandwich Cultures of Intestinal Organoids for Long-Term Live Imaging, Uniform Perturbation, and Statistical Sampling.” Current Protocols 2, no. 1 (2022): e330. [https://doi.org/10.1002/cpz1.330](https://pubmed.ncbi.nlm.nih.gov/35030297/) <br>
