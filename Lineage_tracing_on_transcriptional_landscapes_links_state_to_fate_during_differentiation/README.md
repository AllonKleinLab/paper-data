### Lineage tracing on transcriptional landscapes links state to fate during differentiation (data)

This page contains data from our [manuscript](https://science.sciencemag.org/content/367/6479/eaaw3381), in which we describe a method for simultaneous lineage tracing and single-cell RNA sequencing using expressed barcodes. See below for a description of each data file and a download link. 

#### Experiment 1: In vitro differentiation time course
In this experiment, we isolated hematopoietic progenitor cells from mouse bone marrow, tagged them using the LARRY lentiviral barcode library, and then sampled them at three time points during culture (days 2, 4 and 6). We performed the experiment on two separate starting cell populations: Lin-Kit+Sca1+ (LSK) cells, and Lin-Kit+ (LK) cells. Our final dataset includes single-cell transcriptomes from each time point, as well as clonal labels and a cell type annotation for each mature cell. 

* [SPRING plot](https://kleintools.hms.harvard.edu/tools/springViewer_1_6_dev.html?cgi-bin/client_datasets/SF_all/all_combined)<br/>
This links to an online user-interface where one can plot gene expression and query enriched genes. 

* [stateFate_inVitro_normed_counts.mtx.gz](https://kleintools.hms.harvard.edu/paper_websites/state_fate2020/stateFate_inVitro_normed_counts.mtx.gz)<br/>This matrix reports the number of transcripts (UMIs) for each gene in each cell, after total-counts normalization (i.e. L1 normalization on cells). Rows represent cells and columns represent genes. There are no column or row labels. Gene names and cell metadata are provided in separate files. 

* [stateFate_inVitro_gene_names.txt.gz](https://kleintools.hms.harvard.edu/paper_websites/state_fate2020/gene_names_in_vitro.txt.gz)<br/>List of gene symbols (one per line). The rows of this file correspond to the columns of _counts_matrix_in_vitro_ (above). 

* [stateFate_inVitro_clone_matrix.mtx.gz](https://kleintools.hms.harvard.edu/paper_websites/state_fate2020/stateFate_inVitro_clone_matrix.mtx.gz)<br/>Binary matrix indicating the clonal membership of each cell. The rows of this file represent cells and correspond to the rows of _counts_matrix_in_vitro_ (above). The columns represent clones. Not every cell belongs to a clone. 

* [stateFate_inVitro_metadata.txt.gz](https://kleintools.hms.harvard.edu/paper_websites/state_fate2020/stateFate_inVitro_metadata.txt.gz)<br/>Table of cell metadata. There is a header row followed by one row for each cell. The rows of this file (after the header) correspond to the rows of _counts_matrix_in_vitro_ (above). The headers are: 
  - "Library": The sequencing library that the cell belongs to
  - "Cell barcode": Cell barcode used for demultiplexing transcripts in inDrops
  - "Time point": The time point (in days) that the cells were profiled
  - "Starting population": Either Lin-Kit+Sca1+ (LSK) for Lin-Kit+ (LK)
  - "Cell type annotation": Either "undiff" or one of ten mature cell types that appeared in culture
  - "Well": To evaluate the missing heritability of fate bias, we cultured cells in two different sets of wells after day 2. For all cells from day 2, this column will have a label of "0". For later time points, the label is either "1" or "2", corresponding to the two sets of wells. 
  - "SPRING-x/y": Coordinates for plotting cells in SPRING

* [stateFate_inVitro_neutrophil_pseudotime.txt.gz](https://kleintools.hms.harvard.edu/paper_websites/state_fate2020/stateFate_inVitro_neutrophil_pseudotime.txt.gz)<br/>Pseudotime for neutrophil trajectory cells. There is a header row followed by one row for each neutrophil trajectory cell. The first column contains a cell index (base-0 numbering) and the second column contains a pseudotime value for the cell. 


#### Experiment 2: In vivo differentiation time course
In this experiment, we isolated hematopoietic progenitor cells from mouse bone marrow, tagged them using the LARRY lentiviral barcode library, and then transplanted them into sublethally irradiated host mice. Cells were profiled on day 2 (just before transplantation) and days 9 and 16 (one and two weeks post-transplantation). See above for a detailed description of each file. 

* [SPRING plot](https://kleintools.hms.harvard.edu/tools/springViewer_1_6_dev.html?cgi-bin/client_datasets/IV_post_TP/all_combined)<br/>
* [stateFate_inVivo_normed_counts.mtx.gz](https://kleintools.hms.harvard.edu/paper_websites/state_fate2020/stateFate_inVivo_normed_counts.mtx.gz)
* [stateFate_inVivo_gene_names.txt.gz](https://kleintools.hms.harvard.edu/paper_websites/state_fate2020/stateFate_inVivo_gene_names.txt.gz)
* [stateFate_inVivo_metadata.txt.gz](https://kleintools.hms.harvard.edu/paper_websites/state_fate2020/stateFate_inVivo_metadata.txt.gz)
* [stateFate_inVivo_clone_matrix.mtx.gz](https://kleintools.hms.harvard.edu/paper_websites/state_fate2020/stateFate_inVivo_clone_matrix.mtx.gz) 

#### Experiment 3: In vivo differentiation time course
In this experiment, we isolated hematopoietic progenitor cells from mouse bone marrow, tagged them using the LARRY lentiviral barcode library, and then split the clones between several cytokine conditions for ongoing culture. Cells were profiled on day 2 (just before replating) and days 4 and 6 (i.e. after 2 and 4 days of culture with alternative cytokines).

* [SPRING plot](https://kleintools.hms.harvard.edu/tools/springViewer_1_6_dev.html?cgi-bin/client_datasets/CP2_FINAL/allMerged)<br/>
* [stateFate_cytokinePerturbation_normed_counts.mtx.gz](https://kleintools.hms.harvard.edu/paper_websites/state_fate2020/stateFate_cytokinePerturbation_normed_counts.mtx.gz)
* [stateFate_cytokinePerturbation_gene_names.txt.gz](https://kleintools.hms.harvard.edu/paper_websites/state_fate2020/stateFate_cytokinePerturbation_gene_names.txt.gz)
* [stateFate_cytokinePerturbation_metadata.txt.gz](https://kleintools.hms.harvard.edu/paper_websites/state_fate2020/stateFate_cytokinePerturbation_metadata.txt.gz)
* [stateFate_cytokinePerturbation_clone_matrix.mtx.gz](https://kleintools.hms.harvard.edu/paper_websites/state_fate2020/stateFate_cytokinePerturbation_clone_matrix.mtx.gz) 

