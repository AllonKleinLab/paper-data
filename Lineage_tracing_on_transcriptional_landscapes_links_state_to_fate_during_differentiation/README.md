### Lineage tracing on transcriptional landscapes links state to fate during differentiation (data)

This page contains data from our [manuscript](https://www.biorxiv.org/content/early/2018/11/11/467886), in which we describe a method for simultaneous lineage tracing and single-cell RNA sequencing using expressed barcodes. See below for a description of each data file and a download link. 

#### Experiment 1: In vitro differentiation time course
In this experiment, we isolated hematopoietic progenitor cells from mouse bone marrow, tagged them using the LARRY lentiviral barcode library, and then sampled them at three time points during culture (days 2, 4 and 6). We performed the experiment on two separate starting cell populations: Lin-Kit+Sca1+ (LSK) cells, and Lin-Kit+ (LK) cells. Our final dataset includes single-cell transcriptomes from each time point, as well as clonal labels and a cell type annotation for each mature cell. 

* [SPRING plot](https://kleintools.hms.harvard.edu/tools/springViewer_1_6_dev.html?cgi-bin/client_datasets/SF_all/all_combined)<br/>
This links to an online user-interface where one can plot gene expression and query enriched genes. 

* [counts_matrix_in_vitro.mtx.gz](https://kleintools.hms.harvard.edu/paper_websites/state_fate2018/counts_matrix_in_vitro.mtx.gz) / [counts_matrix_in_vitro.npz.gz](https://kleintools.hms.harvard.edu/paper_websites/state_fate2018/counts_matrix_in_vitro.npz.gz)<br/>This matrix reports the number of transcripts (UMIs) for each gene in each cell, after total-counts normalization (i.e. L1 normalization on cells). Rows represent cells and columns represent genes. There are no column or row labels. Gene names and cell metadata are provided in separate files. 

* [gene_names_in_vitro.txt](www.google.com)<br/>List of gene symbols (one per line). The rows of this file correspond to the columns of _counts_matrix_in_vitro_ (above). 

* [coordinates_in_vitro.txt]()<br/>Coordinates for plotting cells (generated using SPRING). This file has two columns representing the x-coordinate and y-coordinate of each cell respectively. The rows of this file correspond to the rows of _counts_matrix_in_vitro_ (above).

* [clone_annotation_in_vitro.mtx](www.google.com) / [clone_annotation_in_vitro.npz](www.google.com)<br/>Binary matrix indicating the clonal membership of each cell. The rows of this file represent cells and correspond to the rows of _counts_matrix_in_vitro_ (above). The columns represent clones. Not every cell belongs to a clone. 

* [cell_metadata_in_vitro.txt](www.google.com)<br/>Table of cell metadata. There is a header row followed by one row for each cell. The rows of this file (after the header) correspond to the rows of _counts_matrix_in_vitro_ (above). The headers are: 
  - "Time point": The time point (in days) that the cells were profiled
  - "Population": Either Lin-Kit+Sca1+ (LSK) for Lin-Kit+ (LK)
  - "Annotation": Either "undiff" or one of ten mature cell types that appeared in culture
  - "Well": To evaluate the missing heritability of fate bias, we cultured cells in two different sets of wells after day 2. For all cells from day 2, this column will have a label of "0". For later time points, the label is either "1" or "2", corresponding to the two sets of wells. 

