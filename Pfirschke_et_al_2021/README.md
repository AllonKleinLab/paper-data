# ScRNAseq-related analyses for Pfirschke et al. 2021 [1]

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
Large data files used by the code in this repository are stored on [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE161771).
| Description | Filename(s) | 
 ---  | --- 
| Raw count data before filtering | <library_name>.barcodes.tsv, <library_name>.counts.mtx, <libraray_name>.genes.tsv ([GSE161771](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE161771)|
| Anndata object combining all libraries after mitochondrial fraction and total count filtering | mito_total_counts_filt_raw_27563x40930_200517_10h29.h5ad |
| Cell annotation file (i.e. adata.obs) | obs_info_27563x32_201025_14h44.tsv (or .npz [here](data/obs_info_27563x32_201025_14h44.npz)*) |

### Methods

#### Contributors to this repo
Angela E. Zou (AZ)  
Marius Messemaker (MM)  
Nicolas A. Gort-Freitas (NAGF)  
Rapolas Zilionis (RZ) 

#### From reads to counts
The [indrop.py](https://github.com/indrops) pipeline was used for obtain cells x genes matrices. Yaml files detailing the parameters used can be found [here](yaml_files_for_indrops_pipeline).


#### Filtering, doublet removal, visualization and annotation
Notebooks starting with "part*" focus on the analysis steps up to filtered annotated data.
Data filtering involved repeating cell visualization several times. For example, drawing and clustering
a kNN graph of T cells only revealed two distinct double populations within T cells, the removal of which
required repeating the visualization of all CD45+ cells.

| Methods | Figure panel(s) | Comment | Relevant notebooks | Contributions |
 ---  | --- | --- | --- | ---
| Filtering on total counts and mitochondrial fraction | NA | Filter data and combine individual mtx files into one AnnData object | [part1_read_and_filter_data.ipynb](methods_clean_visual_annot/part1_read_and_filter_data.ipynb) | RZ
| Cell annotation using a Bayesian classifier | NA | Classify single cell transcriptomes using a Bayesian classifier as described previously [2,3,4] | [part3_classify_cell_by_published_profiles.ipynb](methods_clean_visual_annot/part3_classify_cell_by_published_profiles.ipynb) | RZ
| Initial visualization of the data using SPRING [5]| [Interactive explorer online](https://kleintools.hms.harvard.edu/tools/springViewer_1_6_dev.html?datasets/pittet/mouse/CSF1Ri/all_cells_w_dblt_1000umi) | Visualize data using a force-directed layout of a kNN graph | [part3c_spring_plot_1000_umi.ipynb](methods_clean_visual_annot/part3c_spring_plot_1000_umi.ipynb) | RZ
| Doublet detection | NA | Obtaining doublet scores using Scrublet [6]. Modified version using precomputed PCA: for consistency, the same PCA-transformed data as for drawing the kNN graph in SPRING visualization is used| [part4_detect_doublets_in_each_emulsion_precomputed_PCA.ipynb](methods_clean_visual_annot/part4_detect_doublets_in_each_emulsion_precomputed_PCA.ipynb) | RZ
| Spectral clustering of the SPRING plot | NA | Divide the kNN graph into a predefined number of clusters| [part5_divide_graph_with_doublets_into_clusters.ipynb](methods_clean_visual_annot/part5_divide_graph_with_doublets_into_clusters.ipynb)
| Identify cells to exclude | NA | Identify and flag for removal clusters enriched in doublets, Krt8 (non-immune), and hemoglobin genes (erythroid) | [part6_decide_which_cells_to_exclude.ipynb](methods_clean_visual_annot/part6_decide_which_cells_to_exclude.ipynb) | RZ |
| Repeat data visualization and spectral clustering after cleanup | [Interactive explorer online](https://kleintools.hms.harvard.edu/tools/springViewer_1_6_dev.html?datasets/pittet/mouse/CSF1Ri/all_cells_clean_iter1_refCSF1Ri) | Visualization exclude cells flagged in the previous notebook. Eigenvectors for PCA calculated on Csf1Ri condition cells only (as a method for batch correction). | [part7_sping_plot_main_iter1.ipynb](methods_clean_visual_annot/part7_sping_plot_main_iter1.ipynb), [part8_divide_graph_into_clusters.ipynb](methods_clean_visual_annot/part8_divide_graph_into_clusters.ipynb) | RZ |
| Identify more cells to exclude | NA | Repeating the visualization and clustering after initial clean up revealed more distinct doublet clusters | [part8_divide_graph_into_clusters.ipynb](methods_clean_visual_annot/part8_divide_graph_into_clusters.ipynb) | RZ |
| Data visualization and clustering after second round of cleanup | [Interactive explorer online](https://kleintools.hms.harvard.edu/tools/springViewer_1_6_dev.html?datasets/pittet/mouse/CSF1Ri/all_cells_clean_iter2) | Second iteration of plotting cleaned up data | [part10_sping_plot_main_iter2.ipynb](methods_clean_visual_annot/part10_sping_plot_main_iter2.ipynb), [part11_divide_graph_into_clusters.ipynb](methods_clean_visual_annot/part11_divide_graph_into_clusters.ipynb)| RZ |
| Cluster annotation using Bayesian classifier results | [Interactive explorer online](https://kleintools.hms.harvard.edu/tools/springViewer_1_6_dev.html?datasets/pittet/mouse/CSF1Ri/all_cells_clean_iter2), colortrack "*population" | Clusters are labeled after the dominant classification result obtained for individual cells. Ambiguous cases are reviewed manually. | [part12_define_populations_based_on_classifier_results.ipynb](methods_clean_visual_annot/part12_define_populations_based_on_classifier_results.ipynb) | RZ |
| Visualize, cluster, and annotated T cells separately | [Interactive explorer online](https://kleintools.hms.harvard.edu/tools/springViewer_1_6_dev.html?datasets/pittet/mouse/CSF1Ri/T_cells_only)| Resolving T cell sub-populations required subclustering T cells. Annotation is based on interactively exploring the resulting SPRING plot, and comparing cluster-enriched gene expression to known marker genes | [part13_sping_plot_of_T_only.ipynb](methods_clean_visual_annot/part13_sping_plot_of_T_only.ipynb), [part14_divide_T_cell_graph_into_clusters.ipynb](methods_clean_visual_annot/part14_divide_T_cell_graph_into_clusters.ipynb), [part15_define_T_subsets.ipynb](methods_clean_visual_annot/part15_define_T_subsets.ipynb) | RZ |
| Final data visualization, clustering, and annotation after removing doublet clusters within T cells | Fig. 2B-C; 4E; 5D and more,  [Interactive explorer online](https://kleintools.hms.harvard.edu/tools/springViewer_1_6_dev.html?datasets/pittet/mouse/CSF1Ri/all_Cd45pos_cells)  | Visualization and annotation of T cells only in the previous step revealed two distinct doublet clusters. Visualization of all CD45+ cells was repeated with the further cleaned up data | [part16_sping_plot_main_iter3.ipynb](methods_clean_visual_annot/part16_sping_plot_main_iter3.ipynb), [part17_divide_graph_into_clusters.ipynb](methods_clean_visual_annot/part17_divide_graph_into_clusters.ipynb), [part18_define_populations_based_on_classifier_results.ipynb](methods_clean_visual_annot/part18_define_populations_based_on_classifier_results.ipynb) | RZ |
| Final annotation of T cells | [Interactive explorer online](https://kleintools.hms.harvard.edu/tools/springViewer_1_6_dev.html?datasets/pittet/mouse/CSF1Ri/T_cells_only_iter2) | Repeat T cell annotation after removing T cell doublets | [part19_sping_plot_of_T_only.ipynb](methods_clean_visual_annot/part19_sping_plot_of_T_only.ipynb), [part20_divide_T_cell_graph_into_clusters.ipynb](methods_clean_visual_annot/part20_divide_T_cell_graph_into_clusters.ipynb), [part21_define_T_subsets.ipynb](methods_clean_visual_annot/part21_define_T_subsets.ipynb) | RZ |
| Clean up annotations | Fig. 2B | Tidy up cell annotations in adata.obs, merge minor populations to main types (e.g. DC1, DC2, DC3, pDC collectively are DCs) | [part22_cleanup_labels.ipynb](methods_clean_visual_annot/part22_cleanup_labels.ipynb) | RZ |

	
#### Analyses using annotated data
| Methods | Figure panel(s) | Comment | Relevant notebooks | Contributions |
 ---  | --- | --- | --- | ---
| Example notebook for loading annotated data and xy coordinates | E.g. Fig. 2B | This notebooks uses a few examples to guide anyone interested through how the annotated data is organized | [example_load_data_plot_something.ipynb](methods_post_annotation/example_load_data_plot_something.ipynb) | RZ |
| Plot a subset of cells from the main SPRING plot, color by gene expression | Multiple figures, the motif of Fig. 4E | Load xy coordinates, select a subset of cells, color by gene expression or population annotation | [Colored_SPRING_plots.ipynb](methods_post_annotation/Colored_SPRING_plots.ipynb) | RZ |
| Challenge annotations by plotting a heatmap of previously identified marker genes | Fig. 2D | Recreate marker gene heatmap from previous study [2] (same gene order) but using the newly defined cell populations | [Annotation_challenging_marker_gene_heatmaps.ipynb](methods_post_annotation/Annotation_challenging_marker_gene_heatmaps.ipynb) | RZ |
| State %CD45 abundance, Arrow gene-expression change, and Differential Expression Analysis volcano Figures | Figs. 2E-F, S2A, S2C-E, S2G-H, S3, 4G, 5F, S5H, S5J-K, S6E-J  | |  [Abundance_and_expression_change_analysis.ipynb](methods_post_annotation/Abundance_and_expression_change_analysis.ipynb) | MM |
| Make dot plots of relative gene expression and % cells expressing genes | Figs 2G, 3A, S4A |  | [for-github_dotplot.ipynb](methods_post_annotation/for-github_dotplot.ipynb) | AZ |
| Perform GSEA on GO:BP terms, make scatterplot of enriched immune activation-related terms | Fig. 2H |  | [for-github_fgsea-scatterplot.ipynb](methods_post_annotation/for-github_fgsea-scatterplot.ipynb) | AZ |
| Heatmap of scores for selected GO:BP terms in MoMac cells | Fig. 2I |  | [for-github_scored-pathway-HM.ipynb](methods_post_annotation/for-github_scored-pathway-HM.ipynb) | AZ |
| Make circos plots for differentially expressed and immune activating/inhibitory interactions | Figs. 3B, 5B |  | [for-github_cell-cell-comms_filter+circos.ipynb](methods_post_annotation/for-github_cell-cell-comms_filter+circos.ipynb) | AZ |
| Make heatmaps depicting selected ligand-receptor interactions | Figs. 3C, S4C, 5C, S6C |  | [for-github_cell-cell-comms_make-HMs.ipynb](methods_post_annotation/for-github_cell-cell-comms_make-HMs.ipynb) | AZ |
| Fold change with respect to the median across states compared (relative expression); Pearson's r correlation; Linear regression | 3D, S4D, 5E, S6D | In this notebook, we compare the relative expression of highlighted ligands & receptors in DCs, NKs, T cells, and Monocyte/Macrophages in non-small cell lung cancer patients (Zilionis et al., 2019) and in vehicle-treated mice to support a cross-species analogy in the behavior of immunity. | [heatmaps_scatter_human_mouse.ipynb](methods_post_annotation/heatmaps_scatter_human_mouse.ipynb) | NAGF |


### References   
[1] To be updated after publication  
[2] Zilionis R, Engblom C, Pfirschke C, et al. Single-Cell Transcriptomics of Human and Mouse Lung Cancers Reveals Conserved Myeloid Populations across Individuals and Species. Immunity. 2019;50(5):1317-1334.e10. doi:10.1016/j.immuni.2019.03.009  
[3] Zemmour, D., Zilionis, R., Kiner, E., Klein, A.M., Mathis, D., and Benoist, C. (2018). Single-cell gene expression reveals a landscape of regulatory T cell phenotypes shaped by the TCR. Nat. Immunol. 19, 291–301  
[4] Jaitin, D.A., Kenigsberg, E., Keren-Shaul, H., Elefant, N., Paul, F., Zaretsky, I., Mildner, A., Cohen, N., Jung, S., Tanay, A., and Amit, I. (2014). Massively par- allel single-cell RNA-seq for marker-free decomposition of tissues into cell types. Science 343, 776–779.  
[5] Weinreb C, Wolock S, Klein AM. SPRING: a kinetic interface for visualizing high dimensional single-cell expression data. Bioinformatics. 2018;34(7):1246-1248. doi:10.1093/bioinformatics/btx792  
[6] Wolock SL, Lopez R, Klein AM. Scrublet: Computational Identification of Cell Doublets in Single-Cell Transcriptomic Data. Cell Syst. 2019;8(4):281-291.e9. doi:10.1016/j.cels.2018.11.005  
