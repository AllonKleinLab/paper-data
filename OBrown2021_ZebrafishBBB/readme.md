### This repository holds the code for the analysis of the single-cell data collected in:
O'Brown et. al. The secreted neuronal signal Spock1 regulates the blood-brain barrier (2023)

### The annotated counts matrix (GSM7208220_210902_PATBROWN_final.h5ad) can be found here:
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM7208220

### Description of the Jupyter notebooks:
- STEP1_Filter_Dataset.ipynb: Annotation of the full single-cell dataset for both the WT and hm41 mutant genotypes. Filter out non-neuronal cell types.
- STEP2_Full_Dataset_Plots.ipynb: Visualize marker genes in the filtered single-cell dataset. Perform differential gene-expression (DGE) analysis for WT ahd hm41 genotypes.
- STEP3_Subcluster_Vasculature: Further filter the data to analyze cells in the vasculature
- STEP4_Vasculature_Plots: Plot Endothelial and Pericyte marker genes in the vasculature. Perform vasculature-specific DGE analysis.
