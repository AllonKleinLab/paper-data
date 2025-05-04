The scripts in this folder were used to run Gene Ontology (GO) enrichment analysis. See methods section "Gene Ontology Enrichment Analysis".

The files are:
- `go_enrichment_analysis.py` is the python script that performs GO analysis. The script outputs a folder `go_enrichment_analysis_output/` containing tsv files with results from Fisher's exact test of GO term enrichment for each cell state's enriched DEGs.
- `go_enrichment_env.txt` contains information on the conda environment used, saved to a text file (by running `conda list --explicit > go_enrichment_env.txt`).
