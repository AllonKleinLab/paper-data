The scripts in this folder were used to run Gene Ontology (GO) enrichment analysis. See methods section "Gene Ontology Enrichment Analysis".

The files are:

_Enrichment analysis on curated gene sets:_
- `curated_gene_lists.csv` is part of Table S3, listing the genes in the manually curated gene sets.
- `curated_gene_sets_fishers.py` is the python script that performs Fisher's exact test enrichment analysis on the manually curated gene sets. The script outputs a folder `kegg_pathway_gsea_output/` containing tsv files with full results.
- `curated_gene_sets_gsea.py` is the python script that performs GSEA on the manually curated gene sets. The script outputs a folder `kegg_pathway_gsea_output/` containing tsv files with full results.

_Enrichment analysis on KEGG Pathways:_
- `kegg_pathway_fishers.py` is the python script that performs Fisher's exact test enrichment analysis on KEGG pathways. The script outputs a folder `kegg_pathway_gsea_output/` containing tsv files with full results.
- `kegg_pathway_gsea.py` is the python script that performs GSEA on KEGG pathways. The script outputs a folder `kegg_pathway_gsea_output/` containing tsv files with full results.

_Supporting files:_
- `gene_enrichment_env.yml` and `gene_enrichment_env.txt` contain information on the conda environment used, saved to a text file (by running `conda env export --no-builds > gene_enrichment_env.yml` and `conda list --explicit > gene_enrichment_env.txt`). The .txt files is for MacOS only.
