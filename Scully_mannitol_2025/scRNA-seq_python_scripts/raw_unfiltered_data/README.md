Download datasets directly from GEO ([GSE292926](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE292926)) into this folder, or run the `download_datasets.sh` script to do it via the command line.

To run all scRNA-seq analysis scripts, you need...

Already in this folder:
- `GSM8869530_Cr_blood_asw_barcodes.tsv`, list of barcodes called as cells by Cell Ranger
- `GSM8869531_Cr_blood_mannitol_barcodes.tsv`, list of barcodes called as cells by Cell Ranger

Need to download from GEO or with `download_datasets.sh`:
- `GSM8869530_Cr_blood_asw_raw_unfiltered_matrix.h5`, unfiltered counts matrix
- `GSM8869531_Cr_blood_mannitol_raw_unfiltered_matrix.h5`, unfiltered counts matrix
- `GSM8869530_Cr_blood_asw_molecule_info.h5`, file with info about reads
- `GSM8869531_Cr_blood_mannitol_molecule_info.h5`, file with info about reads
