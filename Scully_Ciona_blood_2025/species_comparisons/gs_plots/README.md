The scripts in this folder were used to make generalized Sankey (GS) plots to compare cell states across species, and to calculate Jaccard similarity scores. See methods sections "Generalized Sankey (GS) Plots" and "Jaccard Similarity Score".

The files are:
- `gs_plots_Hs_vs_Cr.py` is the python script that both makes GS plots and calculates Jaccard similarity scores for the human vs. _C. robusta_ comparison. The script outputs a folder `gs_plots_Hs_vs_Cr_output/` containing GS plots (with both humand and _C. robusta_ cell states at the top) and Jaccard similarity matrices.
- `gs_plots_Hs_vs_Dr.py` is the python script which repeats these anlyses for the human vs. zebrafish comparison. The script saves outputs in a folder `gs_plots_Hs_vs_Dr_output/` containing containing GS plots (with both humand and zebrafish cell states at the top) and Jaccard similarity matrices.
- `single_cell_analysis_env.yml` and `single_cell_analysis_env.txt` contain information on the conda environment used, saved to a text file (by running `conda env export --no-builds > single_cell_analysis_env.yml` and `conda list --explicit > single_cell_analysis_env.txt`). The .txt files is for MacOS only.

Note that the functions for actually making the GS plots are in the helper_functions folder: `Scully_Ciona_blood_2025/helper_functions/gs_plot_helper_functions.py`. These functions are imported by the python scripts in this folder.
