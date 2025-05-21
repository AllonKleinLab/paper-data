The scripts in this folder were used to assess gene coexpression conservation, following methods introduced by [Crow et al. 2022](https://dx.doi.org/10.1093/nar/gkac276). See methods sections "Co-expression Conservation".

The files are:
- `coexpr_cons_Hs_vs_Cr.py` is the python script that calculates coexpression AUROC scores, as described in Crow et al. 2022. The script outputs a folder `coexpr_cons_Hs_vs_Cr_output/` containing:
    - `coexpr_conservation_all.csv`, a CSV containing AUROC scores. Note that the column `raw_auroc` is used in the paper (Supplementary Fig. 9), and that the column `coexpr_conservation` contains a normalized coexpression conservation score described in [Suresh et al. 2023](http://dx.doi.org/10.1038/s41559-023-02186-7) which was not used in our study. See the section "Note on scaled AUROC score" for more detail.
    - `network1_gene_list.npy` and `network1_weights.npy`, gene coexpression network connectivity matrices
    - PDFs with example AUROC curves for select gene pairs
- `tf_coexpr_Hs_vs_Cr.py` is a python script that looks at coexpression conservation scores for a select set of hematopoiesis-related TFs. The script saves outputs in a folder `tf_coexpr_Hs_vs_Cr_output/` containing graphs in Fig. 6f-h and Supplementary Fig. 9.
- `coexpr_cons_Hs_vs_Dr.py` and `tf_coexpr_Hs_vs_Dr.py` are python scripts which repeat these analyses for the human vs. zebrafish comparison.
- `single_cell_analysis_env.txt` contains information on the conda environment used, saved to a text file (by running `conda list --explicit > single_cell_analysis_env.txt`).

Note that the functions for actually calculating coexpression conservation AUROCs are in the helper_functions folder: `Scully_Ciona_blood_2025/helper_functions/gillis_style_coexpression_hf.py`. These functions are imported by the python scripts in this folder.

## Note on scaled AUROC score
The coexpression conservation score calculation was introduced in a few papers from Jesse Gillis's group. [Crow et al. 2022](https://dx.doi.org/10.1093/nar/gkac276) introduces a coexpression conservation score based on calculation of an AUROC for N-to-M gene homology. [Suresh et al. 2023](http://dx.doi.org/10.1038/s41559-023-02186-7) adds an additional step to get a scaled coexpression conservation score:
1. Calculate AUROC scores for all gene pairs, not just homologous pairs, then
2. Calculate a final coexpression similarity score for each homologous gene pair as ranked with respect to the AUROCs for all possible gene pairs. Normalize to 1 to get the final score.

In our paper, we only use the AUROC score as introduced by Crow et al., but the code supports an approximation of the Suresh et al. ranking score. To improve script run times, we approximate the Suresh et al. ranking score by calculating AUROCs for a large number of non-homologous gene pairs (rather than for all possible pairs). The function `calculate_coexpression_conservation_complex_orthology` takes an input `num_random_pairs`, which sets the number of non-homologous gene pairs to look at. If this parameter is not specified, the functions default to setting `num_random_pairs` to be 100 times the number of homologous gene pairs, such that homologous gene pairs make up only 1% of the final AUROC distribution.
