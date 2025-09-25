There are three python scripts in this folder:

#### 1. Filtering and preprocessing, `filtering_and_preprocessing.py`
This script takes in the raw, unfiltered counts matrices, and outputs the filtered and preprocessed dataset in an .h5ad file. This script assumes the unfiltered counts matrices, [downloadable from GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE292926), are in the subfolder `raw_unfiltered_data/` (see that folder's README file for scripts to download this data using bash). The processed data file is saved in a subfolder which is created by the script.

#### 2. Data quality metrics, `data_quality_metrics.py`
This script takes in the processed counts matrix and [Cell Ranger Molecule Info HDF5 file](https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/outputs/cr-outputs-molecule-info), and outputs a csv file summarizing several data quality metrics, which also appear in [TBD supplemental table]. This script assumes the molecule info files are in the subfolder `raw_unfiltered_data/` (see that folder's README file for scripts to download this data using bash), and that `filtering_and_preprocessing.py` has been run. The csv files are saved in a subfolder which is created by the script.

#### 3. Figure 2 plots, `figure2_plots.py`
This script uses the unfiltered and the processed data and generates many of the plots in Figure 2 of the paper. This script assumes that the unfiltered counts matrices are in the subfolder `raw_unfiltered_data/`, and that both `filtering_and_preprocessing.py` and `data_quality_metrics.py` have been run. The plots are saved in a subfolder created by this script.

#### 4. Helper functions, `helper_functions.py`
This script contains various helper functions used by both of the above scripts.
