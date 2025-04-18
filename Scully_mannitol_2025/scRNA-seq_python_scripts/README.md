There are three python scripts in this folder:

#### 1. Filtering and preprocessing, `filtering_and_preprocessing.py`
This script takes in the raw, unfiltered counts matrices, and outputs the filtered and preprocessed dataset in an .h5ad file. This script assumes the unfiltered counts matrices, [downloadable from GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8869531), are in a subfolder called `raw_unfiltered_data/`. The processed data file is saved in a subfolder which is created by the script.

#### 2. Figure 2 plots, `figure2_plots.py`
This script uses the unfiltered and the processed data and generates many of the plots in Figure 2 of the paper. This script assumes that the unfiltered counts matrices are in a subfolder called `raw_unfiltered_data/`, and that the filtered & processed data is in the subfolder generated by the filtering and preprocessing script. The plots are saved in a subfolder created by this script.

#### 3. Helper functions, `helper_functions.py`
This script contains various helper functions used by both of the above scripts.
