"""
This script fits a Generalized Linear Model (GLM) to analyze gene expression data for multiple cell states 
(EEP, CEP-1, and CEP-2). The model aims to assess the effects of IL-17, EPO, and the organ source on gene expression, 
while accounting for their potential interactions.

## Key Concepts:

1. **Generalized Linear Model (GLM)**:
    - A GLM is a flexible extension of ordinary linear regression that allows the dependent variable (gene expression) 
      to follow different distributions (here we assume Gaussian for log-transformed data).
    - The GLM will fit coefficients to explain how different factors (IL-17, EPO, Organ) and their interactions 
      contribute to changes in gene expression.

2. **Model Terms**:
    - **IL-17**: Binary indicator of whether the treatment includes IL-17 stimulation.
    - **EPO**: Binary indicator of whether the treatment includes EPO stimulation.
    - **Organ**: A categorical variable that represents the organ source of the cells, with `Organ_SP` representing one organ 
      (e.g., Spleen) and the other organs serving as reference.
    - **Interactions**: The model includes interactions between **IL-17**, **EPO**, and **Organ**, allowing the effects of 
      IL-17 and EPO to vary depending on the organ.

## Model Comparisons (Likelihood Ratio Tests):

Each model comparison uses a **Likelihood Ratio Test (LRT)** to compare the full model (with all terms) to reduced models
that omit certain terms. The LRT helps determine whether removing a specific term or interaction significantly worsens
the model's fit, meaning that the term or interaction is important for explaining the gene expression data.

### Model Definitions:

- **Full Model**: 
  `expression ~ Organ_SP * IL17 * EPO`
  - This model includes all main effects (IL-17, EPO, and Organ) and all possible interactions between them.
  - The most complex model, it tests whether gene expression is influenced by all combinations of these factors.

- **Reduced Models**: 
  We fit several reduced models, each omitting specific terms or interactions, and compare them to the full model:

  - **Model 0**: 
    `expression ~ IL17 * EPO`
    - No organ effects. This model tests the influence of **IL-17 and EPO** but assumes that the **organ** does not influence gene expression.
    
  - **Model 1**: 
    `expression ~ Organ_SP + IL17 + EPO + Organ_SP:EPO + EPO:IL17`
    - This model omits the **Organ_SP:IL17 interaction**, testing whether the effect of **IL-17** is consistent across organs.
    
  - **Model 2**: 
    `expression ~ Organ_SP + IL17 + Organ_SP:IL17`
    - This model excludes **EPO**, testing whether **IL-17** and its interaction with **Organ** drive gene expression changes, 
      while ignoring the effects of **EPO**.
    
  - **Model 3**: 
    `expression ~ Organ_SP + EPO + Organ_SP:EPO`
    - This model excludes **IL-17**, testing the effects of **EPO** and its interaction with **Organ**, ignoring the role of IL-17.
    
  - **Model 4**: 
    `expression ~ Organ_SP + EPO + IL17 + Organ_SP:EPO + Organ_SP:IL17`
    - This model excludes the interaction between **EPO and IL-17**, testing whether the effects of these two factors are independent.
    
  - **Model 5**: 
    `expression ~ Organ_SP + EPO + IL17 + EPO:IL17 + Organ_SP:EPO + Organ_SP:IL17`
    - This model omits the **Organ_SP:EPO:IL17** interaction, testing whether the combined effect of **EPO and IL-17** 
      is consistent across organs.
      
  - **Model 6**: 
    `expression ~ Organ_SP + IL17 + EPO + Organ_SP:IL17 + EPO:IL17`
    - This model omits the **Organ_SP:EPO interaction**, testing whether the effect of **EPO** is consistent across organs.


## Procedure:

1. **Fit the GLM**: For each cell state (EEP, CEP-1, CEP-2), we fit the full model and reduced models for each gene.
2. **Likelihood Ratio Test (LRT)**: For each gene, we compute the LRT statistic, which compares the fit of each reduced model 
   to the full model. A significant LRT p-value indicates that the removed terms (or interactions) are important for explaining gene expression.
3. **Multiple Hypothesis Correction**: To adjust for multiple comparisons, we apply the Benjamini-Hochberg (FDR) method to the LRT p-values.
4. **Fold Changes**: For significant terms, we calculate the fold change by exponentiating the coefficients to interpret the effect of each factor on gene expression.
5. **Saving Results**: The results for each cell state are saved to separate CSV files (`results_{cell_state}.csv`) for further analysis.

"""

from multiprocessing import Pool
from tqdm import tqdm
import statsmodels.api as sm
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
import scanpy as sc

# Load the AnnData object
path = '/n/groups/klein/qiu/exp_0070_analysis/Data_Files/240319_adjusted_annotations.h5ad'
adata = sc.read(path)

# List of sufficient cell states you want to run through
cell_states = ['EEP', 'CEP-1', 'CEP-2']

# Subsetting relevant observables
metadata_vals = ['stimulation', 'HT_organ', 'HT_animal', 'HT_stimulation', 'cell_annotations']
adata = adata[adata.obs['cell_annotations'].isin(cell_states) & adata.obs['Valid_Cell'], :]
adata.obs = adata.obs.loc[:, metadata_vals]
adata.X = adata.layers['raw'].toarray()  # Convert sparse matrix to array

# Preprocessing: Filter genes with expression in at least 10 cells
sc.pp.filter_genes(adata, min_cells=10)

# Normalize the data and log-transform it
sc.pp.normalize_total(adata, target_sum=1e6)
adata.layers['tpm'] = adata.X
sc.pp.log1p(adata, layer='tpm')

# Prepare metadata for modeling
metadata = adata.obs.reset_index().rename(columns={
    'HT_organ': 'Organ',
    'stimulation': 'Lane',
    'index': 'SampleID',
    'HT_animal': 'AnimalID',
    'HT_stimulation': 'Treatment',
    'cell_annotations': 'CellType'
}).loc[:, ['SampleID', 'AnimalID', 'Organ', 'Treatment', 'CellType', 'Lane']]
metadata.index = metadata.SampleID

# Create binary indicators for IL17 and EPO
metadata['IL17'] = metadata['Treatment'].apply(lambda x: 1 if 'IL17' in x else 0)
metadata['EPO'] = metadata['Treatment'].apply(lambda x: 1 if 'EPO' in x else 0)
metadata['Organ'] = metadata['Organ'].astype('category')

# One-hot encode Organ column (categorical)
metadata = pd.get_dummies(metadata, columns=['Organ'], drop_first=True)

# Define the formulas for the full model and reduced models
formulas = {
    'full_model': 'expression ~ Organ_SP * IL17 * EPO',
    'model_0': 'expression ~ IL17 * EPO',  # No organ difference
    'model_1': 'expression ~ Organ_SP + IL17 + EPO + Organ_SP:EPO + EPO:IL17',  # No organ difference in IL17 action
    'model_2': 'expression ~ Organ_SP + IL17 + Organ_SP:IL17',  # EPO doesn’t matter
    'model_3': 'expression ~ Organ_SP + EPO + Organ_SP:EPO',  # IL17 doesn’t matter
    'model_4': 'expression ~ Organ_SP + EPO + IL17 + Organ_SP:EPO + Organ_SP:IL17',  # IL17:EPO interaction doesn’t matter
    'model_5': 'expression ~ Organ_SP + EPO + IL17 + EPO:IL17 + Organ_SP:EPO + Organ_SP:IL17',  # IL17:EPO interaction ONLY doesn’t vary between organs
    'model_6': 'expression ~ Organ_SP + IL17 + EPO + Organ_SP:IL17 + EPO:IL17',  # No organ difference in EPO action
}

# Define which terms are excluded in each reduced model
excluded_terms = {
    'model_0': ['Organ_SP', 'Organ_SP:IL17', 'Organ_SP:EPO', 'Organ_SP:EPO:IL17'],
    'model_1': ['Organ_SP:IL17'],
    'model_2': ['EPO', 'EPO:IL17', 'Organ_SP:EPO', 'Organ_SP:EPO:IL17'],
    'model_3': ['IL17', 'IL17:EPO', 'Organ_SP:IL17', 'Organ_SP:EPO:IL17'],
    'model_4': ['IL17:EPO'],
    'model_5': ['Organ_SP:EPO:IL17'],
    'model_6': ['Organ_SP:EPO'],}

adata.obs = metadata

# Function to process a single cell state (this will be parallelized)
def process_cell_state(cell_state):
    print(f"Processing cell state: {cell_state}")
    
    # Subset the data for the current cell state
    current_adata = adata[adata.obs['CellType'] == cell_state]
    expression_data = current_adata.to_df()
    metadata_subset = metadata.loc[expression_data.index]

    # Prepare a list to store results
    results = []

    # Loop through each gene
    for gene in tqdm(expression_data.columns, desc=f'Processing Genes for {cell_state}'):
        data = metadata_subset.copy()
        gene_expression = expression_data[gene].values

        # Skip genes with low or zero expression
        if np.mean(gene_expression) < 0.1 or np.var(gene_expression) == 0:
            continue

        data['expression'] = gene_expression

        try:
            # Fit the full model
            full_model = sm.GLM.from_formula(formulas['full_model'], data, family=sm.families.Gaussian())
            full_result = full_model.fit()

            # Get the log-likelihood of the full model
            full_log_likelihood = full_result.llf

            gene_result = {'Gene': gene}

            # Store full model coefficients and p-values
            for term in full_model.exog_names:
                coef_key = f'{term}_coef'
                pval_key = f'{term}_pval'
                gene_result[coef_key] = full_result.params.get(term, np.nan)
                gene_result[pval_key] = full_result.pvalues.get(term, np.nan)

            # Calculate fold changes for relevant terms
            if 'Organ_SP_coef' in gene_result:
                gene_result['Organ_SP_fold_change'] = np.exp(gene_result['Organ_SP_coef'])
            if 'IL17_coef' in gene_result:
                gene_result['IL17_fold_change'] = np.exp(gene_result['IL17_coef'])
            if 'EPO_coef' in gene_result:
                gene_result['EPO_fold_change'] = np.exp(gene_result['EPO_coef'])
            if 'EPO:IL17_coef' in gene_result:
                gene_result['EPO_IL17_fold_change'] = np.exp(gene_result['EPO:IL17_coef'])

            # Fit reduced models and compute LRTs
            for model_name, formula in formulas.items():
                if model_name == 'full_model':
                    continue
                
                reduced_model = sm.GLM.from_formula(formula, data, family=sm.families.Gaussian())
                reduced_result = reduced_model.fit()

                # Log-likelihood of the reduced model
                reduced_log_likelihood = reduced_result.llf

                # Compute LRT statistic
                lr_stat = 2 * (full_log_likelihood - reduced_log_likelihood)
                df_diff = full_model.df_model - reduced_model.df_model
                p_value_lrt = stats.chi2.sf(lr_stat, df_diff)

                # Store the LRT statistic and p-value
                gene_result[f'{model_name}_lrt_stat'] = lr_stat
                gene_result[f'{model_name}_lrt_pval'] = p_value_lrt
            
            # Summarize overall effect of the excluded terms in this model
                excluded_coef_sum = 0
                excluded_coef_count = 0
                for term in excluded_terms.get(model_name, []):
                    coef_key = f'{term}_coef'
                    if coef_key in gene_result:
                        excluded_coef_sum += gene_result.get(coef_key, 0)
                        excluded_coef_count += 1
                
                if excluded_coef_count > 0:
                    excluded_coef_avg = excluded_coef_sum / excluded_coef_count
                else:
                    excluded_coef_avg = 0
                
                # Determine overall direction (positive/negative/no effect)
                if excluded_coef_avg > 0:
                    gene_result[f'{model_name}_overall_effect'] = 'positive'
                elif excluded_coef_avg < 0:
                    gene_result[f'{model_name}_overall_effect'] = 'negative'
                else:
                    gene_result[f'{model_name}_overall_effect'] = 'no effect'

                # Calculate the fold change for the average excluded effect
                gene_result[f'{model_name}_overall_fold_change'] = np.exp(excluded_coef_avg)

        except Exception as e:
            print(f"Model fitting failed for gene {gene}: {e}")
            gene_result = {'Gene': gene}
            for model_name in formulas.keys():
                gene_result[f'{model_name}_lrt_stat'] = np.nan
                gene_result[f'{model_name}_lrt_pval'] = np.nan
    
    

        results.append(gene_result)

    # Convert the results to a DataFrame
    results_df = pd.DataFrame(results)

    # Multiple hypothesis correction per model comparison across genes
    
    # Collect all LRT p-values across all comparisons (columns that contain 'lrt_pval')
    lrt_pval_columns = [col for col in results_df.columns if 'lrt_pval' in col]

    # Flatten all LRT p-values into a single list
    all_pvals = results_df[lrt_pval_columns].values.flatten()

    # Remove NaN values from the list
    all_pvals = all_pvals[~pd.isna(all_pvals)]

    # Apply Benjamini-Hochberg correction across all comparisons and all p-values
    _, pvals_corrected, _, _ = multipletests(all_pvals, method='fdr_bh')

    # Reshape the corrected p-values back into the same shape as the original DataFrame
    corrected_pvals = pd.Series(pvals_corrected).values.reshape(results_df[lrt_pval_columns].shape)

    # Assign the corrected p-values back to the corresponding columns in the DataFrame
    for i, col in enumerate(lrt_pval_columns):
        results_df[f'{col}_corrected'] = corrected_pvals[:, i]

    results_df.to_csv(f'241003_GLM_Updated_results_{cell_state}.csv')
# Parallelize the cell state processing using Pool
if __name__ == '__main__':
    with Pool(processes=len(cell_states)) as pool:
        pool.map(process_cell_state, cell_states)