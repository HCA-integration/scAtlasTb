# The `batch_analysis` module

This module provides an exploratory framework for understanding the technical effects of batch variables at different hierarchical levels in single-cell datasets. It quantifies the strength of technical covariates by measuring the linear variance they explain using principal component regression (PCR). This systematic evaluation helps to assess the impact of each covariate, making it easier to decide which variable should be treated as the batch for subsequent batch correction or data integration steps.

## Features
- Automated batch effect analysis and permutation-based significance testing
- Optional preprocessing pipeline: normalization, gene filtering, PCA
- Parallelized computation of batch PCR for scalability
- Flexible configuration for covariates, permutations, and sample keys
- Output includes per-covariate PCR scores, permutation z-scores, and summary plots

## Environments

The following conda environments are used by different steps:
- `scanpy` (for preprocessing, PCA, plotting)
- `scib` (for PC regression and batch PCR)

## Configuration

Configure the module under your dataset key using the `batch_analysis` section. Common keys:

- `sample`: column(s) in `.obs` to use as the sample key (can be comma-separated for composite keys). It should represent the smallest common grouping of a technical effect and the covariate of interest. An error will be thrown if this is not the case.
- `covariates`: list of covariate columns in `.obs` to test for batch effects
- `permute_covariates`: (optional) list of covariates to permute for computing a z-score. If not determined, all `covariates` will be used for perturbations. The covariates of interest will be permuted per 
- `n_permutations`: number of permutations for each covariate
- Step-specific overrides (e.g., `normalize`, `highly_variable_genes`, `pca`) for preprocessing

> **Note:** Preprocessing steps such as normalization, and PCA are optional. Each step will only be executed if its corresponding key (`normalize`, `highly_variable_genes`, `pca`) is defined in your configuration. If a key is omitted, that step will be skipped for the dataset.

### Example configuration

```yaml
DATASETS:
  precomputed_pca:
    input:
      batch_analysis: test/input/blood_pca.zarr
    batch_analysis:
      covariates:
        - sample
        - donor
        - assay
        - sex
        - disease
        - self_reported_ethnicity
      permute_covariates:
        - assay
        - sex
        - disease
        - self_reported_ethnicity
      n_permutations: 1000
      sample: sample,donor

  recompute_pca:
    input:
      batch_analysis: test/input/pbmc68k.h5ad
    batch_analysis:
      sample: batch, bulk_labels
      covariates:
        - bulk_labels
        - batch
        - is_cd14_mono
      normalize:
      highly_variable_genes:
      pca:
```

The example configuration above demonstrates how to set up the `batch_analysis` module for two datasets:

- **`precomputed_pca`**: Uses a dataset with precomputed PCA in the input AnnData object. It specifies multiple covariates (e.g., `sample`, `donor`, `assay`, `sex`, `disease`, `self_reported_ethnicity`) to test for batch effects, and defines which covariates to permute for significance testing. The number of permutations is set to 1000, and a composite sample key (`sample,donor`) is used.

- **`recompute_pca`**: Uses a dataset without the necessary PCA information and configures the workflow to perform normalization, highly variable gene selection, and PCA as preprocessing steps. It sets `batch` and `bulk_labels` as the composite sample key, and tests covariates such as `bulk_labels`, `batch`, and `is_cd14_mono` for batch effects.


## Workflow steps

The batch_analysis workflow consists of the following steps:

1. **Preprocessing** (normalize, filter genes, HVG selection, PCA):
   - Uses rules from the preprocessing module, with dataset-specific overrides.
2. **Covariate setup** (`determine_covariates`):
   - Determines which covariates to test and sets up permutation schemes.
3. **Batch PCR** (`batch_pcr`):
   - Runs principal component regression for each covariate and permutation, computes z-scores.
4. **Collect**:
   - Aggregates per-covariate results into a single table.
5. **Plot**:
   - Generates barplots and violin plots summarizing PCR and permutation results.

## Output

- Per-covariate PCR results: `<output_dir>/<dataset>/batch_pcr/{covariate}.tsv`
- Aggregated results: `<output_dir>/<dataset>/batch_pcr.tsv`
- Plots: `<images>/<dataset>/batch_pcr_bar.png`, `<images>/<dataset>/batch_pcr_violin.png`

Each result file contains columns such as:
- `covariate`: tested covariate
- `pcr`: PCR score
- `permuted`: whether the score is from a permutation
- `n_covariates`: number of unique values in the covariate
- `z_score`: z-score of observed PCR vs. permutations

## Testing

Activate the snakemake environment and run the test workflow:

```bash
conda activate snakemake
bash test/run_test.sh -n
```

This test requires precomputed objects from the `preprocessing` tests.
See the `preprocessing` module for more details.