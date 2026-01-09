# Cell Type Prediction

This module enables automated cell type annotation of single-cell RNA-seq datasets using pre-trained models. This allows for:

- Consistent cell type annotation across studies
- Probabilistic cell type assignments with confidence scores
- Majority voting and over-clustering analysis for robust predictions

The module currently supports:

- **CellTypist**: Automated cell type annotation using pre-trained models

## Environments

The following environments are needed for cell type prediction:

- [`celltypist`](https://github.com/HCA-integration/scAtlasTb/blob/main/envs/celltypist.yaml)

## Configuration

```yaml
DATASETS:
  test:
    input:
      celltype_prediction:
        preprocessed: test/input/preprocessing/dataset~all/file_id~pbmc/preprocessed.zarr
    celltype_prediction:
      reference_label: bulk_labels
      counts: layers/counts
      is_normalized: false
      celltypist:
        params:
          majority_voting: true
          over_clustering: bulk_labels
        models:
          - Healthy_COVID19_PBMC
          - Immune_All_Low
```

### Input

The input AnnData object should contain the single-cell RNA-seq data to be annotated.
Reference cell type labels can optionally be provided for evaluation and visualization purposes.

### Configuration Parameters

- **`counts`** (default: `'X'`): Which data layer to use from the AnnData object

- **`is_normalized`** (default: `true`): Boolean flag indicating whether the input data is already normalized
  - CellTypist expects log-normalized data, so this parameter controls preprocessing

- **`reference_label`** (optional): Column name in `.obs` containing reference cell type labels for comparison and visualization

#### CellTypist Parameters (`celltypist`)
Configuration for CellTypist cell type prediction:

- **`models`**: List of pre-trained model names to use for prediction (required)
  - Available models include `Healthy_COVID19_PBMC`, `COVID19_HumanChallenge_Blood`, `Immune_All_Low`, etc.
  - Multiple models can be applied to the same dataset

- **`params`**: Model parameters (optional)
  - **`majority_voting`**: Enable majority voting across over-clustering results (default: `false`)
  - **`over_clustering`**: Column name in `.obs` for over-clustering analysis (optional)

> **Note:** CellTypist models are trained on specific tissue types and cell populations. Choose models appropriate for your data type (e.g., PBMC, immune cells, etc.).

## Output

### CellTypist
The cell type prediction workflow produces the following outputs:

* `<out_dir>/celltype_prediction/dataset~<dataset>/file_id~<file_id>.zarr`: Annotated AnnData object containing:
  - **Direct predictions** (`obs['celltypist_<model>:predicted_labels']`): Primary cell type predictions
  - **Majority voting results** (`obs['celltypist_<model>:majority_voting']`): Consensus predictions (if enabled)
  - **Over-clustering results** (`obs['celltypist_<model>:over_clustering']`): Fine-grained clustering results (if specified)
  - **Confidence scores** (`obs['celltypist_<model>:conf_score']`): Prediction confidence values

* `<out_dir>/images/`: Visualization plots comparing predictions with reference labels (if provided)
