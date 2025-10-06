# Reference mapping

This module enables projection of new single-cell RNA-seq datasets onto pre-trained variational autoencoder (VAE) models or foundation models. This allows for:

- Transfer learning from large reference datasets to smaller query datasets
- Integration of new data with existing atlases
- Consistent cell type annotation across studies
- Leveraging pre-trained embeddings for downstream analysis

## Environments

The following environments are needed for reference mapping:

- [`scarches`](https://github.com/HCA-integration/scAtlasTb/blob/main/envs/scarches.yaml)

## Supported Models

The module currently supports:
- **scvi-tools models**: Various variational inference models (scVI, scANVI, etc.) mapping with scArches
- TODO: **scArches models**: VAE-based models (trVAE, scPoli, Expimap, etc.)
- TODO: **Foundation models**: Large-scale pre-trained models

## Configuration

```yaml
DATASETS:
  test:
    input:
      reference_mapping:
        file_1: test/input/query_data.h5ad
    reference_mapping:
      layer: X  # or layer name like 'counts'
      model_params:
        batch_key: sample_id
        labels_key: cell_type
        categorical_covariate: [donor, condition]
        continuous_covariate: [age]
      train_kwargs:
        max_epochs: 10
        early_stopping: true
        check_val_every_n_epoch: 1
```

### Input

The input AnnData object is the query dataset that should be mapped to the reference model.
The reference model should be defined as a Pytorch model ('model.pt' file) under `scarches > model`.

### Configuration Parameters

- `layer`: Which data layer to use from the query AnnData ('X' for .X, or layer name)
- `model_params`: Model-specific parameters for alignment
  - `batch_key`: Column name in `.obs` for batch information
  - `labels_key`: Column name in `.obs` for cell type labels
  - `categorical_covariate`: List of categorical covariate column names in `.obs`
  - `continuous_covariate`: List of continuous covariate column names in `.obs`
- `train_kwargs`: Training parameters for the mapping process
  - `max_epochs`: Maximum number of training epochs
  - `early_stopping`: Whether to use early stopping
  - `check_val_every_n_epoch`: Validation frequency

> Note: The query data must have genes that overlap with the reference model's training data. Assuming that gene naming is consistent between the pre-trained model and the input AnnData object (`var_key`), the data is subsetting to the overlapping genes before training.

## Output

The reference mapping workflow produces the following outputs:

* `<out_dir>/reference_mapping/dataset~<dataset>/file_id~<file_id>.zarr`: Mapped AnnData object containing:
  - **Latent embedding** (`obsm['X_emb']`): Low-dimensional representations in the reference space

* `<out_dir>/reference_mapping/model/dataset~<dataset>/file_id~<file_id>/`: Updated model directory with query data integrated, ready for further reference mapping tasks
