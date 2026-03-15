# Preprocessing

This module runs a configurable single-cell preprocessing workflow and assembles selected results into one output zarr file.

Implemented steps:

1. Normalize (`normalize.py`)
2. Mark non-zero genes for HVG-safe filtering (`filter_genes.py`)
3. Highly variable genes (`highly_variable_genes.py`)
4. Optional extra HVGs (`extra_hvgs.py`)
5. PCA (`pca.py`)
6. Neighbors graph (`neighbors.py`)
7. UMAP (`umap.py`)
8. PCA/UMAP plots (`plot.py`)
9. Assembly (`assemble.py`)

Rules are declared in `rules/rules.smk`, parameterized in `rules/assemble.smk`, and plotting rules are in `rules/plots.smk`.


## Inputs, outputs, and execution

- The module reads dataset inputs from `DATASETS.<dataset>.input.preprocessing`.
- Inputs may be `.h5ad` or `.zarr`.
- Rule outputs are written as `.zarr`.
- The default target (`rule all` in `Snakefile`) runs:
  - assembly output zarr files
  - PCA and UMAP plots

The assembled output path follows the module parameter space and is written under the configured output directory.


## Global configuration

Common knobs (global defaults, overridable per step):

- `dask`
- `n_threads`
- `resources`
- step-specific argument dictionaries (`normalize`, `filter`, `highly_variable_genes`, `extra_hvgs`, `pca`, `neighbors`, `umap`)

Example:

```yaml
preprocessing:
  dask: true
  n_threads: 10
  resources: gpu
```


## Step details

### Normalize

Script: `scripts/normalize.py`

Behavior:

- Reads raw counts from `raw_counts` (default `X`)
- Ensures sparse representation
- Runs `scanpy.pp.normalize_total` and `scanpy.pp.log1p`
- Stores normalized matrix in both `.X` and `.layers["normcounts"]`
- Preserves raw counts in `.layers["counts"]` and `.raw`
- Writes metadata under `.uns["preprocessing"]` and `.uns["log1p"]`

Important params:

- `raw_counts`
- `gene_id_column`
- `normalize` (args passed to `scanpy.pp.normalize_total`)
- `dask`

Example:

```yaml
DATASETS:
  dataset_name:
    preprocessing:
      raw_counts: X
      gene_id_column: gene_id
      normalize:
        target_sum: 1e4
```


### Filter genes

Script: `scripts/filter_genes.py`

Behavior:

- Computes `.var["nonzero_genes"]` using `_filter_genes`.
- This is a marker step used by HVG/extra-HVG steps; it does not shrink the final feature space by itself.

Important params:

- `filter` (default contains `min_cells: 1`)
- `dask`

Example:

```yaml
DATASETS:
  dataset_name:
    preprocessing:
      filter:
        min_cells: 3
```


### Highly variable genes

Script: `scripts/highly_variable_genes.py`

Behavior:

- Runs HVG selection on filtered cells/genes.
- Maps HVG results back to the full `.var`.
- Always provides `highly_variable`; additionally writes a parameterized variant column
  `highly_variable-...` when args are provided.
- Stores run args in `.uns["preprocessing"]["highly_variable_genes"]`.

Notes:

- `subset` is explicitly removed from args in the script.
- If args are `False`, all genes are marked as highly variable.

Example:

```yaml
DATASETS:
  dataset_name:
    preprocessing:
      highly_variable_genes:
        n_top_genes: 2000
        flavor: seurat_v3
        batch_key: sample
```


### Extra HVGs

Script: `scripts/extra_hvgs.py`

Behavior:

- Computes an additional HVG mask in `extra_hvgs` (or `extra_hvgs-...` when `overwrite_args` is used).
- Supports:
  - union of per-group HVGs via `union_over`
  - adding genes via `extra_genes`
  - removing genes via `remove_genes`
- Stores metadata in `.uns["preprocessing"][<extra_hvg_column>]`.

Example:

```yaml
DATASETS:
  dataset_name:
    preprocessing:
      extra_hvgs:
        union_over: [lineage]
        extra_genes: [CCR7, PTPRC]
        remove_genes: [MALAT1]
        min_cells: 200
        overwrite_args:
          n_top_genes: 3000
          flavor: seurat_v3
```


### PCA

Script: `scripts/pca.py`

Behavior:

- Subsets to HVGs using `mask_var` (default `highly_variable`).
- Optionally scales before PCA (`scale`).
- Writes `.obsm["X_pca"]`, `.uns["pca"]`, and `.varm` loadings.
- Stores preprocessing metadata in `.uns["preprocessing"]["pca"]` and `...["scaled"]`.

Example:

```yaml
DATASETS:
  dataset_name:
    preprocessing:
      scale: true
      pca:
        n_comps: 50
        svd_solver: covariance_eigh
```


### Neighbors

Script: `scripts/neighbors.py`

Behavior:

- Computes or reuses neighbor graph.
- If `neighbors` params are `False`, reuses existing graph from input.
- Defaults to `use_rep="X_pca"` when available, otherwise `X`.
- Writes `.obsp["distances"]`, `.obsp["connectivities"]`, `.uns["neighbors"]`.

Example:

```yaml
DATASETS:
  dataset_name:
    preprocessing:
      neighbors:
        n_neighbors: 15
        metric: cosine
        use_rep: X_pca
```


### UMAP

Script: `scripts/umap.py`

Behavior:

- Computes UMAP from a selected neighbors graph (`neighbors_key`, default `neighbors`).
- If required representation is missing in input, it is loaded from the `rep` input.
- Writes `.obsm["X_umap"]` and updated `.uns`.

Example:

```yaml
DATASETS:
  dataset_name:
    preprocessing:
      umap:
        min_dist: 0.3
        spread: 1.0
```


### Plots

Script: `scripts/plot.py`

Rules:

- `preprocessing_plot_pca`: basis `X_pca`
- `preprocessing_plot_umap`: basis `X_umap`

Behavior:

- Plots obs columns from `colors`.
- Treats non-obs `colors` entries as genes/patterns and creates expression panels.
- Supports centroid overlays for categorical columns via `plot_centroids`.
- Uses `plot_gene_chunk_size` to chunk gene panels.

Example:

```yaml
DATASETS:
  dataset_name:
    preprocessing:
      colors: [batch, cell_type, CCR7, PTPRC]
      plot_centroids: [cell_type]
      plot_gene_chunk_size: 12
```


## Assembly

Script: `scripts/assemble.py`

The `assemble` list controls which step outputs are linked into the final object.

Supported `assemble` entries:

- `normalize`
- `highly_variable_genes`
- `extra_hvgs`
- `pca`
- `neighbors`
- `umap`

Notes:

- `counts` is not a standalone assembly key in the current implementation.
- For HVG-like outputs, assembly links parameterized columns (for example
  `highly_variable-...`, `extra_hvgs-...`) and also fills default slots
  (`highly_variable`, `extra_hvgs`) for the default file of each type.
- Wildcards used to build the assembled output are stored under `uns/wildcards`.

Example:

```yaml
DATASETS:
  dataset_name:
    preprocessing:
      assemble:
        - normalize
        - highly_variable_genes
        - extra_hvgs
        - pca
        - neighbors
        - umap
```


## Minimal config example

```yaml
DATASETS:
  dataset_name:
    input:
      preprocessing: adata.h5ad
    preprocessing:
      raw_counts: X
      dask: true
      n_threads: 10
      highly_variable_genes:
        n_top_genes: 2000
        batch_key: batch
      extra_hvgs:
        union_over: [lineage]
        extra_genes: [CCR7, PTPRC]
      pca:
        n_comps: 50
      neighbors:
        n_neighbors: 15
      umap:
        min_dist: 0.5
      colors: [batch, cell_type, CCR7]
      plot_centroids: [cell_type]
      plot_gene_chunk_size: 12
      assemble:
        - normalize
        - highly_variable_genes
        - extra_hvgs
        - pca
        - neighbors
        - umap
```


## Notes

- GPU-enabled rules use RAPIDS via the configured GPU environment when available.
- Empty datasets are explicitly handled in scripts and written as valid empty zarr outputs.
