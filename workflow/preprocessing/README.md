# Preprocessing


This module provides a flexible preprocessing pipeline for single-cell datasets, with rules for:

1. Normalization
2. Gene filtering
3. Highly variable gene selection
4. PCA
5. Nearest neighbor graph construction
6. UMAP dimensionality reduction
7. Optional extra HVG selection and custom assembly
8. UMAP/PCA plotting (colors, gene sets, centroids)

Rules are defined in `rules/rules.smk` and orchestrated in `rules/assemble.smk`, which assembles outputs from each step into a single AnnData object. The pipeline supports parameterization via config and dynamic resource allocation (threads, memory, partition, GPU) for each rule.


## General input and configuration

Each rule accepts either a `.h5ad` or `.zarr` file as input and always writes output as a `.zarr` file for efficient storage. Input AnnData must include:

- Raw counts in `.X` (default) or `.layers` (as specified by `raw_counts` in config)
- Optional batch/lineage columns in `.obs` for HVG selection

AnnData files from the `load_data` workflow are compatible by default. Each step saves only the changed slots, and the final `assemble` rule collects specified outputs into a single AnnData object. The output is saved as `config["output_dir"] + "/dataset_name/preprocessed.zarr"`.

The pipeline is highly configurable: parameters for each step, resource allocation (threads, memory, partition, GPU), and conda environments can be set via config and rule directives. See `workflow/preprocessing/rules/assemble.smk` for details.

```yaml
DATASETS:
  dataset_name:
    input:
      preprocessing: adata.h5ad
    preprocessing:
      raw_counts: X  # default
      assemble:
        - counts  # raw counts that are provided as input in .layers['counts']
        - normalize # normalised counts
        - highly_variable_genes
        - pca
        - neighbors
        - umap
```

## Resource configuration

Resource settings such as Dask usage, number of threads, and memory/partition/GPU allocation can be set globally in your config file for more scalable implementations. These defaults are applied to all preprocessing steps unless overridden for a specific rule.

- `dask`: Whether to use Dask for parallelization (default: True)
- `n_threads`: Number of threads for parallel operations (default: 10)
- `resources`: A string indicating 

For medium to larger datasets (500k+ cells), using the settings below is recommended.
For smaller datasets, you can disable dask and GPU support to avoid computational overhead.

**Global config example:**
```yaml
preprocessing:
  dask: true
  n_threads: 10
  resources: gpu
```

You can override these whether `dask` for each step of the pipeline individually. If not specified, the pipeline uses the global defaults.

**Step-specific override example (only for `dask`)**
```yaml
preprocessing:
  dask: true
  n_threads: 10
  resources: gpu
  highly_variable_genes:
    dask: false # overrides global default for highly_variable_genes step
```


## Preprocessing steps and rule structure

Each step is implemented as a Snakemake rule, with parameters, resources, and conda environments set via config and rule lambdas. If parameters are not defined, sensible defaults are used.


### Normalize

Transforms the data using `scanpy.pp.normalize_total` and log-transforms with `scanpy.pp.log1p`.

**Output**

- `.X`: normalized and log1p-transformed counts (sparse)
- `.uns["preprocessing"]`: metadata on normalization

**Parameters**

- `raw_counts`: Key for input matrix
- `n_threads`: Number of threads (default: 10)
- `dask`: Use Dask for parallelization (default: True)

**Resources/Conda**
- Partition, GPU, memory, and conda environment are set dynamically via config and rule lambdas.

**Config example**

```yaml
DATASETS:
  dataset_name:
    input:
      preprocessing: 'adata.h5ad'
    preprocessing:
      raw_counts: X
      assemble:
        - normalize
```


### Filter Genes

Filters genes using `scanpy.pp.filter_genes` before HVG selection.

**Output**
- `.zarr` file with filtered genes

**Parameters**
- `dask`: Use Dask (default: True)
- `n_threads`: Number of threads (default: 10)

**Resources/Conda**
- Partition, GPU, memory, and conda environment set via config/rule lambdas.

### Highly Variable Genes Selection

Highly variable genes are calculated after using `scanpy.pp.filter_genes(min_cells=1)` and the `.var` in the unfiltered object is updated.
Note, that the highly variable gene selection will be run on the normalisation output.

**Output**

- `.var["highly_variable"]` boolean column with highly variable gene status.
- `.uns["preprocessing"]["highly_variable_genes"]` containing all the arguments passed to the `scanpy.pp.highly_variable_genes` function (`batch_key` included)

**Parameters**

- `batch`: batch column to account when getting the HVGs.
- `lineage`: if you want to calculate lineage specific genes alone
  or combined with batch.
- `highly_variable_genes`: any additional arguments to pass to the `scanpy.pp.highly_variable_genes` function.

**Config example**

```yaml
DATASETS:
  dataset_name:
    input:
      ...
    preprocessing:
      ...
      highly_variable_genes:
        n_top_genes: 2000
        subset: true  # will subset the feature space of the object, resulting object will be smaller
      assemble:  # only include highly_variable_genes output here, normalized count matrix won't be saved
        - highly_variable_genes
```


### PCA

If `scale=True` scaling of `.X` will be performed, otherwise the counts in `.X` are used directly.
Subsequently, the PCA is calculated using the highly variable genes.
Note, that this step requires the normalized and log-transformed counts in `.X` as well as the highly variable genes information in `.var`.
As part of the pipeline, the `normalization` and `highly_variable_genes` rules will be used as input.

**Output**

- `.obsm["X_pca"]` PCA embedding
- `.uns["preprocessing"]["pca"]` containing all the arguments passed to the `scanpy.pp.pca` function (`use_highly_variable=True` included)
- `.uns["preprocessing"]["scale"]`: whether the matrix was scaled or not before PCA

**Parameters**

- `scale`: wether or not to scale counts before calculating
the PCA embedding.
- `pca`: any additional arguments to pass to the `scanpy.pp.pca` function.

**Config example**

```yaml
DATASETS:
  dataset_name:
    input:
      ...
    preprocessing:
      ...
      pca:
        n_comps: 50
      assemble:
        - pca
```



### K-nearest neighbor graph

It will attempt to use the RAPIDS[^1] implementation but will default, if it fails, to the UMAP implementation
[arXiv:1802.03426v3](https://arxiv.org/abs/1802.03426v3).

By default, the rule will compute the neighbors based on the PCA distances in `.obsm["X_pca"]`, if available, or on `.X` directly otherwise.
In order to use a different representation, `use_rep` must be specified under the additional parameters.

**Output**

- `.obsp["neighbors"]["distances"]` distance matrix
- `.obsp["neighbors"]["connectivities"]` adjacency matrix
- `.uns['preprocessing']['scale']`: whether the matrix was scaled or not before PCA

**Parameters**

- `neighbors`: any additional arguments for the `scanpy.pp.neighbors` function.

**Config example**

```yaml
DATASETS:
  dataset_name:
    input:
      ...
    preprocessing:
      ...
      neighbors:
        use_rep: X_pca
      assemble:
        - neighbors
```


### UMAP

UMAP dimensionality reduction is calculated from the PCA output. It can also use the RAPIDS[^1] implementation.

**Output**

- `.obsm["X_umap"]` UMAP embedding; it becoms .obsm["X_umap_{key}"] if multiple `neighbors_key` are given in the config file.

**Parameters**

- `umap`: any additional arguments for the `scanpy.tl.umap` function.

**Config example**

```yaml
DATASETS:
  dataset_name:
    input:
      ...
    preprocessing:
      ...
      umap:
        neighbors_key: neighbors
      assemble:
        - umap
```

## Assembled output and custom assembly


The `assemble` rule collects user-specified preprocessing outputs and saves them in a single AnnData object (zarr file). You can customize which slots to include via the `assemble` parameter in your config. The rule supports dynamic resource allocation and conda environments.

**Output**

- `.layers["counts"]`: raw counts (from input) if `counts` is present under `assemble`
- `.X` and `.layers["normcounts"]`: normalized counts if `normalize` is present under `assemble`
- `.var[["highly_variable", "means", "dispersions", "dispersions_norm", "highly_variable_nbatches", "highly_variable_intersection"]]`: highly variable gene information if `highly_variable_genes` is present under `assemble`
- `.obsm["X_pca"]`: PCA representation if `pca` is present under `assemble`
- `.uns["neighbors"]`, `.obsp["distances"]`, `.obsp["connectivities"]`: kNN graph if `neighbors` is present under `assemble`
- `.obsm["X_umap"]`: UMAP representation if `umap` is present under `assemble`
- `.uns["preprocessing"]` containing the metadata on how the different preprocessing steps were run. See each preprocessing step for more detailed descriptions of the preprocessing metadata.


**Parameters**
- `assemble`: list of preprocessing outputs (including raw counts) to assemble in the final output


**Config example**

```yaml
DATASETS:
  dataset_name:
    input:
      ...
    preprocessing:
      assemble:
        - counts  # raw counts that are provided as input in .layers['counts']
        - normalize # normalised counts
        - highly_variable_genes
        - pca
        - neighbors
        - umap
```


### UMAP/PCA Plots

Generate publication-ready scatter plots for PCA and UMAP embeddings. You can plot categorical/numeric columns from `.obs` and visualize gene expression panels, including predefined gene sets.
Plotting will happen by default when `preprocessing_all` is called.

Plots are produced by rules in [workflow/preprocessing/rules/plots.smk](workflow/preprocessing/rules/plots.smk) and saved under your configured image directory (per dataset and step) in `pca/` and `umap/` folders.

**Output**

- `images/.../pca/*.png`: PCA scatter plots for requested colors
- `images/.../umap/*.png`: UMAP scatter plots for requested colors and gene panels

**Parameters**
The plot configuration is independent of `assembly`.
Even if no `pca` or `umap` are defined, the PCA and UMAP workflows will be triggered in order to create the plots.

- `colors`: list of `.obs` columns and/or genes to plot. Entries not found in `.obs` are treated as genes or gene name patterns and matched against `.var_names` to create expression panels.
- `plot_centroids`: list of categorical `.obs` columns for which to overlay category numbers at centroid positions (enabled for up to 102 categories). Legends are annotated to map numbers back to labels.
- `plot_gene_chunk_size`: number of genes per panel when plotting many genes together (default: 12).

**Config example**

```yaml
DATASETS:
  dataset_name:
    input:
      ...
    preprocessing:
      # Plot settings
      colors:
        - batch
        - dataset
        - cell_type
        - study
        - development_stage
        - CCR7
        - PTPRC
      plot_centroids:
        - cell_type
      plot_gene_chunk_size: 12
```

Notes on plotting genes:
- Genes/patterns listed in `colors` are matched to `.var_names` (regexes will be evaluated).
- When many genes are requested, they are grouped into panels of `plot_gene_chunk_size`, with `ncols` controlling the panel layout.


---

**Notes:**
- RAPIDS implementation for neighbors/UMAP is used if available and configured.
- All resource and environment settings are controlled via config and rule lambdas in `assemble.smk`.
