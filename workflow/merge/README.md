# Merging

This workflow merges multiple AnnData files into a single consolidated dataset.
It supports various merge strategies and can handle large datasets using Dask and backed mode for memory efficiency.

## Configuration

```yaml
DATASETS:
  Lee2020:
    input:
      merge:
        file_1: test/input/load_data/harmonize_metadata/Lee2020.zarr
        file_1.1: test/input/load_data/harmonize_metadata/Lee2020.zarr
        file_1.2: test/input/load_data/harmonize_metadata/Lee2020.zarr
        file_1.3: test/input/load_data/harmonize_metadata/Lee2020.zarr
        file_1.4: test/input/load_data/harmonize_metadata/Lee2020.zarr
        file_2: test/input/load_data/harmonize_metadata/SchulteSchrepping2020.zarr
        file_2.1: test/input/load_data/harmonize_metadata/SchulteSchrepping2020.zarr
        file_2.2: test/input/load_data/harmonize_metadata/SchulteSchrepping2020.zarr
        file_2.3: test/input/load_data/harmonize_metadata/SchulteSchrepping2020.zarr
        file_2.4: test/input/load_data/harmonize_metadata/SchulteSchrepping2020.zarr
        file_2.5: test/input/load_data/harmonize_metadata/SchulteSchrepping2020.zarr
    merge:
      merge_strategy: inner
      keep_all_columns: true
      allow_duplicate_obs: true
      allow_duplicate_var: false
      threads: 5
      stride: 500_000
      dask: true
      backed: true
      slots:
        X: X
        obs: obs
        var: var
        layers: layers
  test:
    input:
      merge:
        file_1: test/input/pbmc68k.h5ad
        file_2: test/input/pbmc68k.h5ad
    merge:
      merge_strategy: outer
      allow_duplicate_obs: true
```

## Parameters

- **`merge_strategy`**: Strategy for handling overlapping genes/variables
  - `"inner"` (default): Keep only genes present in all datasets
  - `"outer"`: Keep all genes from all datasets, filling missing values with zeros

- **`keep_all_columns`**: Boolean flag for handling observation metadata
  - `true`: Merge all observation columns from all datasets
  - `false` (default): Only keep columns present in all datasets

- **`allow_duplicate_obs`**: Allow duplicate observation (cell) names in the merged object
  - `true`: Duplicates are allowed
  - `false` (default): Duplicates are not allowed and will raise an error

- **`allow_duplicate_var`**: Allow duplicate variable (gene) names in the merged object
  - `true`: Duplicates are allowed
  - `false` (default): Duplicates are not allowed and will raise an error

- **`threads`**: Number of threads for Dask processing
  - Integer value (default: system dependent)
  - Used when `dask: true`

- **`stride`**: Chunk size for processing large datasets
  - Integer value (default: 500,000)
  - Controls memory usage and processing efficiency

- **`dask`**: Boolean flag to enable Dask for distributed processing
  - `true`: Use Dask arrays for memory-efficient processing of large datasets
  - `false` (default): Use standard in-memory processing

- **`backed`**: Boolean flag to enable backed mode for AnnData
  - `true`: Use backed mode with AnnCollection for memory-efficient merging
  - `false` (default): Load all data into memory

- **`slots`**: Dictionary specifying which data components to read from zarr files
  - Keys: slot names (`X`, `obs`, `var`, `layers`, etc.)
  - Values: corresponding zarr group names
  - Only relevant for zarr input files

## Processing Modes

The workflow supports three different processing modes:

### 1. Dask Mode (`dask: true`)
- Loads all files using Dask arrays
- Memory-efficient for very large datasets
- Supports parallel processing
- Best for datasets that don't fit in memory

### 2. Backed Mode (`backed: true, dask: false`)
- Uses AnnCollection for efficient merging
- Files remain on disk during processing
- Good balance between memory efficiency and performance
- Suitable for moderately large datasets

### 3. Standard Mode (`dask: false, backed: false`)
- Loads files sequentially into memory
- Merges using scanpy.concat
- Fastest for datasets that fit in memory
- Default mode for smaller datasets

## Behavior

- **Empty Dataset Handling**: Automatically skips empty datasets during merging
- **Single File**: If only one file is provided, creates a symbolic link instead of merging
- **Gene Annotation Merging**: Combines gene metadata from all datasets
- **Index Generation**: Creates new cell identifiers in format `{dataset}-{index}`
- **Metadata Preservation**: Stores original cell names in `obs_names_before_{dataset}` column
- **Dataset Tracking**: Adds dataset identifier to `adata.obs['file_id']` and `adata.uns['dataset']`
- **Duplicate Handling**: Duplicate cell or gene names are checked and controlled by `allow_duplicate_obs` and `allow_duplicate_var` parameters

## Input/Output

- **Input**: Multiple AnnData files in `.zarr` or `.h5ad` format
- **Output**: Single merged AnnData file in `.zarr` format
- **File Naming**: Input files are identified by their keys in the `input.merge` dictionary

## Memory Optimization

The workflow includes several memory optimization features:
- Automatic garbage collection after processing each file
- Slot removal to free memory
- Sparse matrix preservation
- Chunked processing for large datasets
- Progress tracking with tqdm integration
