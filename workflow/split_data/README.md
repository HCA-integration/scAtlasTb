# Split Data

This module splits single-cell data (AnnData objects) into multiple files based on categorical values in the `.obs` metadata. Each unique value in the specified column becomes a separate output file.

## Functionality

The script:
- Reads an AnnData file (`.h5ad` or `.zarr` format)
- Splits cells based on a categorical column in `.obs`
- Outputs each subset as a separate `.zarr` file
- Supports both memory-based copying and efficient linked subsets
- Adds wildcard annotations to track split metadata

## Input Parameters

- `input`: Path to input AnnData file
- `split_key`: Column name in `.obs` to split by (passed as wildcard)
- `values`: List of specific values to extract (uses sanitized filenames)
- `dask`: Whether to keep arrays as dask arrays before writing copy (default: False)
- `write_copy`: Whether to write full copies vs linked subsets (default: False, auto-enabled for .h5ad inputs)
- `slots`: Optional mapping of slots to read/write

## Example Config

```yaml
output_dir: test/out
images: test/images

DATASETS:
  test:
    input:
      split_data:
        pbmc: test/input/pbmc68k.h5ad
    split_data:
      key: bulk_labels
      values:
        - CD4+_CD45RA+_CD25-_Naive_T
        - Dendritic
        - CD14+_Monocyte
        - CD19+_B
```

## Output Structure

The script creates files with sanitized names (spaces and slashes replaced with underscores):

```shell
test/out/split_data
├── dataset~test
│   └── file_id~pbmc
│       └── key~bulk_labels
│           ├── value~CD14+_Monocyte.zarr -> ../../../splits/dataset~test/file_id~pbmc/key~bulk_labels/value~CD14+_Monocyte.zarr
│           ├── value~CD19+_B.zarr -> ../../../splits/dataset~test/file_id~pbmc/key~bulk_labels/value~CD19+_B.zarr
│           ├── value~CD4+_CD45RA+_CD25-_Naive_T.zarr -> ../../../splits/dataset~test/file_id~pbmc/key~bulk_labels/value~CD4+_CD45RA+_CD25-_Naive_T.zarr
│           └── value~Dendritic.zarr -> ../../../splits/dataset~test/file_id~pbmc/key~bulk_labels/value~Dendritic.zarr
├── input_files.tsv
└── splits
    └── dataset~test
        └── file_id~pbmc
            └── key~bulk_labels
```

Each output file contains:
- Subset of cells matching the split value
- All original `.var` data
- Added wildcard annotations in `.uns` tracking split metadata
- Either full data copies or efficient links to original file (depending on `write_copy` parameter)

## Performance Notes

- For `.zarr` inputs with `write_copy=False`: Creates efficient linked subsets
- For `.h5ad` inputs: Always creates full copies due to format limitations
- Uses Dask for memory-efficient processing of large datasets when reading file (relevant when `write_copy=True` or `.h5ad` as input file)