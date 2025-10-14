# Filtering

This workflow filters the input data based on a set of conditions.
The conditions can be defined in the config under the `filter` module.

## Configuration

```yaml
DATASETS:
  Lee2020:
    input:
      filter: test/input/load_data/harmonize_metadata/Lee2020.zarr
    filter:
      subset: false ## whether to subset the file according to the filters
      remove_by_column:
        sample: # column name
          # entries from column to exclude (will be treated as String)
          - Schulte-Schrepping_C2P01H_d0
          - Schulte-Schrepping_C2P05F_d0
          - Schulte-Schrepping_C2P07H_d0
          - Schulte-Schrepping_C2P10H_d0
          - Schulte-Schrepping_C2P13F_d0
          - Schulte-Schrepping_C2P15H_d0
          - Schulte-Schrepping_C2P16H_d0
          - Schulte-Schrepping_C2P19H_d0
        donor: # column name
          # entries from column to exclude (will be treated as String)
          - C19-CB-0008 
        disease:
          - influenza
      remove_by_query:
        # pandas query strings for complex filtering conditions
        - 'random < 3'
  test:
    input:
      filter: test/input/pbmc68k.h5ad
    filter:
      subset: true
      remove_by_column:
        phase:
          - G1
        is_cd14_mono:
          - true
```

## Parameters

- **`subset`**: Boolean flag that determines whether to physically subset the file according to the filters
  - `true` (default): Creates a new subset file with filtered data
  - `false`: Only applies filters without modifying the original file structure

- **`remove_by_column`**: Dictionary defining filtering conditions
  - Keys: Column names in the AnnData observation metadata (`adata.obs`)
  - Values: List of entries to exclude from that column
  - All entries are treated as strings for comparison

- **`remove_by_query`**: List of pandas query strings for complex filtering conditions
  - Each query string follows pandas DataFrame query syntax
  - Cells matching ANY query will be removed
  - Allows for complex numerical and logical operations

## Behavior

- The different conditions in `remove_by_column` are combined with an **AND** operation
- The conditions in `remove_by_query` are combined with an **OR** operation
- All filtering conditions (`remove_by_column` and `remove_by_query`) are combined with **OR** logic
- Cells matching ANY of the specified values in ANY column OR ANY query will be removed
- Example: A cell will be excluded if it has `phase == "G1"` OR `is_cd14_mono == "true"` OR satisfies `random < 3`

## Input/Output

- **Input**: AnnData files in `.zarr` or `.h5ad` format
- **Output**: Filtered AnnData file (subset if `subset: true`, otherwise original structure with filter metadata)
