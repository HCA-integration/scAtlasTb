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
      subset: true ## whether to subset the file according to the filters
      remove_by_column:
        donor: # column name
          # entries from column to exclude (will be treated as String)
          - C19-CB-0008 # value not in data
          - Normal 1
          - nCoV 1
          - Flu 1
      keep_by_column:
        sex:
          - female
      remove_by_query:
        # pandas query strings for complex filtering conditions
        - 'random < 3'
      keep_by_query:
        # pandas query strings for cells to keep
        - 'disease == "influenza"'
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

- **`write_copy`**: Boolean flag for output format when subsetting
  - `true`: Writes a full copy of the data
  - `false` (default): Writes linked zarr format when possible

- **`remove_by_column`**: Dictionary defining filtering conditions
  - Keys: Column names in the AnnData observation metadata (`adata.obs`)
  - Values: List of entries to exclude from that column
  - All entries are treated as strings for comparison

- **`keep_by_column`**: Dictionary defining columns and entries to retain
  - Keys: Column names in the AnnData observation metadata (`adata.obs`)
  - Values: List of entries to keep from that column
  - All entries are treated as strings for comparison

- **`remove_by_query`**: List of pandas query strings for complex filtering conditions
  - Each query string follows pandas DataFrame query syntax
  - Cells matching these queries will be removed

- **`keep_by_query`**: List of pandas query strings for cells to retain
  - Each query string follows pandas DataFrame query syntax
  - Only cells matching these queries will be kept

## Behavior

- All `remove_by` filtering conditions are combined with **AND** logic - a cell must pass ALL filters to be retained
  - Within `remove_by_column`: cells are excluded if they match ANY value in ANY specified column
  - Within `remove_by_query`: cells are excluded if they match ANY of the query conditions
- All `keep_by` filtering conditions are combined with **OR** logic - a cell only needs to pass ONE filter to be retained
  - Within `keep_by_column`: cells are retained if they match ANY value in ANY specified column
  - Within `keep_by_query`: cells are retained only if they match ANY of the query conditions
- The final mask keeps cells that:
  1. Do NOT match any values in `remove_by_column` columns, AND
  2. Do NOT match any `remove_by_query` conditions, AND  
  3. DO match all `keep_by_column` conditions (if specified)
  4. DO match all `keep_by_query` conditions (if specified)

## Input/Output

- **Input**: AnnData files in `.zarr` or `.h5ad` format
- **Output**: Filtered AnnData file (subset if `subset: true`, otherwise original structure with filter metadata)
- A `filtered` column is added to `adata.obs` indicating which cells passed the filters
