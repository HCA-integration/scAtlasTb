# Label Harmonisation

## Testing

Activate the snakemake environment and call `test/run_test.sh` with run specific Snakemake parameters.

```
conda activate snakemake
bash test/run_test.sh -n  # dry run
bash test/run_test.sh -c2  # actual run with max 2 cores
```

## Input

### Config file

```yaml
DATASETS:
  dataset_name:  # can replace with a representative name
    input:
      relabel:
        file_1: test/input/pbmc68k.h5ad
        file_2: test/input/pbmc68k.h5ad
    relabel:
      new_columns:  # mapping setup for new columns
        file:  test/input/mapping.tsv  # TSV file with cell label mapping
        order:  # order of which label column should be matched to which
          - cell_type  # first entry MUST be an existing column in the anndata object
          - harmonized_label  # first remapping level (mapping 'cell_type' to 'harmonized_label')
          - lineage  # second remapping level (mapping 'harmonized_label' to 'lineage')
      merge_columns:  # mapping setup for merging existing columns
        file:  test/input/merge_test.tsv  # TSV file with cell label mapping
        sep: '-'  # separator for merging columns
      
      selective_update:
        base_column: bulk_labels
        new_column: bulk_labels_new
        update_map:
          louvain:
            '6': Dendritic
            '8': 'CD19+ B'
```

### Mapping files

#### New columns

The new column mapping file must be in TSV format with at least one column in the anndata object.
The columns must include all columns that are specified in the `order` list.

Example for `test/input/pbmc68k.h5ad`:

```
bulk_labels	lineage
CD14+ Monocyte	myeloid
Dendritic	lymphoid
CD56+ NK	lymphoid
CD4+/CD25 T Reg	lymphoid
```

#### Merge existing columns

The merge column mapping file must be in TSV format with the following columns:

* `file_id`: input file id that is used in the pipeline to match which file to add the new column to
* `column_name`: new column name
* `columns`: columns to merge, separated by `,` (no whitespaces)

Example for `test/input/pbmc68k.h5ad`:

```
file_id	column_name	columns
file_1	cc_joint	S_score,G2M_score
file_1	ct_joint	bulk_labels,louvain
file_2	counts_joint	n_genes,n_counts
```

### Selective update of an existing column: `selective_update`

The selective update allows you to create a new column based on an existing column and a mapping of values.

```yaml
...
      selective_update:
        base_column: bulk_labels # column for which values should be updated
        new_column: bulk_labels_new # new column to be created
        update_map: # mapping of values to be updated based on column
          louvain: # column to take values from 
            '6': Dendritic # selected values for determing cells that should be remapped
            '8': 'CD19+ B'
```

* `base_column`: the starting column of which the values should be modified
* `new_column`: the name of the new column to be created. If it is not defined, it will be the same as `base_column`.
* `update_map`: a mapping of values that should be updated. The keys are the columns in the AnnData object, and the values are dictionaries with the values to be updated. Alternatively, the values can be a file path to a TSV file with the mapping.

Below is an example where one update enrty is defined by a TSV file instead of a dictionary:

```yaml
...
      selective_update:
        base_column: bulk_labels
        new_column: lineage
        update_map:
          lineage: test/input/mapping_test.tsv # mapping defined in a TSV file
          louvain:
            '6': myeloid
            '8': lymphoid
```

When using a mapping TSV file for `update_map`, the file must contain the `base_column` (in this example `bulk_labels`) as well as the `update_map` key column (in this example `lineage`).
The following example `test/input/mapping_test.tsv` shows a valid TSV definition of the `update_map` key `lineage`: 

```
bulk_labels	lineage
CD14+ Monocyte	myeloid
Dendritic	myeloid
CD56+ NK	lymphoid
CD4+/CD25 T Reg	lymphoid
CD19+ B	lymphoid
CD8+ Cytotoxic T	lymphoid
CD4+/CD45RO+ Memory	lymphoid
CD8+/CD45RA+ Naive Cytotoxic	lymphoid
CD4+/CD45RA+/CD25- Naive T	lymphoid
CD34+	HSCs
```


## Testing

Activate the snakemake environment and call `test/run_test.sh` with run specific Snakemake parameters.

```
conda activate snakemake
bash test/run_test.sh -n  # dry run
bash test/run_test.sh -c2  # actual run with max 2 cores
```
