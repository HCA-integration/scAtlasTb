# Doublet detection

This module runs per-batch doublet detection using one or more callers (currently `doublet_detection` and `scrublet`). The pipeline runs doublet callers within each library/batch and writes scores and predictions back into the AnnData object.

## Environments

The following environments are useful for running the module. Install only the ones you need.

- [`scanpy`](https://github.com/HCA-integration/scAtlasTb/blob/main/envs/scanpy.yaml)
- (caller-specific envs may be required for `scrublet` or `doublet_detection`)

# Configuration

```yaml
DATASETS:
  Lee2020:
    input:
      doublets: 
        Lee2020: test/input/load_data/harmonize_metadata/Lee2020.zarr
    doublets:
      counts: X
      batch: donor
      chunk_size: 10_000

  test:
    input:
      doublets:
        test: test/input/pbmc68k.h5ad
        test2: test/input/pbmc68k.h5ad
    doublets:
      counts: layers/counts
      methods:
        - scrublet
        - doubletdetection

defaults:
  datasets:
    - test
    - Lee2020
```

* `counts`: Slot in anndata that contains raw (unnormalized) counts. Examples: `X`, `raw/X`, or `layers/<layer_name>`.
* `batch`: Column in `obs` that contains batch/library IDs. Doublet detection is executed separately per batch.
* `chunk_size`: Number of cells used to group batches for more efficient parallel processing. Default: `100_000`.
* `methods`: Optional list of doublet callers to run for this dataset (e.g. `['scrublet', 'doubletdetection']`). If omitted, pipeline defaults apply.

> Note: `counts` is resolved relative to the input object (e.g., anndata.X or anndata.layers). Methods are run per-batch; choose `batch` to reflect library-level grouping so cross-library doublets are not considered.

## Output

Results are written into the AnnData `.obs` and to summary files/plots:

* `<out_dir>/doublets/dataset~<datasets>/file_id~<file_id>.zarr` — Updated AnnData with caller-specific score and prediction columns in `.obs`, for example:
  - scrublet:
    - `doublet_score`
    - `predicted_doublet`
  - doublet_detection:
    - `doublet_score_doubletdetection`
    - `doubletdetection_prediction`
* `<image_dir>/doublets/dataset~<datasets>/file_id~<file_id>/doublet_summary.tsv` — Summary table with per-batch and per-method statistics.
* `<image_dir>/doublets/dataset~<datasets>/file_id~<file_id>/plots/` — Visualizations of score distributions and per-batch results (e.g., histograms, per-batch summary plots).

If additional callers are enabled, corresponding score/prediction columns will be added to `.obs` using method-specific names to avoid collisions.