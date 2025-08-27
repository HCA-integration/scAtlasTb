# Doublet detection

Currently this module implements [`doublet_detection`](https://github.com/JonathanShor/DoubletDetection) and [`scrublet`](https://github.com/swolock/scrublet), but in future, there will be more flexibility on which method the user wants to run.

## Config

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
```

* `counts`: Slot in anndata that contains raw (unnormalized) counts. Can be e.g. `X`, `raw/X`, or `layers/<layer_name>`.
* `batch`: Column in `obs` that contains batch information. This should be ideally be the library ID, since doublets only occur within 1 library, but not across libraries. The doublet algorithm will run separately for each batch.
* `chunk_size`: Number of cells to process in each chunk when the number of batches is large. This is used internally to more efficiently parallel doublet runs per batch, where as many batches are grouped until the `chunk_size` is reached. Default: 100_000.

## Output
The output is stored in the `obs` slot of the anndata object, with the following columns:

* `doublet_score`: The doublet score predicted by `scrublet`
* `predicted_doublet`: Whether the cell is predicted to be a doublet by `scrublet` (default threshold used)
* `doublet_score_doubletdetection`: The doublet score predicted by `doublet_detection`
* `doubletdetection_prediction`: Whether the cell is predicted to be a doublet by `doublet_detection`