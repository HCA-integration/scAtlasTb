# The `metrics` module

This module computes metrics for quantitatively evaluating data integration (scib, Moran's I metrics).
The metrics measure either the biological conservation or the batch correction of a cellular representation.
Under the hood it parallelises the computation of each metric and allows for robust inclusion or exclusion of metrics without needing to recompute other metrics.
The final output are a dataframe of the different metric scores in the `AnnData` object and as a `TSV` file, as well as a funkyheatmap plot (see #output).

## Environments

The following environments are needed for the different metrics. Depending on which metrics you want to run, you do not need to install all environments.

- `scanpy` (default; CPU or GPU via `rapids_singlecell`)
- `scib` for metrics implemented in `scib`
- `scib_metrics` for metrics implemented in `scib-metrics`
- `pegasus` for `kbet_pg`
- `funkyheatmap` for summary heatmaps of results (needed in all cases)

You can find a complete mapping of all available metrics to their environments in `params.tsv` (see #metrics-parameters).

## Configuration

Configure the module under your dataset key using the `metrics` section. Common keys:

- `label`: biological label column in `.obs` (e.g. cell type) used for bio metrics
- `batch`: batch column in `.obs`
- `unintegrated`: slot for unintegrated counts/embedding (e.g. `layers/norm_counts` or `X`)
- `corrected`: slot for corrected representation when evaluating a single file
- `raw_counts`: slot for raw (unnormalised) counts, only needed by metrics that use gene sets
- `var_mask`: `.var` mask name to subset features (default `highly_variable`)
- `recompute_neighbors`: whether to recompute neighbors even if they exist in the AnnData object (default `false`). Will only be recomputed for non kNN outputs
- `clustering`: optional args for cluster-based metrics; supports `kwargs` and `overwrite` (see `clustering` module for more details)
- `gene_sets`: optional gene sets (dict of gene set name and gene list OR str with comma-separated keys of gene set mapping in `MARKER_GENES`) for gene-based metrics
- `covariates`: optional covariates (str or list) used by Moran's I and PC_regression metric
- `n_permutations`, `n_quantiles`: settings for gene-scoring utilities
- `weight_batch`: weight in [0, 1] of aggregated batch score used in funkyheatmap. The aggregated bio score gets a weight of 1 - `weight_batch`
- `metrics`: list of metric names to compute; must exist in `params.tsv`

### Example configuration

Below is an example of a metrics workflow for the test input file `data/pbmc68k.h5ad`.
Note that the file_id (under `input: metrics:`) can be parsed, if you set `overwrite_file_id` to `true`.

```yaml
output_dir: data/out
images: images

DATASETS:
  compare_reprs:
    input:
      metrics:
        # parsable names
        model1--example_param=a: data/pbmc68k.h5ad
        model2--example_param=b: data/pbmc68k.h5ad
        unintegrated: data/pbmc68k.h5ad
    metrics:
      unintegrated: layers/normcounts
      raw_counts: layers/counts
      corrected: X
      overwrite_file_id: true  # parse file_id, so that an additional column called `example_param` gets created in results table
      label: bulk_labels  # cell type label
      batch: batch
      output_type: full  # will use .X as representation for metrics
      metrics:  # list of metrics to be computed
        - nmi
        - ari
        - asw_label
        - asw_batch
```

### Cell representation type `output_type`

The type of a representation that should be used by the module needs to either be specified under `output_type` or in `.uns['output_type']`. If `.uns['output_type']` does not exist in an input file, the module will default to what is set in the module config, which is `'embed'` if nothing is set.

```yaml
...
DATASETS:
  ...
    metrics:
      ...
      output_type: embed # will be used as default when .uns['output_type'] is not set in file
```

The different `output_type` options are:

1. corrected graph (`knn`)
    + `.obsp['connectivities']`, `.obsp['distances']` integrated kNN graph returned by integration method
2. corrected embedding (`embed`)
    + `.obsm['X_emb']` integrated embedding returned by integration method
3. corrected features (`full`)
    + `.X` corrected feature counts returned by integration method
    + `.obsm['X_pca']` PCA on corrected feature counts


### Gene sets

The `gene_sets` parameters is only needed by the `morans_i_genes` and `pcr_genes` metrics.
If none of these metrics are needed, you also don't need to configure `gene_sets`.

```yaml
...

DATASETS:
  ...
    metrics:
      ...
      metrics:
        - pcr_comparison # does NOT use gene_sets
        - pcr # does NOT use gene_sets, but uses covariates
        - pcr_genes # uses gene_sets
        - morans_i_genes # gene_sets
      gene_sets:
        ifn_signature: [  # list can be written as object...
          "IRF7", "XAF1", "UBE2L6", "TRIM22", "STAT1",
          "SP110", "SAMD9L", "SAMD9", "PLSCR1", "PARP9",
          "OAS2", "OAS1", "MX2", "MX1", "LY6E",
          "ISG15", "IFIT3", "IFI6", "IFI44L", "IFI35",
          "HERC5", "EPSTI1", "EIF2AK2", "CMPK2", "BST2"
        ]
        platelets: # ... or as yaml list
          - GP1BB
          - ITGA2B
          - PF4
          - PPBP
          - TUBB1
```

Gene sets can also be provided as comma-separated strings of gene set keys that are defined under `MARKER_GENES`.

```yaml
...

DATASETS:
  ...
    metrics:
      ...
      gene_sets: ifn_signature, platelets

MARKER_GENES:
  ifn_signature: [
    "IRF7", "XAF1", "UBE2L6", "TRIM22", "STAT1",
    "SP110", "SAMD9L", "SAMD9", "PLSCR1", "PARP9",
    "OAS2", "OAS1", "MX2", "MX1", "LY6E",
    "ISG15", "IFIT3", "IFI6", "IFI44L", "IFI35",
    "HERC5", "EPSTI1", "EIF2AK2", "CMPK2", "BST2"
  ]
  platelets: ["GP1BB", "ITGA2B", "PF4", "PPBP", "TUBB1"]
```

## Output

The output directory are `<out_dir>/metrics` for the metrics and `<out_dir>/images/metrics` for the plots.

For each input file, there exists one output file containing an `AnnData` object in `<out_dir>/metrics/dataset~{dataset}/file_id~{file_id}.zarr`, where the metrics DataFrame is stored in `.uns['metrics']`.

In addition, merged result tables are written to:

- `results/metrics.tsv` (all datasets/files)
- `results/per_dataset/{dataset}_metrics.tsv`
- `results/per_batch/{batch}_metrics.tsv`
- `results/per_label/{label}_metrics.tsv`

Plots are written to `images/metrics/`:

- Barplots per all/dataset/file: `{metric}-barplot.png`
- Swarmplots per all/dataset/file: `{metric}-swarmplot.png`
- Funkyheatmap summaries: `all/per_dataset/per_batch/per_label/funky_heatmap.(pdf|tsv)`

E.g.

```tsv
metric      method   output_type   metric_type        score
nmi         scvi     embed         bio_conservation   0.876
asw_label   scvi     embed         bio_conservation   0.876
nmi         scanvi   embed         bio_conservation   0.605
asw_label   scanvi   embed         bio_conservation   0.605
```

# Contributing to the module

## Adding a new metric

1. Create a new script in `scripts`, choose a good name for the script.
2. Implement the metric according to the input/output specifications below.
3. Add a new entry in the `params.tsv` specifying the metric name, metric type and output types (#metrics-parameters)
4. Ensure that any additional dependencies are specified in the environment file used for calling the script
5. Add a new entry for the metric in `tests/config.yaml` and test the module.

If needed, this module can be extended to use different environment files for different metrics (analogous to the
integration module).

## Metric paramters

The metrics requirements are encoded in `params.tsv`:

- `metric`, `metric_type` (`bio_conservation` or `batch_correction`)
- `output_types` (`knn`, `embed`, `full`) the output types (representation types) for which the metric will be computed. If an input file has an output type that does not apply to a metric, no metric will be computed and `np.nan` will be return instead.
- `input_type` the representation type used by the metric after computing e.g. PCA or kNN in the prepare step. E.g. if a metric requires a low-dimensional embedding, it is not limited to embedding integrations, but can also be applied on the PCA of a feature output integration.
- `comparison` whether the metric compares to the unintegrated representation
- `needs_clustering` whether the metrics requires Leiden clustering for multiple resolutions [0.1 - 2.0]
- `use_covariate` whether the metric uses custom covariates. Currently all covariates are assumed to contain biological signal and the resulting metric is considered to be a bio conservation metric
- `use_gene_set` whether the metric uses a gene set
- `env` which environment to use
- `resources` (cpu/gpu)
- `threads` number of threads to provide for the metric call

These features inform the `scripts/metrics/run.py` on the required preparation steps as what parts of the object need to be read for efficient computation.

## Testing

Activate the snakemake environment and call `test/run_test.sh` with run specific Snakemake parameters.

```
conda activate snakemake
bash test/run_test.sh -n
bash test/run_test.sh -c
```

The script will call the pipeline and run a test script.
If the input files don't yet exist, the test script might throw a `FileNotFoundError`.
You can ignore it, if you haven't yet executed the pipeline completely and therefore don't yet have all outputs.
