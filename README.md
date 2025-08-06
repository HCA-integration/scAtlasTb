# Single Cell Atlasing Toolbox 🧰

[![Documentation][badge-docs]][documentation]

[badge-tests]: https://img.shields.io/github/actions/
[badge-docs]: https://img.shields.io/readthedocs/scAtlasTb-utils

**Toolbox of Snakemake pipelines for easy-to-use analyses and benchmarks for building integrated atlases**

This toolbox provides multiple modules that can be easily combined into custom workflows that leverage the file management of [Snakemake](https://snakemake.readthedocs.io/en/v7.31.1/).
This allows for an efficient and scalable way to run analyses on large datasets that can be easily configured by the user.

## Getting started

Please refer to the [documentation][].


## 🧰 Which Modules does the Toolbox Support?

The modules are located under `workflow/` and can be run independently or combined into a more complex workflow.

| Module                 | Description                                                               |
|------------------------|---------------------------------------------------------------------------|
| `load_data`            | Loading datasets from URLs and converting them to AnnData objects         |
| `exploration`          | Exploration and quality control of datasets                               |
| `batch_analysis`       | Exploration and quality control of batches within datasets                |
| `qc`                   | Semi-automaetd quality control of datasets using [sctk AutoQC](https://teichlab.github.io/sctk/notebooks/automatic_qc.html) |
| `doublets`             | Identifying and handling doublets in datasets                             |
| `merge`                | Merging datasets                                                          |
| `filter`               | Filtering datasets based on specified criteria                            |
| `subset`               | Creating subsets of datasets                                              |
| `relabel`              | Relabeling data points in datasets                                        |
| `split_data`           | Splitting datasets into training and testing sets                         |
| `preprocessing`        | Preprocessing of datasets (normalization, feature selection, PCA, kNN graph, UMAP) |
| `integration`          | Running single cell batch correction methods on datasets                  |
| `metrics`              | Calculating [scIB metrics](https://scib.readthedocs.io/en/latest/), mainly for benchmarking of integration methods  |
| `label_harmonization`  | Providing alignment between unharmonized labels using [CellHint](https://cellhint.readthedocs.io/en/latest) |
| `label_transfer`       | Transfer annotations of annotated cells to annotated cells e.g. via majority voting |
| `sample_representation`| Methods for aggregating cells to sample level e.g. pseudobulk             |
| `collect`              | Collect multiple input anndata objects into a single anndata object       |
| `uncollect`            | Distribute slots of an anndata object to multiple anndata objects         |
| `cell_type_prediction` | Predict cell types from reference model e.g. celltypist                   |

## 👀 TL;DR What does a full workflow look like?

The heart of the configuration is captured in a YAML (or JSON) configuration file.
Here is an example of a workflow configuration in `configs/example_config.yaml` containing the `preprocessing`, `integration` and `metrics` modules:

```yaml
output_dir: data/out
images: images

os: intel
use_gpu: true

DATASETS:

  my_dataset: # custom task/workflow name

    # input specification: map of module name to map of input file name to input file path
    input:
      preprocessing:
        file_1: data/pbmc68k.h5ad
        # file_2: ... # more files if required
      integration: preprocessing # all outputs of module will automatically be used as input
      metrics: integration
    
    # module configuration
    preprocessing:
      highly_variable_genes:
        n_top_genes: 2000
      pca:
        n_comps: 50
      assemble:
        - normalize
        - highly_variable_genes
        - pca
    
    # module configuration
    integration:
      raw_counts: raw/X
      norm_counts: X
      batch: batch
      methods:
        unintegrated:
        scanorama:
          batch_size: 100
        scvi:
          max_epochs: 10
          early_stopping: true

    # module configuration
    metrics:
      unintegrated: layers/norm_counts
      batch: batch
      label: bulk_labels
      methods:
        - nmi
        - graph_connectivity
```

Which allows you to call the pipeline as follows:

```
snakemake --configfile configs/example_config.yaml --snakefile workflow/Snakefile --use-conda -nq
```

giving you the following dryrun output:

```
Job stats:
job                                    count
-----------------------------------  -------
integration_all                            1
integration_barplot_per_dataset            3
integration_benchmark_per_dataset          1
integration_compute_umap                   6
integration_plot_umap                      6
integration_postprocess                    6
integration_prepare                        1
integration_run_method                     3
preprocessing_assemble                     1
preprocessing_highly_variable_genes        1
preprocessing_normalize                    1
preprocessing_pca                          1
total                                     31

Reasons:
    (check individual jobs above for details)
    input files updated by another job:
        integration_all, integration_barplot_per_dataset, integration_benchmark_per_dataset, integration_compute_umap, integration_plot_umap, integration_postprocess, integration_prepare, integration_run_method, preprocessing_assemble, preprocessing_highly_variable_genes, preprocessing_pca                                                                                             
    missing output files:
        integration_benchmark_per_dataset, integration_compute_umap, integration_postprocess, integration_prepare, integration_run_method, preprocessing_assemble, preprocessing_highly_variable_genes, preprocessing_normalize, preprocessing_pca

This was a dry-run (flag -n). The order of jobs does not reflect the order of execution.
```

💖 Beautiful, right? Chek out the [documentation][] to learn how to set up your own workflow!

## Release notes

See the [changelog][].

## Contact

If you found a bug, please use the [issue tracker][].

## Citation

> t.b.a

[issue tracker]: https://github.com/HCA-integration/scAtlasTb/issues
[documentation]: https://scatlastb.readthedocs.io
[changelog]: https://scatlastb-utils.readthedocs.io/en/latest/changelog.html
