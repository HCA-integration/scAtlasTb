# 🚀 Quickstart

Activate the snakemake environment

```
conda activate snakemake
```

Before running the pipeline, install the required conda environments.
`install_all_environments.sh` installs every environment in the repository — most of which you won't need for this quickstart. Instead, use the quickstart-specific installer:

```
bash envs/install_quickstart_environments.sh
```

This installs 9 out of 20+ environments: just what `run_example.sh` actually needs.

> 📝 **Note:** This only needs to be done once. See the [installation guide](installation.md) for more details on environment management, including how to install only the environments you need.

## What `run_example.sh` does

`run_example.sh` calls Snakemake with `configs/quickstart.yaml`. The outline of the configuration is as follows:

```yaml
input:
    ...
preprocessing:
    ...
integration:
    ...
    methods:
        unintegrated:
        bbknn:
        combat:
        harmony_pytorch:
        harmonypy:
        scanorama:
        scvi:
        drvi:
        sysvi:
metrics:
    # the metrics listed here are run in this pipeline

```

`preprocessing` reads the raw file, `integration` takes its output, and `metrics` takes the integration output — a linear chain.

## Running the pipeline

Call the pipeline with `-n` for a dry run and `-q` for reduced output.
Here's the command for running preprocessing, integration and metrics

```
bash run_example.sh preprocessing_all integration_all metrics_all -nq
```

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

If the dryrun was successful, you can let Snakemake compute the different steps of the workflow with e.g. 10 cores:

```
bash run_example.sh preprocessing_all integration_all metrics_all -c 10
```

> You have now successfully called the example pipeline! 🎉  
> Read on to learn how to [configure your own workflow](configure_workflow.rst).
