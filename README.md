# Single Cell Atlasing Toolbox 🧰

**Toolbox of Snakemake pipelines for easy-to-use analyses and benchmarks for building integrated atlases**


This toolbox provides multiple modules that can be easily combined into custom workflows that leverage the file management of [Snakemake](https://snakemake.readthedocs.io/en/v7.31.1/).
This allows for an efficient and scalable way to run analyses on large datasets that can be easily configured by the user.

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

💖 Beautiful, right? Chek out the [documentation]() to learn how to set up your own workflow!


## ⚙️ Advanced configuration
<a name="advanced-configuration"></a>

### Set defaults

You can set module-specific defaults that will be used for all tasks (under `configs['DATASETS']`), if the parameters have not been specified for those tasks.
This can shorten the configuration file, make it more readable and help avoid misconfiguration if you want to reuse the same configurations for multiple tasks.

Under the `defaults` directive, you can set the defaults in the same way as the task-specific configuration.

<details>
  <summary>Example defaults for modules</summary>

  ```yaml
  defaults:
    preprocessing:
      highly_variable_genes:
        n_top_genes: 2000
      pca:
        n_comps: 50
      assemble:
        - normalize
        - highly_variable_genes
        - pca
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
    metrics:
      unintegrated: layers/norm_counts
      batch: batch
      label: bulk_labels
      methods:
        - nmi
        - graph_connectivity
  ```

</details>

Additionaly to the module defaults, you can set which datasets you want to include in your workflow, without having to remove or comment out any entries in `configs['DATASETS']`.

```yaml
defaults:
...
  datasets:
  # list of dataset/task names that you want your workflow to be restricted to
    - test
    - test2
```


### Automatic environment management
Snakemake supports automatically creating conda environments for each rule.

```yaml
env_mode: from_yaml
```

You can trigger Snakemake to install all environments required for your workflow in advance by adding the following parameters

```
<snakemake_cmd> --use-conda --conda-create-envs-only --cores 1
```

### Snakemake profiles

Snakemake profiles help you manage the many flags and options of a snakemake command in a single file, which will simplify the Snakemake call considerably.
The toolbox provides some example Snakemake profiles are provided under `.profiles`, which you can copy and adapt to your needs.

To use a profile e.g. the local profile, add `--profile .profiles/<profile_name>` to your Snakemake command.
You can read more about profiles in [Snakemake's documentation](https://snakemake.readthedocs.io/en/v7.31.1/executing/cli.html#profiles).

### Cluster execution

Snakemake supports scheduling rules as jobs on a cluster.
If you want your workflow to use your cluster architecture, create a Snakemake profile under `.profiles/<your_profile>/config.yaml`.

<details>
  <summary>Example profile for SLURM</summary>

  Adapted from https://github.com/jdblischak/smk-simple-slurm
  ```yaml
  cluster:
    mkdir -p logs/{rule} &&
    sbatch
      --partition={resources.partition}
      --qos={resources.qos}
      --gres=gpu:{resources.gpu}
      --cpus-per-task={threads}
      --mem={resources.mem_mb}
      --job-name={rule}
      --output=logs/%j-{rule}.out
      --parsable
  default-resources:
    - partition=cpu
    - qos=normal
    - gpu=0
    - mem_mb=90000
    - disk_mb=20000
  restart-times: 0
  max-jobs-per-second: 10
  max-status-checks-per-second: 1
  local-cores: 1
  latency-wait: 30
  jobs: 20
  keep-going: True
  rerun-incomplete: True
  printshellcmds: True
  scheduler: ilp
  use-conda: True
  cluster-cancel: scancel
  rerun-triggers:
    - mtime
    - params
    - input
    - software-env
    - code
  show-failed-logs: True
  ```
</details>

In order to specify the actual cluster parameters such as memory requirements, nodes or GPU, you need to specify the resources in your config file.
The toolbox requires different settings for CPU and GPU resources.

```yaml
resources:
  cpu:
    partition: cpu
    qos: normal
    gpu: 0
    mem_mb: 100000
  gpu:
    partition: gpu
    qos: normal
    gpu: 1
    mem_mb: 100000
```

If you don't have have GPU nodes, you can configure the gpu resources to be the same as the cpu resources.

You can find detailed information on cluster execution in the [Snakemake documentation](https://snakemake.readthedocs.io/en/v7.31.1/executing/cluster.html).

## 🛠️ Troubleshooting
<a name="troubleshooting"></a>

### Conda environment installation fails
For some environments (e.g. `envs/rapids_singlecell.yaml`) you need to set the channel priority to flexible in order for conda to properly resolve the environment.

```
conda config --set channel_priority flexible
```

After setting the flag and calling the installation command, the environment should resolve.

### Working with GPUs

Some scripts can run faster if their dependencies are installed with GPU support.
Currently, whether the GPU version of a package with GPU support is installed, depends on the architecture of the system that you install you **install** the environment on.
If you work on a single computer with GPU, GPU-support should work out of the box.
However, if you want to your code to recognize GPUs when working on a cluster, you need to make sure you install the conda environments from a node that has access to a GPU.

Environments that support GPU are:

* `rapids_singlecell` (only installs when GPU is available)
* `scarches`
* `scib_metrics`
* `scvi-tools`

If you have already installed a GPU environment on CPU, you need to remove and re-install it on node with a GPU.

```
conda env remove -n <env_name>
mamba env create -f envs/<env_name>.yaml
```

In case you are working with `env_mode: from_yaml`, gather the environment name from the Snakemake log, remove the environment manually.
The next time you call your pipeline again, Snakemake should automatically reinstall the missing environment.

### Working with CPUs only

If your system doesn't have any GPUs, you can set the following flag in your config.

```yaml
use_gpu: false
```

This will force Snakemake to use the CPU versions of an environment.

### FAQs

Below are some scenarios that can occur when starting with the pipeline.
If you have any additional questions or encounter any bugs, please open up a [github issue](https://github.com/HCA-integration/scAtlasTb/issues).
If you want to contribute to improving the pipeline, check out the [contribution guidelines](CONTRIBUTING.md).

#### I configured my pipeline and the dry run doesn't fail, but it doesn't want to run the modules I configured. What do I do?

This likely happens when you don't specify which rule you want Snakemake to run. By default, Snakemake will try create a visualisation of the modules you configured. If you want it to run the modules themselves, you will need to add the rule name with your Snakemake command. For each rule, there is a `<module>_all`, but you can view all possible rules through `snakemake -l`

