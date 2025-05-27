# Single Cell Atlasing Toolbox :toolbox:

**Toolbox of Snakemake pipelines for easy-to-use analyses and benchmarks for building integrated atlases**

- [:rocket: Getting Started](#getting-started)
- [:gear: Configure Your Workflow](#configure-your-workflow)
- [:gear: Advanced Configuration](#advanced-configuration)
- [:hammer_and_wrench: Troubleshooting](#troubleshooting)

This toolbox provides multiple modules that can be easily combined into custom workflows that leverage the file management of [Snakemake](https://snakemake.readthedocs.io/en/v7.31.1/).
This allows for an efficient and scalable way to run analyses on large datasets that can be easily configured by the user.

## :toolbox: Which Modules does the Toolbox Support?

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

## :eyes: TL;DR What does a full workflow look like?

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

:sparkling_heart: Beautiful, right? Read on to learn how to set up your own workflow!

## :rocket: Getting started
<a name="getting-started"></a>

### 1. Installation

#### Clone the repository

Depending on whether you have set up SSH or HTTPS with [PAT](https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/managing-your-personal-access-tokens), you can clone the repository with

SSH:
```
git clone git@github.com:HCA-integration/scAtlasTb.git
```

HTTPS:
``` clone
git clone https://github.com/HCA-integration/scAtlasTb.git
```

#### Requirements

* Linux (preferred) or MacOS on Intel (not rigorously tested, some bioconda dependencies might not work)
* conda e.g. via [miniforge](https://github.com/conda-forge/miniforge)(recommended) or [miniconda](https://docs.anaconda.com/free/miniconda/index.html)


The modules are tested and developed using task-specific conda environments, which should be quick to set up when using [mamba](https://mamba.readthedocs.io).

> :memo: **Note** If you use conda, but have never used mamba, consider installing the mamba package into your base environment and use it for all installation commands.
You can still replace all mamba commands with conda commands if you don't want to install mamba.

#### Install conda environments

The different parts of the workflow (modules, rules) require specific conda environments.
The simplest way to install all environments is to run the following script:

```
bash envs/install_all_environments.sh -h # help message for customization
bash envs/install_all_environments.sh
```

> :memo: **Notes**
> 1. The script will create new environments for each file in the `envs` directory if they don't yet exist and update any pre-existing environments.
> 2. The environment names correspond the their respective file names and are documented under the `name:` directive in the `envs/<env_name>.yaml` file.
> 3. If an environment creation fails, the script will skip that environment and you might need to troubleshoot the installation manually.
> 4. Some environments require the channel priority to be set to `flexible`.
> If your installation command fails, try setting `conda config --set channel_priority flexible` before restarting the command.

If you know you only need certain environments (you can get that information from the README of the module you intend to use), you can install that environment directly.
You will at least require the snakemake environment.

```
mamba env create -f envs/snakemake.yaml
mamba env create -f envs/<env_name>.yaml
```

### 2. Call the example pipeline

Activate the snakemake environment

```
conda activate snakemake
```

Call the pipeline with  `-n` for a dry run and `-q` for reduced output.
Here's the command for running preprocessing, integration and metrics

```
bash run_example.sh preprocessing_all integration_all metrics_all -nq
```

If the dryrun was successful, you can let Snakemake compute the different steps of the workflow with e.g. 10 cores:

```
bash run_example.sh preprocessing_all integration_all metrics_all -c 10
```

> You have now successfully called the example pipeline! :tada:
> Read on to learn how to configure your own workflow.

## :gear: Configure Your Workflow
<a name="configure-your-workflow"></a>

Configuring your workflow requires configuring global settings as well as subworkflows consisting of modules.
The global configuration allows you to set output locations, computational resources and other settings that are used across all modules, while module settings affect the behaviour of a module in the scope of a given task

> :memo: **Note** The recommended way to manage your workflow configuration files is to save them outside of the toolbox directory in a directory dedicated to your project. That way you can guarantee the separatation of the toolbox and your own configuration.

You can find example configuration files under `configs/`.

### 1. Global configuration: Output settings

You can specify pipeline output as follows.
Intermediate and large files will be stored under `output_dir`, while images and smaller outputs that are used for understanding the outputs will be stored under `images`.
If you use relative paths, you need to make them relative to where you call the pipeline (not the config file itself).
The directories will be created if they don't yet exist.

```yaml
# Note: relative paths must be relative to the project root, not the directory of the config file.
output_dir: data/out
images: images
```

Another setting is the output file pattern map.
By default, the final output pattern of a rule follows the pattern of
`<out_dir>/<module>/<wildcard>~{<wildcard>}/<more wildcards>.zarr`.
For some modules the final output pattern differs from that default and needs to be specified explicitly in the `output_map`.
In future, this shouldn't be necessary.

```yaml
output_map:
  sample_representation: data/out/sample_representation/dataset~{dataset}/file_id~{file_id}/pseudobulk.h5ad
  subset: data/out/subset/dataset~{dataset}/file_id~{file_id}/by_sample.zarr
  pca: data/out/preprocessing/dataset~{dataset}/file_id~{file_id}/pca.zarr
  neighbors: data/out/preprocessing/dataset~{dataset}/file_id~{file_id}/neighbors.zarr
  preprocessing: data/out/preprocessing/dataset~{dataset}/file_id~{file_id}/preprocessed.zarr
  metrics: data/out/metrics/results/per_dataset/{dataset}_metrics.tsv
```

The default output settings under `configs/outputs.yaml` should work out of the box.

### 2. Global configuration: Computational settings

Depending on the hardware you have available, you can configure the workflow to make use of them.
If you have a GPU, you can set `use_gpu` to `true` and the pipeline will try to use the GPU for all modules that support it.
The same applies if you have an Intel CPU.
In the backend, this affects which conda environment Snakemake uses, whenever hardware-accelerated environments are specified in a rule.

```yaml
os: intel
use_gpu: true
```

### 3. Input configuration

You can select and combine modules to create a custom workflow by specifying the input and module configuration in a YAML file.
Each instance of a workflow needs a unique task name and it can take any number of inputs consist of modules.

```yaml
DATASETS: # TODO: rename to TASKS

  my_dataset: # custom task/workflow name
    # input specification: map of module name to map of input file name to input file path
    input:
      preprocessing:
        file_1: data/pbmc68k.h5ad
        # file_2: ... # more files if required
      integration: preprocessing # all outputs of module will automatically be used as input
      metrics: integration

  another_dataset:
    ...
 ```

> :warning: **Warning** There can only be one instance of a module as a key in the input mapping (in the backend this is a dictionary). But you can reuse the same module output as input for multiple other modules. The order of the entries in the input mapping doesn't matter. 

### 4. Module configuration

You can configure the behaviour of each module by specifying their parameters under the same dataset name.
 ```yaml
DATASETS:
  my_dataset:
    input:
      ...

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

Each module has a specific set of parameters that can be configured.
Read more about the specific parameters in the README of the module you want to use.

### 5. Create a wrapper script (recommended)

Next, you can create a wrapper script that will call the pipeline with the correct profile and configuration file(s).
This way, it is easier to call the pipeline and you can avoid having to remember all the flags and options.

Below is an example of a wrapper script that you can use to call the pipeline.

```bash
#!/usr/bin/env bash
set -e -x

# assuming that the toolbox is on the same level of the directory you're calling the script from
pipeline="$(realpath ../scAtlasTb)"   # adjust depending on location of wrapper script

snakemake \
  --configfile <my_config_file>.yaml \
  --snakefile $pipeline/workflow/Snakefile \
  --use-conda \
  --rerun-incomplete \
  --keep-going \
  --printshellcmds \
    $@
```

You must set the flag `--use-conda` to ensure that Snakemake uses the conda environments specified in the rules.
If your config file becomes very big, you can split the workflows into separate config files and include them to the `configfile` in the wrapper script.

```bash
#!/usr/bin/env bash
set -e -x

pipeline="$(realpath ../scAtlasTb)"

snakemake \
  --configfile \
    <my_config_file>.yaml \
    <my_config_file_2>.yaml \
  ...
```

> :bulb: **Tip** Check out the [snakemake documentation](https://snakemake.readthedocs.io/en/v7.31.1/executing/cli.html) for more commandline arguments.

### 6. Call the Snakemake pipeline

Before running the pipeline, you need to activate your Snakemake environment.

```
conda activate snakemake
```

<details>
  <summary>How does Snakemake work?</summary>

  > The general command for running a pipeline is:
  >
  > ```
  > snakemake <snakemake args>
  > ```
  >
  > The most relevant snakemake arguments are:
  > 
  > + `-n`: dryrun
  > + `--use-conda`: use rule-specific conda environments to ensure all dependencies are met
  > + `-c`: maximum number of cores to be used
  > + `--configfile`: specify a config file to use. The overall workflow already defaults to the config file under `configs/config.yaml`
  
  >  For more information on how Snakemake works, please refer to [Snakemake's extensive documentation](https://snakemake.readthedocs.io/en/v7.31.1/index.html).

</details>

#### First dry run

When you execute the script (say, we call it `run_pipeline.sh`), you can treat it like a snakemake command and add any additional snakemake arguments you want to use.
A dryrun would be:

```
bash run_pipeline.sh -n
```

This will show you what Snakemake wants to run.
Without specifying any rule, the default rules that the pipeline will request are `common_dag` and `common_rulegraph`.
You can ignore these for now.

```
...
Building DAG of jobs...
Job stats:
job                 count
----------------  -------
all                     1
common_dag              1
common_rulegraph        1
total                   3
```

#### List all available rules

The pipeline will only run the target that you explicitly tell it to run.
A target can be either the name of a Snakemake rule or a file that can be generated by Snakemake (as defined by the Snakefiles).
You can list all possible rules with:

```
bash run_pipeline.sh -l
```

Which should give you something like this:

```
all                            
batch_analysis_all             
batch_analysis_batch_pcr 
batch_analysis_collect           
batch_analysis_dependency_graph
batch_analysis_determine_covariates
batch_analysis_plot        
clustering_all            
clustering_cluster            
clustering_compute_neighbors
clustering_compute_umap            
clustering_dependency_graph 
clustering_merge
...
split_data_all
split_data_dependency_graph
split_data_link
split_data_split
subset_all
subset_dependency_graph
subset_subset
```

All the rules ending with `_all` are callable, i.e. you can use them to specify that their workflow should be run.
The rest are needed by the pipeline, but can't be called by the user, you can just ignore them.

#### Specify which workflow/rule you want to run

Given the [config above](#example_config), you can call the integration workflow by specifying the `integration_all` target:

```
bash run_pipeline.sh integration_all -n
```

This should list all the rules with details such as inputs, outputs and parameters, as well as the following summary:

```
...

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

Notice, that this also includes preprocessing rules that were defined in the config as input to the integration.
Since Snakemake only computes the files that are direct dependencies of the target that you define, the workflow does not include preprocessing-specific rules are not used by the integration module.
If you want to include all preprocessing rules, you need to include it in the command:

```
bash run_pipeline.sh preprocessing_all integration_all -n
```

Following the same principle, you can call the metrics by including the `metrics_all` rule to the target list:

```
bash run_pipeline.sh preprocessing_all integration_all metrics_all -n
```

#### Execute the workflow

If you are happy with the dryrun, you can dispatch the workflow by specifying the number of cores you want to provide for the pipeline.

```
bash run_pipeline.sh preprocessing_all integration_all metrics_all -c 10
```

You can also use [Snakemake profiles](#snakemake_profiles) to dispatch your pipeline with extra configurations such as Snakemake presets or [cluster execution](#cluster_execution).

> You have now successfully set up and configured your pipeline!
> Give it a spin and feel free to edit the configs to your custom workflow! :tada:


## :gear: Advanced configuration
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

## :hammer_and_wrench: Troubleshooting
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

