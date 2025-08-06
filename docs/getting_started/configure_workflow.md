# 🔧 Configure Your Workflow
<a name="configure-your-workflow"></a>

Configuring your workflow requires configuring global settings as well as subworkflows consisting of modules.
The global configuration allows you to set output locations, computational resources and other settings that are used across all modules, while module settings affect the behaviour of a module in the scope of a given task

> 📝 **Note** The recommended way to manage your workflow configuration files is to save them outside of the toolbox directory in a directory dedicated to your project. That way you can guarantee the separatation of the toolbox and your own configuration.

You can find example configuration files under `configs/`.

## 1. Global configuration: Output settings

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

## 2. Global configuration: Computational settings

Depending on the hardware you have available, you can configure the workflow to make use of them.
If you have a GPU, you can set `use_gpu` to `true` and the pipeline will try to use the GPU for all modules that support it.
The same applies if you have an Intel CPU.
In the backend, this affects which conda environment Snakemake uses, whenever hardware-accelerated environments are specified in a rule.

```yaml
os: intel
use_gpu: true
```

## 3. Input configuration

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

> ⚠️ **Warning** There can only be one instance of a module as a key in the input mapping (in the backend this is a dictionary). But you can reuse the same module output as input for multiple other modules. The order of the entries in the input mapping doesn't matter. 

## 4. Module configuration

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
