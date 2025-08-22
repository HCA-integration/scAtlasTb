# ⚙️ Advanced configuration
<a name="advanced-configuration"></a>

## Set defaults

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


## Automatic environment management
Snakemake supports automatically creating conda environments for each rule.

```yaml
env_mode: from_yaml
```

You can trigger Snakemake to install all environments required for your workflow in advance by adding the following parameters

```
<snakemake_cmd> --use-conda --conda-create-envs-only --cores 1
```

## Snakemake profiles

Snakemake profiles help you manage the many flags and options of a snakemake command in a single file, which will simplify the Snakemake call considerably.
The toolbox provides some example Snakemake profiles are provided under `.profiles`, which you can copy and adapt to your needs.

To use a profile e.g. the local profile, add `--profile .profiles/<profile_name>` to your Snakemake command.
You can read more about profiles in [Snakemake's documentation](https://snakemake.readthedocs.io/en/v7.31.1/executing/cli.html#profiles).

## Cluster execution

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
