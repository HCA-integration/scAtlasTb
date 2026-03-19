.. _advanced-configuration:

⚙️ Advanced configuration
=========================


Set defaults
------------

You can set module-specific defaults that will be used for all tasks (under ``configs['DATASETS']``), if the parameters have not been specified for those tasks.
This can shorten the configuration file, make it more readable and help avoid misconfiguration if you want to reuse the same configurations for multiple tasks.

Under the ``defaults`` directive, you can set the defaults in the same way as the task-specific configuration.

.. dropdown:: Example defaults for modules
   :icon: beaker

   .. code-block:: yaml

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

Additionally to the module defaults, you can set which datasets you want to include in your workflow, without having to remove or comment out any entries in ``configs['DATASETS']``.

.. code-block:: yaml

   defaults:
   ...
     datasets:
     # list of dataset/task names that you want your workflow to be restricted to
       - test
       - test2

Automatic environment management
--------------------------------

Snakemake supports automatically creating conda environments for each rule.

.. code-block:: yaml

   env_mode: from_yaml

You can trigger Snakemake to install all environments required for your workflow in advance by adding the following parameters:

.. code-block:: bash

   <snakemake_cmd> --use-conda --conda-create-envs-only --cores 1

.. _snakemake_profiles:

Snakemake profiles
------------------

Snakemake profiles help you manage the many flags and options of a snakemake command in a single file, which will simplify the Snakemake call considerably.
The toolbox provides some example Snakemake profiles under ``.profiles``, which you can copy and adapt to your needs.

To use a profile (e.g. the local profile), add ``--profile .profiles/<profile_name>`` to your Snakemake command.
You can read more about profiles in `Snakemake's documentation <https://snakemake.readthedocs.io/en/v7.31.1/executing/cli.html#profiles>`_.


.. _cluster_execution:

Cluster execution
-----------------

Snakemake supports scheduling rules as jobs on a cluster and ``scAtlasTb`` has been developed and tested to work with SLURM.
In order for Snakemakek  to deploy jobs through SLURM create a Snakemake profile under ``.profiles/<your_profile>/config.yaml``.

.. dropdown:: Example profile for SLURM
   :icon: beaker

   Adapted from https://github.com/jdblischak/smk-simple-slurm

   .. code-block:: yaml

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


You can find detailed information on cluster execution in the `Snakemake documentation <https://snakemake.readthedocs.io/en/v7.31.1/executing/cluster.html>`_.

Additionally, you need to configure ``cpu`` and ``gpu`` resource settings in your workflow config file (NOT the snakemake profile config), under ``resources``.
These settings will be used by scAtlasTb to determine how to schedule the different rules on the cluster, depending on whether they require GPU or not.
For each resource profile, you need to set the memory requirements (``mem_mb``), ``partition``, ``qos``, and ``gpu`` count.
If your system does define any ``qos``, you can set it to ``normal``, which is the default.

.. code-block:: yaml

  use_gpu: true  # assuming you have and want to use GPU nodes
  resources:
    cpu: # set CPU resource settings for rules that do not require GPU; must be called "cpu"
      partition: cpu
      qos: normal
      gpu: 0
      mem_mb: 100000
    gpu: # set GPU resource settings for rules that require GPU; must be called "gpu"
      partition: gpu
      qos: normal
      gpu: 1
      mem_mb: 100000

In order for jobs to make use of the GPU, make sure that your GPU environments are installed correctly (see `Working with GPUs <troubleshooting.html#working-with-gpus>`_).
If you don't have GPU nodes, you can configure the gpu resources to be the same as the cpu resources.

.. code-block:: yaml

  use_gpu: false
  resources:
    cpu:
      partition: cpu
      qos: normal
      gpu: 0
      mem_mb: 100000
    gpu:
      partition: cpu
      qos: normal
      gpu: 1
      mem_mb: 100000

Finally, include the snakemake profile in your Snakemake command:

.. code-block:: bash

   <snakemake_cmd> --profile .profiles/<your_profile>

with ``<snakemake_cmd>`` being either ``snakemake`` or  your runner script (recommended, e.g. ``run.sh``).