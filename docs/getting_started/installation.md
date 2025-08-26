# 📦 Installation

## Clone the repository

Depending on whether you have set up SSH or HTTPS with [PAT](https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/managing-your-personal-access-tokens), you can clone the repository with

SSH:
```
git clone git@github.com:HCA-integration/scAtlasTb.git
```

HTTPS:
``` clone
git clone https://github.com/HCA-integration/scAtlasTb.git
```

## Requirements

* Linux (preferred) or MacOS (not rigorously tested, some bioconda dependencies might not work out-of-the-box)
* Conda e.g. via [miniforge](https://github.com/conda-forge/miniforge)(recommended) or [miniconda](https://docs.anaconda.com/free/miniconda/index.html)


The modules are tested and developed using task-specific conda environments, which should be quick to set up when using [libmamba](https://www.anaconda.com/blog/a-faster-conda-for-a-growing-community).

> 📝  **Note** If you use conda version 22.11 or above, make sure you set the conda solver to [`libmamba`](https://www.anaconda.com/blog/a-faster-conda-for-a-growing-community) for significantly faster installation. For newer versions or if you are using mamba directly, `libmamba` should already be the default.

## Install Conda environments

All the conda environments used by the toolbox are under `envs/*.yaml`.
You will at least require the snakemake environment.

```
conda env create -f envs/snakemake.yaml
```

The different parts of the workflow (modules, rules) require specific conda environments.
The toolbox has 2 modes of dealing with conda evironments:

1. **EITHER:** For maximum reproducibility, but potentially more computational overhead over time, configure `env_mode: from_yaml` and let Snakemake create only the environments for the specific jobs you want to run. In this mode, **do NOT pre-install the environments locally** (option 2), if you want to avoid redundant copies of environments. Beware that any updates to the environment yaml files will trigger a new environment to be installed and you might need to clean up old environments regularly.

2. **OR:** Pre-install conda environments locally and set `env_mode: local` (default). This requires you to manually update conda environments when the YAML specifications change, but gives you full control on whether a new environment needs to be built and does not create any duplicate environments. In order to keep environment management overhead minimal, consider creating the environments you need for your specific workflow. Each module should have a section on the environments needed by it.

## Option 1: Managing environments for `env_mode: from_yaml`

Set the global parameter in your configuration file

```yaml
env_mode: from_yaml
```

> **Note**
> Do NOT use `envs/install_all_environments.sh` in this mode!


### Creating environments
If you are using `env_mode: from_yaml`, just add a [`--conda-create-envs-only`](https://snakemake.readthedocs.io/en/v7.31.1/executing/cli.html#conda) to your snakemake command (see below for example workflow configuration and runner script).  Environments will be saved under `.snakemake/conda` from wherever you call the snakemake commands, so make sure that directory has sufficient space or create a symlink for `.snakemake/conda` to a different location (e.g. scratch) or refer to the [Snakemake documentation](https://snakemake.readthedocs.io/en/v7.31.1/snakefiles/deployment.html#integrated-package-management) on configuring the location of the resulting conda environments.

```
snakemake <target_rule> --conda-create-envs-only
```

### Updating environments

Any updates to the environment specification will trigger creating a new environment under `.snakemake/conda/envs`, but the old environment persists. You might need to clean up old environments once in a while, which is possbile via [`--conda-cleanup-envs`](https://snakemake.readthedocs.io/en/v7.31.1/executing/cli.html#conda).

```
snakemake <target_rule> --conda-cleanup-envs
```

Read the [Snakemake documentation](https://snakemake.readthedocs.io/en/v7.31.1/snakefiles/deployment.html#integrated-package-management) on more information on package management, which includes pre-building environments or removing old environments.


## Option 2: Managing environments for `env_mode: local`

Optional: set the global parameter in your configuration file. This is the default, so it will be used, even when `env_mode` is not configured.

```yaml
env_mode: local
```

### Installing environments

If you want to only install specific environments (which is recommended if you're starting out with small workflows), you can install the enviroment with conda:

```
conda env create -f envs/<env_name>.yaml
```

Alternatively, you can use the `install_environment.sh`, which automatically install or update the environment

```
bash install_environment -h  # help message
bash install_environment -f envs/<env_name>.yaml
```

If, instead, you want to pre-install all environments, `envs/install_all_environments.sh` provides a convenient wrapper.

```
bash envs/install_all_environments.sh -h  # help message
bash envs/install_all_environments.sh -n  # dry run
bash envs/install_all_environments.sh
```

> 📝 **Notes on `install_all_environments.sh`**
> 1. The script will create new environments for each file in the `envs` directory if they don't yet exist and update any pre-existing environments.
> 2. The environment names correspond the their respective file names and are documented under the `name:` directive in the `envs/<env_name>.yaml` file.
> 3. If an environment creation fails, the script will skip that environment and you might need to troubleshoot the installation manually.
> 4. Some environments require the channel priority to be set to `flexible`.
> If your installation command fails, try setting `conda config --set channel_priority flexible` before restarting the command.

### Updating environments

In cases where the environments have updated, you might want to have a clean install of the new environment, instead of updating existing environments.

You can manually remove conda enviroments via

```
conda env remove -n <env_name>
```

If you want to remove all toolbox-related environments, you can use

```
install_all_environments.sh -r -n  # dry run (recommended)
install_all_environments.sh -r
```

This will remove all environments that are defined under `envs/*.yaml`.
After removing all environments, recreate your environments as needed.
