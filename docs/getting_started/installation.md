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

The different parts of the workflow (modules, rules) require specific conda environments.
The toolbox has 2 modes of dealing with conda evironments:

1. **EITHER:** For maximum reproducibility, but potentially more computational overhead over time, configure `env_mode: from_yaml` and let Snakemake create only the environments for the specific jobs you want to run. In this mode, **do NOT pre-install the environments locally** (option 2), if you want to avoid redundant copies of environments. Beware that any updates to the environment yaml files will trigger a new environment to be installed and you might need to clean up old environments regularly.

2. **OR:** Pre-install conda environments locally and set `env_mode: local` (default). This requires you to manually update conda environments when the YAML specifications change, but gives you full control on whether a new environment needs to be built and does not create any duplicate environments. In order to keep environment management overhead minimal, consider creating the environments you need for your specific workflow. Each module should have a section on the environments needed by it.

### Pre-installing environments
You will at least require the snakemake environment.

```
conda env create -f envs/snakemake.yaml
```

For all other environments:

1. If you are using `env_mode: from_yaml`, just add a [`--conda-create-envs-only`](https://snakemake.readthedocs.io/en/v7.31.1/executing/cli.html#conda) to your snakemake command (see below for example workflow configuration and runner script). Do NOT use `envs/install_all_environments.sh` in this mode. Environments will be saved under `.snakemake/conda` from wherever you call the snakemake commands, so make sure that directory has sufficient space or create a symlink for `.snakemake/conda` to a different location (e.g. scratch) or refer to the [Snakemake documentation](https://snakemake.readthedocs.io/en/v7.31.1/snakefiles/deployment.html#integrated-package-management) on configuring the location of the resulting conda environments.

2. If you are using `env_mode: local` and want to pre-install all environments, you can use `bash envs/install_all_environments.sh` (check out the help message with `-h` first to understand the input specifications). If you want to only install specific environments, use `conda env create -f envs/<env_name>.yaml` (or `bash install_environment -f envs/<env_name>.yaml` to automatically install or update the environment).

> 📝 **Notes on `install_all_environments.sh`**
> 1. The script will create new environments for each file in the `envs` directory if they don't yet exist and update any pre-existing environments.
> 2. The environment names correspond the their respective file names and are documented under the `name:` directive in the `envs/<env_name>.yaml` file.
> 3. If an environment creation fails, the script will skip that environment and you might need to troubleshoot the installation manually.
> 4. Some environments require the channel priority to be set to `flexible`.
> If your installation command fails, try setting `conda config --set channel_priority flexible` before restarting the command.

### Removing environments

In cases where the environments have updated, you might want to have a clean install of the new environment.

1. For `env_mode: from_yaml`, any updates to the environment specification will trigger creating a new environment under `.snakemake/conda/envs`, but the old environment persists. You might need to clean up old environments once in a while, which is possbile via [`--conda-cleanup-envs`](https://snakemake.readthedocs.io/en/v7.31.1/executing/cli.html#conda). Read the [Snakemake documentation](https://snakemake.readthedocs.io/en/v7.31.1/snakefiles/deployment.html#integrated-package-management) on more information on package management, which includes pre-building environments or removing old environments.

2. For `env_mode: local`, you can manually remove conda enviroments via `conda env remove -n <env_name>`. The `install_all_environments.sh` scripts also offers an option to remove all environments that are defined under `envs/*.yaml` by providing the `-r` options. Consider running `bash install_all_environments -r -n` for a dry run first to avoid surprises. Afterwards, recreate your environments as needed.
