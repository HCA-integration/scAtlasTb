# Installation

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

* Linux (preferred) or MacOS on Intel (not rigorously tested, some bioconda dependencies might not work)
* conda e.g. via [miniforge](https://github.com/conda-forge/miniforge)(recommended) or [miniconda](https://docs.anaconda.com/free/miniconda/index.html)


The modules are tested and developed using task-specific conda environments, which should be quick to set up when using [mamba](https://mamba.readthedocs.io).

> 📝  **Note** If you use conda, but have never used mamba, consider installing the mamba package into your base environment and use it for all installation commands.
You can still replace all mamba commands with conda commands if you don't want to install mamba.

## Install conda environments

The different parts of the workflow (modules, rules) require specific conda environments.
The simplest way to install all environments is to run the following script:

```
bash envs/install_all_environments.sh -h # help message for customization
bash envs/install_all_environments.sh
```

> 📝 **Notes**
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
