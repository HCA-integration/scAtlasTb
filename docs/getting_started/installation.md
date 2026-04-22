# 📦 Installation

## Clone the repository

Depending on whether you have set up SSH or HTTPS with [PAT](https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/managing-your-personal-access-tokens), you can clone the repository.

SSH:
```
git clone git@github.com:HCA-integration/scAtlasTb.git
```

HTTPS:
```
git clone https://github.com/HCA-integration/scAtlasTb.git
```

## Requirements

* Linux (preferred) or MacOS (not rigorously tested, some bioconda dependencies might not work out-of-the-box)
* Conda, e.g., via [miniforge](https://github.com/conda-forge/miniforge) (recommended) or [miniconda](https://docs.anaconda.com/free/miniconda/index.html)


The modules are tested and developed using task-specific conda environments, which should be quick to set up when using [libmamba](https://www.anaconda.com/blog/a-faster-conda-for-a-growing-community).

> 📝  **Note** If you use conda version 22.11 or above, make sure you set the conda solver to [`libmamba`](https://www.anaconda.com/blog/a-faster-conda-for-a-growing-community) for significantly faster installation. For newer versions or if you are using mamba directly, `libmamba` should already be the default.

## Install dependencies

All the conda environments used by the toolbox are under `envs/*.yaml`.
You will at least require the snakemake environment.

```
conda env create -f envs/snakemake.yaml
```

You can install the other environments as needed, for different parts of the workflow (modules, rules).