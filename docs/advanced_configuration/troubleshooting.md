# 🛠️ Troubleshooting
<a name="troubleshooting"></a>

## Conda environment installation fails
For some environments (e.g. `envs/rapids_singlecell.yaml`) you need to set the channel priority to flexible in order for conda to properly resolve the environment.

```
conda config --set channel_priority flexible
```

After setting the flag and calling the installation command, the environment should resolve.

## Working with GPUs

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

## Working with CPUs only

If your system doesn't have any GPUs, you can set the following flag in your config.

```yaml
use_gpu: false
```

This will force Snakemake to use the CPU versions of an environment.

## FAQs

Below are some scenarios that can occur when starting with the pipeline.
If you have any additional questions or encounter any bugs, please open up a [github issue](https://github.com/HCA-integration/scAtlasTb/issues).
If you want to contribute to improving the pipeline, check out the [contribution guidelines](CONTRIBUTING.md).

### I configured my pipeline and the dry run doesn't fail, but it doesn't want to run the modules I configured. What do I do?

This likely happens when you don't specify which rule you want Snakemake to run. By default, Snakemake will try create a visualisation of the modules you configured. If you want it to run the modules themselves, you will need to add the rule name with your Snakemake command. For each rule, there is a `<module>_all`, but you can view all possible rules through `snakemake -l`

