# Clustering

This module provides a workflow to perform clustering at multiple resolutions and at multiple subclustering levels (iterative sub-clustering).
The code is optimised for memory efficiency and speed, allowing for clustering of large datasets.

## Configuration

Example configuration for clustering can be found in `workflow/clustering/config.yaml` showcasing the setup for hierarchical (iterative) sub-clustering.

```yaml
DATASETS:
  test_hierarchical:
    input:
      clustering: test/input/pbmc68k.h5ad
    clustering:
      neighbors_key: neighbors
      algorithm: leiden
      resolutions:
        - 0.5
      hierarchy:
        1: 0.1
        3: 0.2
      umap_colors:
        - bulk_labels
        - batch
        - n_genes
```

## Calling the clustering script from the command line

Optionally, you can call the clustering script from the command line.
The script is located in the `workflow/clustering/scripts/clustering.py`.
It assumes that you have already computed neighbors on your data, which you can specify with the `--neighbors_key` argument.

Make sure you are using the `envs/scanpy.yaml` environment, which has all the necessary dependencies installed.

Check the help message for the script to see all available options:

```bash
python workflow/scripts/clustering.py -h
```

```
usage: clustering.py [-h] [--resolution RESOLUTION] [--algorithm {louvain,leiden}] [--level LEVEL] [--threads THREADS] [--overwrite] [--max_cluster_factor MAX_CLUSTER_FACTOR] [--clustering_args CLUSTERING_ARGS] [--neighbors_key NEIGHBORS_KEY] [--neighbors_args NEIGHBORS_ARGS]
                     input_file output_file

Run clustering on anndata file.

positional arguments:
  input_file            Input anndata file
  output_file           Output anndata file

options:
  -h, --help            show this help message and exit
  --resolution RESOLUTION
                        Clustering resolution
  --algorithm {louvain,leiden}
                        Clustering algorithm
  --level LEVEL         Hierarchical clustering level
  --threads THREADS     Number of threads to use
  --overwrite           Overwrite existing clustering results
  --max_cluster_factor MAX_CLUSTER_FACTOR
                        Maximum cluster factor for heuristic check (GPU only)
  --clustering_args CLUSTERING_ARGS
                        Additional clustering arguments
  --neighbors_key NEIGHBORS_KEY
                        Key for neighbors in adata.uns
  --neighbors_args NEIGHBORS_ARGS
                        Additional arguments for neighbors computation
```

In the following is an example on how to run the clustering for 2 different levels:

```bash
# run level 1 with resolution 1.0
python workflow/scripts/clustering.py \
    <input zarr or h5ad> \
    output_level=1.zarr \
    --resolution 1.0 \
    --level 1 \

# run level 2 with resolution 0.2 on level 1 output
python workflow/scripts/clustering.py \
    output_level=1.zarr \
    output_level=2.zarr \
    --resolution 1.0 \
    --level 2 \
    --clustering_args '{"resolution": 0.2}' \
    --threads 5
```

`output_level=1.zarr` will contain an `.obs` column called `leiden_1.0_1` with the clustering results at level 1 with resolution 1.0, while `output_level=2.zarr` will contain an `.obs` column called `leiden_1.0_2` with the clustering results at level 2 with resolution 0.2.
Note, that the resolution in the `.obs` column is always the resolution of level 1.
You can of course adjust the output path as you like.

Each call will only output 1 clustering result, so you need to collect the results for all clustering resolutions.

```python
import pandas as pd
from utils.io import read_anndata

obs = pd.concat(
    [
        read_anndata(file, obs='obs', verbose=False).obs 
        for file in ['output_level=1.zarr', 'output_level=2.zarr']
    ], axis=1
)
```