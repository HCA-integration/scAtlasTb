import warnings
warnings.filterwarnings("ignore")
import logging
logging.basicConfig(level=logging.INFO)
import numpy as np
import scanpy as sc
USE_GPU = False
try:
    import subprocess
    assert subprocess.run('nvidia-smi', shell=True, stdout=subprocess.DEVNULL).returncode == 0
    from rapids_singlecell.tl import leiden, louvain
    from rapids_singlecell.pp import neighbors
    USE_GPU = True
except Exception as e:
    logging.info(f'Importing rapids failed, using scanpy implementation\n{e}')
    from scanpy.tools import leiden, louvain
    from scanpy.preprocessing import neighbors

from utils.io import read_anndata, write_zarr_linked, get_store, read_slot


def apply_clustering(
    adata,
    resolution: float,
    cpu_kwargs: dict = None,
    n_cell_cpu: int = 100_000,
    max_cluster_factor: int = 100,
    recompute_neighbors: bool = False,
    neighbors_args: dict = {},
    use_gpu: bool = USE_GPU,
    **kwargs,
):
    """
    :param adata: anndata object
    :param cpu_kwargs: clustering parameters for CPU implementation
    :param n_cell_cpu: number of cells to use CPU implementation
    :param max_cluster_factor: factor to calculate maximum number of clusters to determine if number of clusters are correct
    """
    algorithm_map = {
        'louvain': louvain,
        'leiden': leiden,
    }
    alt_algorithm_map = {
        'louvain': sc.tl.louvain,
        'leiden': sc.tl.leiden,
    }
    
    kwargs |= dict(resolution=resolution)
    key_added = kwargs.get('key_added', algorithm)
    
    if not cpu_kwargs:
        cpu_kwargs = dict()
    
    if not use_gpu:
        kwargs |= cpu_kwargs
    
    if recompute_neighbors:
        neighbors(adata, **neighbors_args)
    
    if use_gpu:
        # following observations from https://github.com/rapidsai/cugraph/issues/4072#issuecomment-2074822898
        adata.obsp['connectivities'] = adata.obsp['connectivities'].astype('float64')

    if adata.n_obs < n_cell_cpu:
        # switch to CPU implementation for smaller numbers of cells
        algorithm_map = alt_algorithm_map
        
    # logging.info(f'{algorithm} clustering with {kwargs} for {adata.n_obs} cells...')
    cluster_func = algorithm_map.get(algorithm, KeyError(f'Unknown clustering algorithm: {algorithm}'))
    cluster_func(adata, **kwargs)
    
    # heuristic to check if number of clusters and sizes are reasonable
    max_clusters = max(1, int(max_cluster_factor * resolution))
    n_clusters = adata.obs[key_added].nunique()
    # check if smallest cluster size is reasonable
    min_cluster_size = adata.obs[key_added].value_counts().min()
    quantile_cluster_size = adata.obs[key_added].value_counts().quantile(0.1)
    size_ok = (quantile_cluster_size > adata.n_obs / max_clusters) or (min_cluster_size > 10)
    
    if use_gpu and n_clusters > max_clusters and not size_ok:
        # fallback when too many clusters are computed (assuming this is a bug in the rapids implementation)
        logging.info(
            f'Cluster {key_added} has {n_clusters} clusters, which is more than {max_clusters}. '
            'Falling back to scanpy implementation...'
        )
        cluster_func = alt_algorithm_map[algorithm]
        kwargs |= cpu_kwargs
        cluster_func(adata, **kwargs)
    return adata


def cluster_subset(
    adata,
    resolution,
    key_added,
    prev_cluster_key,
    prev_cluster_value,
    neighbors_args, # TODO: custom parameters
    use_gpu,
    n_cell_cpu,
    cpu_kwargs,
    max_cluster_factor,
):
    """
    Wrapper for calling apply_clustering function in parallel
    """
    import warnings
    warnings.filterwarnings("ignore")
    
    adata = adata[adata.obs[prev_cluster_key] == prev_cluster_value].copy()
    
    if adata.n_obs < 2 * neighbors_args.get('n_neighbors', 15):
        prev_cluster_clean = adata.obs[prev_cluster_key].astype(str).apply(lambda x: x.split('_')[-1])
        return adata.obs[prev_cluster_key].str.cat(prev_cluster_clean, sep='_')
    
    adata = apply_clustering(
        adata,
        resolution=resolution,
        key_added=key_added,
        cpu_kwargs=cpu_kwargs,
        max_cluster_factor=max_cluster_factor,
        neighbors_args=neighbors_args,
        recompute_neighbors=True,
        use_gpu=use_gpu,
        n_cell_cpu=n_cell_cpu,
    )
    
    return adata.obs[[prev_cluster_key, key_added]].agg('_'.join, axis=1)


if 'snakemake' in globals():
    input_file = snakemake.input[0]
    output_file = snakemake.output[0]
    resolution = float(snakemake.wildcards.resolution)
    algorithm = snakemake.wildcards.algorithm
    level = int(snakemake.wildcards.level)
    threads = snakemake.threads
    overwrite = snakemake.params.get('overwrite', False)
    max_cluster_factor = snakemake.params.get('max_cluster_factor', 50)
    clustering_args = snakemake.params.get('clustering_args', {})
    neighbors_key = snakemake.params.get('neighbors_key', 'neighbors')
    neighbors_args = snakemake.params.get('neighbors_args', {})
    n_cell_cpu = snakemake.params.get('n_cell_cpu', 100_000)
    
else:
    import argparse
    import json
    
    parser = argparse.ArgumentParser(description='Run clustering on anndata file.')
    parser.add_argument('input_file', type=str, help='Input anndata file')
    parser.add_argument('output_file', type=str, help='Output anndata file')
    parser.add_argument('--resolution', type=float, default=1.0, help='Clustering resolution')
    parser.add_argument('--algorithm', type=str, choices=['louvain', 'leiden'], default='leiden', help='Clustering algorithm')
    parser.add_argument('--level', type=int, default=1, help='Hierarchical clustering level')
    parser.add_argument('--threads', type=int, default=1, help='Number of threads to use')
    parser.add_argument('--overwrite', action='store_true', help='Overwrite existing clustering results')
    parser.add_argument('--max_cluster_factor', type=int, default=50, help='Maximum cluster factor for heuristic check (GPU only)')
    parser.add_argument('--clustering_args', type=json.loads, default={}, help='Additional clustering arguments')
    parser.add_argument('--neighbors_key', type=str, default='neighbors', help='Key for neighbors in adata.uns')
    parser.add_argument('--neighbors_args', type=json.loads, default={}, help='Additional arguments for neighbors computation')
    parser.add_argument('--n_cell_cpu', type=int, default=100_000, help='Number of cells for which to force CPU computation')
    args = parser.parse_args()
    
    input_file = args.input_file
    output_file = args.output_file
    resolution = args.resolution
    algorithm = args.algorithm
    level = args.level
    threads = args.threads
    overwrite = args.overwrite
    max_cluster_factor = args.max_cluster_factor
    clustering_args = args.clustering_args
    neighbors_key = args.neighbors_key
    neighbors_args = args.neighbors_args
    n_cell_cpu = args.n_cell_cpu

# set parameters for clustering
cluster_key = f'{algorithm}_{resolution}_{level}'
kwargs = dict(resolution=resolution, key_added=cluster_key) | clustering_args

if 'flavor' in kwargs:
    # if user defines flavor, force clustering on CPU
    USE_GPU = False

# use user kwargs for cpu_kwargs, if provided
cpu_kwargs = dict(flavor=kwargs.pop('flavor', 'igraph'))
if algorithm == 'leiden':
    cpu_kwargs |= dict(n_iterations=kwargs.pop('n_iterations', 2))

logging.info(f'Using GPU: {USE_GPU}')

# check if clusters have already been computed
store = get_store(input_file)
cluster_key_exists = cluster_key in store['obs'].keys()


if cluster_key_exists and not overwrite:
    logging.info(f'Read anndata file {input_file} and skip clustering...')
    adata = read_anndata(input_file, obs='obs')
else:
    read_kwargs = dict(obs='obs')
    
    if level <= 1:        
        logging.info(f'Read anndata file {input_file}...')
        if input_file.endswith('.h5ad'):
            read_kwargs['obsm'] = 'obsm'  # preserve obsm for later hierarchical clustering
        adata = read_anndata(input_file, uns='uns', **read_kwargs)
        
        # select neighbors
        neighbors = adata.uns.get(neighbors_key, {})
        conn_key = neighbors.get('connectivities_key', 'connectivities')
        dist_key = neighbors.get('distances_key', 'distances')
        adata.obsp = {
            'connectivities': read_slot(input_file, store, f'obsp/{conn_key}'),
            'distances': read_slot(input_file, store, f'obsp/{dist_key}'),
        }
        
        logging.info(f'{algorithm} clustering with {kwargs} for {adata.n_obs} cells...')
        adata = apply_clustering(
            adata,
            cpu_kwargs=cpu_kwargs,
            recompute_neighbors=False,
            use_gpu=USE_GPU,
            max_cluster_factor=max_cluster_factor,
            n_cell_cpu=n_cell_cpu,
            **kwargs,
        )
    else:
        from joblib import Parallel, delayed
        from tqdm import tqdm
        import gc
        
        use_rep = neighbors_args.get('use_rep', 'X_pca')
        neighbors_args['use_rep'] = use_rep
        
        # check if use_rep is present in obsm
        if use_rep not in store['obsm'].keys():
            params = read_slot(input_file, store, 'uns/neighbors/params')
            if 'use_rep' not in params:
                raise ValueError(f'use_rep not defined in neighbors_args for {params}, consider recomputing neighbors')
            use_rep = params['use_rep']
            neighbors_args['use_rep'] = use_rep
            assert use_rep in store['obsm'].keys(), f'obsm key {use_rep} not found in {input_file}'
        
        # read data
        adata = read_anndata(input_file, **read_kwargs)
        adata.obsm[use_rep] = read_slot(input_file, store, f'obsm/{use_rep}')
        
        prev_cluster_key = f'{algorithm}_{resolution}_{level-1}'
        cluster_labels = adata.obs[prev_cluster_key].unique()
        
        logging.info(f'Will recompute neighbors with {neighbors_args}')
        logging.info(f'{algorithm} clustering with {kwargs} for clusters from {cluster_key}...')
        
        results = list(
            tqdm(
                Parallel(return_as='generator', n_jobs=threads)(
                    delayed(cluster_subset)(
                        adata,
                        prev_cluster_key=prev_cluster_key,
                        prev_cluster_value=prev_cluster,
                        neighbors_args=neighbors_args,
                        use_gpu=USE_GPU,
                        max_cluster_factor=max_cluster_factor,
                        n_cell_cpu=n_cell_cpu,
                        cpu_kwargs=cpu_kwargs,
                        **kwargs
                    ) for prev_cluster in cluster_labels
                ),
                desc=f'Cluster with {threads} threads',
                total=len(cluster_labels),
            )
        )
        gc.collect()

        for clusters in results:
            adata.obs.loc[clusters.index, cluster_key] = clusters


logging.info(f'Write {cluster_key} to {output_file}...')
adata.obs = adata.obs[[cluster_key]]
write_zarr_linked(
    adata,
    in_dir=input_file,
    out_dir=output_file,
    files_to_keep=['obs'],
)
