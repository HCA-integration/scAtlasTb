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


def check_and_set_neighbors_key(adata, neighbors_key):
    neighbors = adata.uns.get(neighbors_key, {})
    neighbors |= {
        'connectivities_key': neighbors.get('connectivities_key', 'connectivities'),
        'distances_key': neighbors.get('distances_key', 'distances'),
        'params': neighbors.get('params', {'use_rep': 'X_pca', 'method': None}),
    }
    adata.uns[neighbors_key] = neighbors
    
    adata.obsp['connectivities'] = adata.obsp[neighbors['connectivities_key']]
    adata.obsp['distances'] = adata.obsp[neighbors['distances_key']]


def apply_clustering(
    adata,
    resolution: float,
    cpu_kwargs: dict = None,
    n_cell_cpu: int = 300_000,
    max_cluster_factor: int = 50,
    recompute_neighbors: bool = False,
    neighbors_args: dict = {},
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
    
    if not cpu_kwargs:
        cpu_kwargs = dict()
    
    if not USE_GPU:
        kwargs |= cpu_kwargs
    
    if recompute_neighbors:
        neighbors(adata, **neighbors_args)
    
    if USE_GPU:
        # following observations from https://github.com/rapidsai/cugraph/issues/4072#issuecomment-2074822898
        adata.obsp['connectivities'] = adata.obsp['connectivities'].astype('float64')

    if adata.n_obs < n_cell_cpu:
        # switch to CPU implementation for smaller numbers of cells
        algorithm_map = alt_algorithm_map
        
    # logging.info(f'{algorithm} clustering with {kwargs} for {adata.n_obs} cells...')
    cluster_func = algorithm_map.get(algorithm, KeyError(f'Unknown clustering algorithm: {algorithm}'))
    cluster_func(adata, **kwargs)
    
    # heuristic to check if number of clusters is reasonable
    max_clusters = max(1, int(max_cluster_factor * resolution))
    n_clusters = adata.obs[cluster_key].nunique()
    
    if USE_GPU and n_clusters > max_clusters:
        # fallback when too many clusters are computed (assuming this is a bug in the rapids implementation)
        logging.info(
            f'Cluster {cluster_key} has {n_clusters} custers, which is more than {max_clusters}. '
            'Falling back to scanpy implementation...'
        )
        cluster_func = alt_algorithm_map[algorithm]
        kwargs |= cpu_kwargs
        cluster_func(adata, **kwargs)
    return adata


input_file = snakemake.input[0]
output_file = snakemake.output[0]
resolution = float(snakemake.wildcards.resolution)
algorithm = snakemake.wildcards.algorithm
level = int(snakemake.wildcards.level)
threads = snakemake.threads
overwrite = snakemake.params.get('overwrite', False)
max_cluster_factor = snakemake.params.get('max_cluster_factor', 50)

# set parameters for clustering
cluster_key = f'{algorithm}_{resolution}_{level}'
kwargs = dict(
    resolution=resolution,
    key_added=cluster_key,
) | snakemake.params.get('clustering_args', {})
cpu_kwargs = dict(flavor='igraph')
if algorithm == 'leiden':
    cpu_kwargs |= dict(n_iterations=2)

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
        neighbors_key = snakemake.params.get('neighbors_key', 'neighbors')
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
            **kwargs,
        )
    else:
        from joblib import Parallel, delayed
        from tqdm import tqdm
        
        def cluster_subset(
            adata,
            resolution,
            key_added,
            prev_cluster_key,
            prev_cluster_value,
            neighbors_args, # TODO: custom parameters
        ):
            # logging.info(f'Subsetting to {prev_cluster_key}={prev_cluster_value}...')
            sub_adata = adata[adata.obs[prev_cluster_key] == prev_cluster_value].copy()
            
            if sub_adata.n_obs < 2 * neighbors_args.get('n_neighbors', 15):
                prev_cluster_clean = sub_adata.obs[prev_cluster_key].astype(str).apply(lambda x: x.split('_')[-1])
                return sub_adata.obs[prev_cluster_key].str.cat(prev_cluster_clean, sep='_')
            
            sub_adata = apply_clustering(
                sub_adata,
                resolution=resolution,
                key_added=key_added,
                cpu_kwargs=cpu_kwargs,
                max_cluster_factor=max_cluster_factor,
                neighbors_args=neighbors_args,
                recompute_neighbors=True,
            )
            return sub_adata.obs[[prev_cluster_key, key_added]].agg('_'.join, axis=1)
        
        neighbors_args = snakemake.params.get('neighbors_args', {})
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
        
        adata = read_anndata(input_file, **read_kwargs)
        adata.obsm[use_rep] = read_slot(input_file, store, f'obsm/{use_rep}')
        
        prev_cluster_key = f'{algorithm}_{resolution}_{level-1}'
        cluster_labels = adata.obs[prev_cluster_key].unique()
        
        logging.info(f'Will recompute neighbors with {neighbors_args}')
        logging.info(f'{algorithm} clustering with {kwargs} for clusters from {cluster_key}...')
            
        results = Parallel(n_jobs=threads)(
            delayed(cluster_subset)(
                adata,
                prev_cluster_key=prev_cluster_key,
                prev_cluster_value=prev_cluster,
                neighbors_args=neighbors_args,
                **kwargs
            )
            for prev_cluster in tqdm(cluster_labels, desc=f'Cluster with {threads} threads', miniters=1)
        )

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