import traceback
import logging
logging.basicConfig(level=logging.INFO)
import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse
from tqdm import tqdm
try:
    import subprocess
    if subprocess.run('nvidia-smi', shell=True).returncode != 0:
        logging.info('No GPU found...')
        raise ImportError()
    import rapids_singlecell as sc
    import cupy as cp
    import rmm
    from rmm.allocators.cupy import rmm_cupy_allocator
    rmm.reinitialize(
        managed_memory=True,
        pool_allocator=False,
    )
    cp.cuda.set_allocator(rmm_cupy_allocator)
    logging.info('Using rapids_singlecell...')
    USE_GPU = True
except ImportError as e:
    import scanpy as sc
    logging.info('Importing rapids failed, using scanpy...')
    USE_GPU = False

from .assertions import assert_neighbors
from .io import to_memory, csr_matrix_int64_indptr
from .misc import dask_compute
from .accessors import _filter_genes


def compute_neighbors(adata, output_type=None, force=False, check_n_neighbors=False, **kwargs):
    """ Compute kNN graph based on output type.

    :param adata: integrated anndata object
    :param output_type: string of output type
    :param force: force re-computation of kNN graph
    :param kwargs: additional arguments for sc.pp.neighbors
    :return: anndata with kNN graph based on output type
    """
    
    if not output_type:
        output_type = 'full'

    try:
        logging.info(adata.__str__())
        assert not force, 'force neighbors computation'
        assert_neighbors(adata, check_n_neighbors=check_n_neighbors)
        logging.info(f'kNN graph already computed for {output_type}. Using pre-computed {output_type} kNN graph')
        return
    except AssertionError as e:
        logging.info(e.__str__())
        logging.info(f'Compute kNN graph for {output_type}...')

    if output_type == 'knn':
        assert_neighbors(adata, check_params=False)
        if 'params' not in adata.uns['neighbors']:
            adata.uns['neighbors']['params'] = adata.uns['neighbors'].get('params', {})
        adata.uns['neighbors']['params'] |= dict(use_rep='X')
    
    elif output_type == 'embed':
        assert 'X_emb' in adata.obsm, 'Embedding key "X_emb" not found'
        kwargs |= dict(use_rep='X_emb')
        sc.pp.neighbors(adata, **kwargs)
    
    elif output_type == 'full':
        if 'X_emb' not in adata.obsm:
            logging.info('Load corrected counts...')
            adata_tmp = adata
            if 'highly_variable' in adata_tmp.var.columns:
                adata_tmp = adata_tmp[:, adata_tmp.var['highly_variable']]
            adata_tmp = dask_compute(adata_tmp.copy(), layers='X')
            
            logging.info('Compute PCA on corrected counts...')
            sc.pp.pca(adata_tmp)
            adata.obsm['X_emb'] = adata_tmp.obsm['X_pca']
            del adata_tmp
        kwargs |= dict(use_rep='X_emb')
        logging.info('Neighbors...')
        sc.pp.neighbors(adata, **kwargs)
    
    else:
        raise ValueError(f'Invalid output type {output_type}')
    
    assert_neighbors(adata, check_n_neighbors=False)


def filter_genes(adata, batch_key=None, compute=True, **kwargs):
    """
    Filter anndata based on .X matrix
    """
    from dask import array as da
    import dask.config
    
    # filter cells from batches with too few cells
    if batch_key is not None:
        cells_per_batch = adata.obs[batch_key].value_counts()
        if cells_per_batch.min() < 2:
            adata = adata[adata.obs[batch_key].isin(cells_per_batch[cells_per_batch > 1].index)]
        
    # apply filter
    gene_mask = _filter_genes(adata, **kwargs)
    if any(gene_mask == False):
        logging.info(f'Subset to {sum(gene_mask)}/{adata.n_vars} filtered genes...')
        adata = adata[:, gene_mask]
    
    if adata.is_view:
        adata = adata.copy()    
    
    # convert dask array to csr matrix
    if compute and isinstance(adata.X, da.Array):
        with dask.config.set(**{'array.slicing.split_large_chunks': True}):
            adata = dask_compute(adata, layers='X')

    return adata


def get_pseudobulks(adata, group_key, agg='sum', dtype='float32'):
    from dask import array as da
    
    def aggregate(x, agg):
        if agg == 'sum':
            return x.sum(0)
        elif agg == 'mean':
            return x.mean(0)
        else:
            raise ValueError(f'invalid aggregation method "{agg}"')
    
    def _get_pseudobulk_matrix(adata, group_key, agg, dtype=dtype):
        X = adata.X
        value_counts = adata.obs[group_key].value_counts()
        
        # filter groups for at least 2 replicates
        value_counts = value_counts[value_counts >= 2]
        groups = value_counts.index
        
        print(f'Aggregate {len(groups)} pseudobulks...', flush=True)
        
        if isinstance(X, da.Array):
            from dask import config as dask_config
            from tqdm.dask import TqdmCallback
            
            df = adata.obs[[group_key]].reset_index(drop=True).query(f'{group_key} in @groups')
            df[group_key] = pd.Categorical(df[group_key], categories=groups, ordered=True)
            
            # # determine shuffling index
            # obs_chunks = X.chunks[0]
            # df['dask_chunk'] = np.repeat(np.arange(len(obs_chunks)), obs_chunks)
            # df['chunk_reorder'] = np.empty(df.shape[0], dtype=int)
            # df.loc[df[group_key].sort_values().index, 'chunk_reorder'] = np.arange(df.shape[0], dtype=int)
            # per_chunk_order = df.groupby('dask_chunk')['chunk_reorder'].apply(list).tolist()
            
            # # shuffle data by group_key
            # X = shuffle(X, indexer=per_chunk_order, axis=0)
            # X = X.rechunk((tuple(value_counts.values), -1))
            
            print(f'Sort and rechunk dask array by "{group_key}"...', flush=True)
            with dask_config.set(**{'array.slicing.split_large_chunks': False}):
                # sort dask array by group_key and rechunk by size of groups
                # the result should be a dask chunk per pseudobulk group
                X = X[df.sort_values(by=group_key).index.values] # also subsets for exclued groups
                X = X.rechunk((tuple(value_counts.values), -1))
                pseudobulks = X.map_blocks(lambda x: aggregate(x, agg), dtype=dtype)
            # else:
            #     import sparse
            #
            #     def sparse_indicator(groubpy):                
            #         n = len(groubpy)
            #         A = sparse.COO(
            #             coords=[groubpy.cat.codes, np.arange(n)],
            #             data=np.broadcast_to(1, n), # weight
            #             shape=(len(groubpy.cat.categories), n)
            #         )
            #         return da.from_array(A)

            #     print(f'Subset and convert dask array to Pydata sparse matrix...', flush=True)
            #     with dask_config.set(**{'array.slicing.split_large_chunks': False}):
            #         X = X.map_blocks(sparse.COO, dtype=dtype)
            #     indicator_matrix = sparse_indicator(groubpy=adata.obs[group_key])

            #     if agg == 'sum':
            #         pseudobulks = indicator_matrix @ X
            #     elif agg == 'mean':
            #         _sum = indicator_matrix @ X
            #         pseudobulks = _sum / _sum.sum(axis=1)
            #     else:
            #         raise ValueError(f'invalid aggregation method "{agg}"')
                
            #     pseudobulks = pseudobulks.map_blocks(lambda x: x.tocsr(), dtype=dtype)
            
            print(f'Compute aggregation...', flush=True)
            with TqdmCallback(
                bar_format='{percentage:3.0f}% |{bar}| {elapsed}<{remaining}\n',
                mininterval=10,
                delay=10,
            ):
                pseudobulks = pseudobulks.compute()
                        
        elif isinstance(X, (scipy.sparse.spmatrix, np.ndarray)):
            import scanpy as sc

            pbulk_adata = sc.get.aggregate(
                adata,
                by=group_key,
                func=[agg]
            )
            pbulk_adata = pbulk_adata[pbulk_adata.obs[group_key].isin(groups)].copy()
            pseudobulks = scipy.sparse.csr_matrix(pbulk_adata.layers[agg])
            groups = pbulk_adata.obs_names
            # miniters = max(10, len(value_counts) // 100)
            # pseudobulks = []
            # for group in tqdm(value_counts.index, desc='Aggregate groups', miniters=miniters):
            #     row_agg = aggregate(adata[adata.obs[group_key] == group].X, agg)
            #     row_agg = row_agg.A1 if isinstance(row_agg, np.matrix) else row_agg
            #     pseudobulks.append(row_agg)
            # pseudobulks = np.stack(pseudobulks, axis=0)
        else:
            raise ValueError(f'invalid type "{type(X)}"')
        
        return pseudobulks, groups

    pbulks, groups = _get_pseudobulk_matrix(adata, group_key=group_key, agg=agg)

    return ad.AnnData(
        pbulks,
        obs=aggregate_obs(adata, group_key, groups),
        var=adata.var,
        dtype='float32'
    )


def aggregate_obs(adata: ad.AnnData, group_key: str, groups: list):

    print(f'Aggregate metadata by {group_key}...', flush=True)

    obs = pd.DataFrame(groups, columns=[group_key]).merge(
        adata.obs.drop_duplicates(subset=[group_key]),
        on=group_key,
        how='left'
    )
    obs.index = obs[group_key].values

    # convert category columns to string
    for col in obs.columns:
        if obs[col].dtype.name == 'category':
            obs[col] = obs[col].astype(str).astype('category')

    return obs
