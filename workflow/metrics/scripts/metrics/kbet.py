import numpy as np
from .utils import rename_categories, select_neighbors, scanpy_to_neighborsresults


def kbet_y(adata, output_type, batch_key, label_key, **kwargs):
    import scib_metrics
    
    adata = select_neighbors(adata, output_type)
    batches = rename_categories(adata, batch_key)

    score, _ , _ = scib_metrics.kbet(
        X=scanpy_to_neighborsresults(adata),
        batches=batches,
    )

    return score


def kbet_pg(adata, output_type, batch_key, label_key, n_threads=1, **kwargs):
    from tqdm import tqdm
    import pegasusio
    from scipy.sparse import csr_matrix
    from joblib import Parallel, delayed
    
    
    def compute_kbet(mmdata, *args, **kwargs):
        import pegasus as pg    
        stat_mean, pvalue_mean, accept_rate = pg.calc_kBET(mmdata, *args, **kwargs)
        return accept_rate
    
    
    if output_type == 'knn':
        return np.nan
    
    # prepare adata for pegasusio
    adata.X = csr_matrix(adata.shape)
    
    cell_types = adata.obs[label_key].unique()
    scores = []
    mmdata_per_cell_type = []

    # collect all adata subsets
    for cell_type in cell_types:
        ad_sub = adata[adata.obs[label_key] == cell_type]
        if ad_sub.obs[batch_key].nunique() <= 1:
            print(f'Skipping cell type {cell_type} because it\'s present in only one batch', flush=True)
            continue
        _mmdata = pegasusio.MultimodalData(ad_sub.copy())
        mmdata_per_cell_type.append(_mmdata)

    # compute kBET scores
    scores = list(tqdm(
        Parallel(n_jobs=n_threads)(
            delayed(compute_kbet)(
                _mmdata,
                attr=batch_key,
                rep='emb',
                K=50,
                use_cache=False,
                n_jobs=1,
            ) for _mmdata in mmdata_per_cell_type
        ),
        desc=f'Calculate kBET per cell type with {n_threads} threads',
        miniters=1,
    ))
    return np.nanmean(scores)