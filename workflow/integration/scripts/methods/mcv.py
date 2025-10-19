from tqdm import tqdm
from pathlib import Path
from pprint import pformat
import numpy as np
from scipy.sparse import csr_matrix, issparse
from sklearn.metrics import mean_squared_error
import anndata as ad
from dask import array as da
# from utils.misc import ensure_dense, dask_compute
from utils.processing import sc, USE_GPU
rsc = sc
import scanpy as sc


def plot_mcv_pca(mcv_summary, figdir, prefix='', figsize=(10, 3)):
    """
    plot_mcv_pca(mcv_summary, figdir=figdir)
    """
    from matplotlib import pyplot as plt
    
    plt.figure(figsize=figsize)
    plt.plot(mcv_summary['k_range'], mcv_summary['mcv_loss'])
    optimal_k = mcv_summary['k_range'][np.argmin(mcv_summary['mcv_loss'])]
    plt.scatter([optimal_k], [mcv_summary['mcv_loss'][np.argmin(mcv_summary['mcv_loss'])]], c='k')
    plt.xlabel('Number of PCs')
    plt.ylabel('MCV Loss')
    plt.title(f'Optimal PCs = {optimal_k}')
    plt.tight_layout()
    plt.savefig(Path(figdir) / f"{prefix}pca_mcv.png")
    plt.close()


def mse_sparse(a, b):
    if a.shape != b.shape:
        raise ValueError("Shapes of input matrices must match.")

    if issparse(a) and issparse(b):
        # both sparse
        diff = a - b
        # Compute squared error for all entries
        # For sparse matrices, .power(2) squares all entries
        return diff.power(2).sum() / (a.shape[0] * a.shape[1])
    elif issparse(a):
        # a is sparse, b is dense
        diff = a.toarray() - b
        return np.mean(np.square(diff))
    elif issparse(b):
        # b is sparse, a is dense
        diff = a - b.toarray()
        return np.mean(np.square(diff))
    return mean_squared_error(a, b)


def mcv_optimal_pcs_scanpy(adata, raw_name, max_pcs=100, scale_for_pca=True):
    """
    Molecular cross-validation to select optimal number of PCs with Scanpy.
    Optimized for speed and memory by avoiding densification where possible.
    """
    X = adata.layers[raw_name]

    print("MCV: Splitting train/test data", flush=True)
    np.random.seed(42)
    if hasattr(X, "compute"):  # X is a dask array
        # Avoid full densification; work with chunks
        data = X.map_blocks(lambda arr: arr.data.astype(np.int32), dtype=np.int32)
        data1 = data.map_blocks(lambda arr: np.random.binomial(arr, 0.5), dtype=np.int32)
        data2 = data - data1
        del data
        indices = X.map_blocks(lambda arr: arr.indices, dtype=np.int32)
        indptr = X.map_blocks(lambda arr: arr.indptr, dtype=np.int64)

        # Compute only the final arrays
        data1 = data1.compute()
        data2 = data2.compute()
        indices = indices.compute()
        indptr = indptr.compute()
    else:
        data = X.data.astype(np.int32)
        data1 = np.random.binomial(data, 0.5)
        data2 = data - data1
        del data
        indices, indptr = X.indices, X.indptr

    adata1 = ad.AnnData(
        csr_matrix((data1, indices, indptr), shape=X.shape),
        var=adata.var
    )
    adata2 = ad.AnnData(
        csr_matrix((data2, indices, indptr), shape=X.shape),
        var=adata.var
    )

    for ad_sub in tqdm(
        [adata1, adata2],
        miniters=1,
        desc=f"MCV:Preprocessing, scale={scale_for_pca}"
    ):
        sc.pp.normalize_total(ad_sub, target_sum=1e4)
        sc.pp.log1p(ad_sub)
        if scale_for_pca:
            # ad_sub.X = da.from_array(ad_sub.X)
            sc.pp.scale(ad_sub, zero_center=False)

    print(f"MCV: Calculating {max_pcs} PCs", flush=True)
    sc.pp.pca(adata1, n_comps=max_pcs, svd_solver='covariance_eigh')

    # define search space
    k_range = np.concatenate([
        np.arange(2, 10, 1),
        np.arange(10, 30, 2),
        np.arange(30, max_pcs + 1, 5)
    ])
    k_range.sort()
    print(f'MCV: k_range = {pformat(k_range)}, n={len(k_range)}', flush=True)

    prev_k = 0
    recon = csr_matrix(adata1.shape, dtype=np.float32)
    # recon = np.zeros(adata1.shape, dtype=np.float32)
    mcv_loss = np.zeros(len(k_range), dtype=np.float32)
    rec_loss = np.zeros(len(k_range), dtype=np.float32)
    
    # X1_dense = adata1.X.toarray() if issparse(adata1.X) else adata1.X
    # X2_dense = adata2.X.toarray() if issparse(adata2.X) else adata2.X
    
    X_pca = csr_matrix(adata1.obsm['X_pca'])  # shape (n_obs, max_pcs)
    PCs = csr_matrix(adata1.varm['PCs'])     # shape (n_vars, max_pcs)

    message = 'MCV: Computing reconstruction losses for k_range'
    for i, k in enumerate(tqdm(k_range, desc=message, miniters=1)):
    
        # Incrementally build from previous k to current k
        recon += X_pca[:, prev_k:k] @ PCs[:, prev_k:k].T
            
        # compute losses
        mcv_loss[i] = mse_sparse(recon, adata2.X)
        rec_loss[i] = mse_sparse(recon, adata1.X)
        # mcv_loss[i] = mean_squared_error(recon, X2_dense)
        # rec_loss[i] = mean_squared_error(recon, X1_dense)
        
        # update previous k
        prev_k = k

    optimal_k = k_range[np.argmin(mcv_loss)]
    mcv_summary = {
        'k_range': k_range,
        'mcv_loss': mcv_loss,
        'rec_loss': rec_loss
    }

    return optimal_k, mcv_summary
