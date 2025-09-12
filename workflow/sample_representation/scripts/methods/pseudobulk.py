import logging
import warnings

from dask import array as da
import patpy as pr
import pandas as pd
import scanpy as sc

from utils.io import read_anndata, write_zarr_linked
from utils.misc import dask_compute
from utils.processing import get_pseudobulks
from utils.annotate import add_wildcards

warnings.simplefilter("ignore", UserWarning)
logging.basicConfig(level=logging.INFO)

input_file = snakemake.input.zarr
bulk_file = snakemake.input.bulks
output_file = snakemake.output.zarr

use_rep = snakemake.params.get('use_rep')
var_mask = snakemake.params.get('var_mask')

logging.info(f'Read "{input_file}"...')
if use_rep in [None, 'X'] or use_rep.startswith('layers/'):
    adata = read_anndata(
        bulk_file,
        obs='obs',
        X='X',
        var='var',
        uns='uns'
    )
    # subset HVGs
    if var_mask is not None:
        adata.var = read_anndata(input_file, var='var').var
        adata = adata[:, adata.var[var_mask]].copy()
    sc.pp.pca(adata)
    adata.obsm['X_emb'] = adata.obsm['X_pca']
else:
    n_obs, n_vars = read_anndata(input_file, X=use_rep, dask=True, backed=True, verbose=False).shape
    dask = n_obs > 1e6 and n_vars > 500
    adata = read_anndata(
        input_file,
        X=use_rep,
        obs='obs',
        backed=dask,
        dask=dask,
        stride=int(n_obs / 5),
    )
    # # make sure all values are positive integers
    # _min = adata.X.min()
    # if isinstance(_min, da.Array):
    #     _min = _min.compute()
    # if _min < 0:
    #     adata.X -= _min

    logging.info('Aggregating"...')
    adata = get_pseudobulks(adata, group_key='group', agg='mean', group_cols=[])
    adata.obsm['X_emb'] = adata.X
    sc.pp.pca(adata)

    samples = read_anndata(bulk_file, obs='obs').obs_names
    adata = adata[samples].copy()

# compute kNN graph
# sc.pp.neighbors(adata, use_rep='distances', metric='precomputed', transformer='sklearn')
# sc.pp.neighbors(adata, use_rep='X_emb', key_added='X_emb')
sc.pp.neighbors(adata, use_rep='X_emb')

adata.uns['output_type'] = 'embed'
add_wildcards(adata, snakemake.wildcards, 'sample_representation')
logging.info(f'Write "{output_file}"...')
logging.info(adata.__str__())
write_zarr_linked(
    adata,
    in_dir=bulk_file,
    out_dir=output_file,
    files_to_keep=['obsm', 'obsp', 'uns']
)
