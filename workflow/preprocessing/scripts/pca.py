"""
PCA on highly variable genes
"""
from pathlib import Path
import numpy as np
import logging
logging.basicConfig(level=logging.INFO)
from dask import array as da
from dask import config as da_config
da_config.set(num_workers=snakemake.threads)
import anndata as ad

from utils.accessors import subset_hvg
from utils.io import read_anndata, write_zarr_linked
from utils.misc import dask_compute
from utils.processing import sc, USE_GPU


input_file = snakemake.input[0]
output_file = snakemake.output[0]
scale = snakemake.params.get('scale', False)
args = snakemake.params.get('args', {})
dask = snakemake.params.get('dask', True) # get global dask flag
dask = args.pop('dask', dask) # overwrite with pca-specific dask flag

logging.info(f'Read "{input_file}"...')
adata = read_anndata(
    input_file,
    X='X',
    var='var',
    obs='obs',
    uns='uns',
    dask=True,
    backed=True,
)
logging.info(adata.__str__())

# parse arguments
if 'zero_center' in args:
    if args['zero_center'] == 'None':
        args['zero_center'] = None
args['svd_solver'] = args.get('svd_solver', 'covariance_eigh')

# add preprocessing arguments
if 'preprocessing' not in adata.uns:
    adata.uns['preprocessing'] = {}

adata.uns['preprocessing']['scaled'] = scale
adata.uns['preprocessing']['pca'] = args

if adata.n_obs == 0:
    logging.info('No data, write empty file...')
    adata.obsm['X_pca'] = np.zeros((0, 30))
    adata.write_zarr(output_file)
    exit(0)

logging.info('Subset to highly variable genes...')
hvg_key = args.pop('mask_var', 'highly_variable')
var = adata.var.copy()
adata_pca, _ = subset_hvg(adata, var_column=hvg_key, compute_dask=not dask)
adata = ad.AnnData(obs=adata.obs, var=var)

if USE_GPU:
    sc.get.anndata_to_GPU(adata_pca)

# scaling TODO: move to separate rule
if scale:
    logging.info('Scale data...')
    sc.pp.scale(adata_pca)

logging.info(f'PCA with parameters: {args}')
sc.pp.pca(adata_pca, **args)
adata.obsm['X_pca'] = adata_pca.obsm['X_pca']
adata.uns['pca'] = adata_pca.uns['pca']

# Map all entries in varm back to all genes
for key, varm in adata_pca.varm.items():
    full_arr = np.full(
        (adata.n_vars, varm.shape[1]),
        np.nan,
        dtype=varm.dtype
    )
    full_arr[adata.var_names.isin(adata_pca.var_names)] = varm
    adata.varm[key] = full_arr

del adata_pca

logging.info(f'Write to {output_file}...')
logging.info(adata.__str__())
adata = dask_compute(adata, layers='X_pca')
write_zarr_linked(
    adata,
    input_file,
    output_file,
    files_to_keep=['obsm', 'uns', 'varm']
)
