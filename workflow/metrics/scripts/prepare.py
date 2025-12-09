from pathlib import Path
import logging
from pprint import pformat
import numpy as np
import anndata as ad
import scanpy as sc

from utils.assertions import assert_pca
from utils.accessors import subset_hvg
from utils.io import read_anndata, write_zarr_linked, write_zarr
from utils.processing import compute_neighbors, _filter_genes
from utils.misc import dask_compute, apply_layers


def compute_pca(adata, matrix):
    X_pca, _, variance_ratio, variance = sc.tl.pca(matrix, return_info=True)
    adata.obsm['X_pca'] = X_pca
    adata.uns['pca'] = {
        'variance_ratio': variance_ratio,
        'variance': variance,
    }


input_file = snakemake.input[0]
output_file = snakemake.output[0]
params = snakemake.params
# label_key = params.get('label_key')
neighbor_args = params.get('neighbor_args', {})
unintegrated_layer = params.get('unintegrated_layer', 'X')
corrected_layer = params.get('corrected_layer', 'X')
var_key = params.get('var_mask', 'highly_variable')
output_type = params.get('output_type', 'embed')

PERSIST_MATRIX_THRESHOLD = params.get('PERSIST_MATRIX_THRESHOLD', 5e5)

files_to_keep = ['raw', 'uns', 'var']
slot_map = {'raw/X': unintegrated_layer}

# determine output types
# default output type is 'full'
output_type = read_anndata(input_file, uns='uns').uns.get('output_type', output_type)

logging.info(f'Output type: {output_type}')
kwargs = dict(
    obs='obs',
    obsm='obsm',
    obsp='obsp',
    var='var',
    uns='uns',
)

is_h5ad = input_file.endswith('.h5ad')
if is_h5ad:
    kwargs |= dict(
        X=corrected_layer,
        layers='layers',
        dask=True,
        backed=True
    )
    files_to_keep.extend(['obs', 'X', 'layers'])

if output_type == 'full':
    kwargs |= dict(X=corrected_layer, dask=True, backed=True)
    slot_map |= dict(X=corrected_layer)

logging.info(f'Read {input_file}...')
adata = read_anndata(input_file, **kwargs)

if is_h5ad:
    adata.layers['unintegrated'] = read_anndata(
        input_file,
        X=unintegrated_layer,
        dask=True,
        backed=True,
    ).X

all_obs_names = adata.obs_names.copy()
all_var_names = adata.var_names.copy()

# remove cells without labels
n_obs = adata.n_obs
# logging.info('Filtering out cells without labels')
# TODO: only for metrics that require labels?
# logging.info(f'Before: {adata.n_obs} cells')
# adata = adata[(adata.obs[label_key].notna() | adata.obs[label_key] != 'NA') ]
# logging.info(f'After: {adata.n_obs} cells')
# if adata.is_view:
#     adata = adata.copy()
force_neighbors = n_obs > adata.n_obs

# set HVGs
new_var_column = 'metrics_features'
if var_key in adata.var.columns.values:
    adata.var[new_var_column] = adata.var[var_key]
else:
    logging.info(f'"{var_key}" not in adata var, setting all to True\n{pformat(adata.var.columns)}')
    adata.var[new_var_column] = True
logging.info(f'Set "{new_var_column}" from "{var_key}" in adata.var: {adata.var[new_var_column].sum()} HVGs')

# logging.info('Filter all zero genes...')
# all_zero_genes = _filter_genes(adata, min_cells=1)
# adata.var[new_var_column] = adata.var[new_var_column] & ~adata.var_names.isin(all_zero_genes)

logging.info('Compute PCA...')
if output_type == 'full':
    sc.pp.pca(
        adata,
        mask_var=new_var_column,
        svd_solver='covariance_eigh',
    )
    adata = dask_compute(adata, layers='X_pca')
    adata.obsm['X_emb'] = adata.obsm['X_pca']
    force_neighbors = True
    files_to_keep.append('obsm/X_pca')
elif output_type == 'embed':
    logging.info('Run PCA on embedding...')
    compute_pca(adata, matrix=adata.obsm['X_emb'])
    files_to_keep.append('obsm/X_pca')

if force_neighbors:
    files_to_keep.append('obsp')

logging.info(f'Computing neighbors for output type {output_type} force_neighbors={force_neighbors}...')
compute_neighbors(
    adata,
    output_type,
    force=force_neighbors,
    check_n_neighbors=False,
    **neighbor_args
)

# unintegrated for comparison metrics
# if output_type != 'knn':  # assuming that there aren't any knn-based metrics that require PCA
logging.info(f'Prepare unintegrated data from layer={unintegrated_layer}...')
adata_raw = read_anndata(
    input_file,
    X=unintegrated_layer,
    dask=True,
    backed=True,
)
adata_raw.var = adata.var
logging.info(f'Unintegrated data shape: {adata_raw.shape}') 

if adata_raw.n_obs > PERSIST_MATRIX_THRESHOLD:
    logging.info('Persist matrix...')
    adata_raw = apply_layers(adata_raw, lambda x: x.persist(), layers='X')
else:
    dask_compute(adata_raw, layers='X')

logging.info('Run PCA on unintegrated data...')
sc.pp.pca(
    adata_raw,
    mask_var=new_var_column,
    svd_solver='covariance_eigh',
)

# add raw PCA adata
adata.obsm['X_pca_unintegrated'] = adata_raw.obsm['X_pca']
del adata_raw.obsm['X_pca']
adata = dask_compute(adata, layers='X_pca_unintegrated')
adata.uns['pca_unintegrated'] = adata_raw.uns['pca']
files_to_keep.append('obsm/X_pca_unintegrated')

# add raw data for comparison
adata.raw = adata_raw

logging.info(f'Write to {output_file}...')
logging.info(adata.__str__())
if is_h5ad:
    write_zarr(adata, output_file)
else:
    write_zarr_linked(
        adata,
        in_dir=input_file,
        out_dir=output_file,
        files_to_keep=files_to_keep,
        slot_map=slot_map,
    )
