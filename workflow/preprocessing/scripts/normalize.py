"""
Normalisation
"""
from pathlib import Path
import numpy as np
import logging
logging.basicConfig(level=logging.INFO)
from dask import array as da
from dask import config as da_config
da_config.set(num_workers=snakemake.threads)
import sparse
try:
    import rapids_singlecell as sc
    import cupy as cp
    logging.info('Using rapids_singlecell...')
    rapids = True
except ImportError as e:
    import scanpy as sc
    logging.info('Importing rapids failed, using scanpy...')
    rapids = False

from utils.io import read_anndata, write_zarr_linked, csr_matrix_int64_indptr
from utils.misc import apply_layers, ensure_sparse


input_file = snakemake.input[0]
output_file = snakemake.output[0]
layer = snakemake.params.get('raw_counts', 'X')
gene_id_column = snakemake.params.get('gene_id_column')
args = snakemake.params.get('args', {})
dask = snakemake.params.get('dask', False) and not rapids
backed = snakemake.params.get('backed', False) and dask and not rapids

logging.info(f'Read {input_file}...')
adata = read_anndata(
    input_file,
    X=layer,
    obs='obs',
    var='var',
    uns='uns',
    backed=backed,
    dask=dask,
)
logging.info(adata.__str__())

if gene_id_column is not None:
    adata.var_names = adata.var[gene_id_column]
adata.var_names = adata.var_names.astype(str).values
adata.var_names_make_unique()

# deal with empty files
if adata.n_obs == 0:
    logging.info('No data, write empty file...')
    adata.X = np.zeros(adata.shape)
    adata.write_zarr(output_file)
    exit(0)

if input_file.endswith('.h5ad'):
    logging.info('Copy counts to layers...')
    adata.layers['counts'] = adata.X
    adata.raw = adata

if isinstance(adata.X, da.Array):
    logging.info('Convert dask chunks to sparse chunks...')
    adata.X = adata.X.map_blocks(lambda x: x.toarray(), dtype=adata.X.dtype)
    logging.info(adata.X)

# make sure data is on GPU for rapids_singlecell
if rapids:
    logging.info('Transfer to GPU...')
    # adata.X = adata.X.astype('float32')
    sc.get.anndata_to_GPU(adata)

logging.info('normalize_total with args={args}...')
sc.pp.normalize_total(adata, **args)
logging.info('log-transform...')
sc.pp.log1p(adata)

# make sure data is sparse
logging.info('ensure sparse...')
ensure_sparse(adata, sparse_type=csr_matrix_int64_indptr)

if rapids:
    logging.info('Transfer to CPU...')
    sc.get.anndata_to_CPU(adata)

adata.layers['normcounts'] = adata.X

# add preprocessing metadata
if 'preprocessing' not in adata.uns:
    adata.uns['preprocessing'] = {}

adata.uns['preprocessing']['normalization'] = 'default'
adata.uns['preprocessing']['log-transformed'] = True
# scanpy.pp.log1p was supposed to add it but it's not saved
adata.uns["log1p"] = {"base": None}

logging.info(f'Write to {output_file}...')
logging.info(adata.__str__())
write_zarr_linked(
    adata,
    input_file,
    output_file,
    files_to_keep=['uns', 'var', 'layers/normcounts'],
    slot_map={
        'X': 'layers/normcounts',
        'layers/counts': layer,
        'raw/X': layer,
        'raw/var': 'var',
    },
    in_dir_map={
        'layers/normcounts': output_file,
    },
)