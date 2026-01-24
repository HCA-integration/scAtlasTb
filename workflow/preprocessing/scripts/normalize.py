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

from utils.processing import sc, USE_GPU
from utils.io import read_anndata, write_zarr_linked
from utils.misc import ensure_sparse


input_file = snakemake.input[0]
output_file = snakemake.output[0]
layer = snakemake.params.get('raw_counts', 'X')
if not layer:
    layer = 'X'
gene_id_column = snakemake.params.get('gene_id_column')
args = snakemake.params.get('args', {})
dask = snakemake.params.get('dask', True) # get global dask flag
dask = args.pop('dask', dask) # overwrite with normalize-specific dask flag

logging.info(f'Read {input_file} with X={layer}...')
adata = read_anndata(
    input_file,
    X=layer,
    obs='obs',
    var='var',
    uns='uns',
    backed=dask,
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

# make sure data is sparse
logging.info('ensure sparse...')
ensure_sparse(adata)

if input_file.endswith('.h5ad'):
    logging.info('Copy counts to layers...')
    adata.layers['counts'] = adata.X
    adata.raw = adata

# make sure data is on GPU for rapids_singlecell
if USE_GPU:
    logging.info('Transfer to GPU...')
    # adata.X = adata.X.astype('float32')
    sc.get.anndata_to_GPU(adata)

logging.info(f'normalize_total with args={args}...')
sc.pp.normalize_total(adata, **args)
logging.info('log-transform...')
sc.pp.log1p(adata)

if USE_GPU:
    logging.info('Transfer to CPU...')
    sc.get.anndata_to_CPU(adata)

adata.layers['normcounts'] = adata.X

# add preprocessing metadata
if 'preprocessing' not in adata.uns:
    adata.uns['preprocessing'] = {}

adata.uns['preprocessing']['normalization'] = args
adata.uns['preprocessing']['log-transformed'] = True
# scanpy.pp.log1p was supposed to add it but it's not saved
adata.uns["log1p"] = {"base": None}

logging.info(f'Write to {output_file}...')
logging.info(adata.__str__())

if not input_file.endswith('.h5ad'):
    del adata.X

write_zarr_linked(
    adata,
    input_file,
    output_file,
    files_to_keep=['uns', 'var', 'layers/normcounts'],
    slot_map={
        'layers/counts': layer,
        'X': 'layers/normcounts',
        'raw/X': layer,
        'raw/var': 'var',
    },
    in_dir_map={
        'layers/normcounts': output_file,
    },
)