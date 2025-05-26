"""
Highly variable gene selection with Scanpy
"""
from pathlib import Path
import logging
logging.basicConfig(level=logging.INFO)
import warnings
warnings.filterwarnings("ignore", message="The frame.append method is deprecated and will be removed from pandas in a future version.")
import anndata as ad
from dask import array as da
from dask import config as da_config
da_config.set(num_workers=snakemake.threads)

from utils.io import read_anndata, write_zarr_linked
from utils.processing import _filter_batch, sc, USE_GPU
from utils.misc import dask_compute

input_file = snakemake.input[0]
output_file = snakemake.output[0]
args = snakemake.params.get('args', {})

if args is None:
    args = {}
elif isinstance(args, dict):
    args.pop('subset', None) # don't support subsetting
# subset_to_hvg = isinstance(args, dict) and args.get('subset', False)
logging.info(str(args))

logging.info(f'Read {input_file}...')
kwargs = dict(
    X='X',
    obs='obs',
    var='var',
    uns='uns',
    backed=True,
    dask=True,
)
# if subset_to_hvg:
#     kwargs |= dict(layers='layers')
adata = read_anndata(input_file, **kwargs)
logging.info(adata.__str__())
var = adata.var.copy()

# add metadata
if 'preprocessing' not in adata.uns:
    adata.uns['preprocessing'] = {}
adata.uns['preprocessing']['highly_variable_genes'] = args

if adata.n_obs == 0:
    logging.info('No data, write empty file...')
    adata.var['highly_variable'] = True
    adata.write_zarr(output_file)
    exit(0)

if args == False:
    logging.info('No highly variable gene parameters provided, including all genes...')
    var['highly_variable'] = True
else:
    # filter genes and cells that would break HVG function
    batch_mask = _filter_batch(adata, batch_key=args.get('batch_key'))
    adata = adata[batch_mask, adata.var['nonzero_genes']].copy()
    adata = dask_compute(adata, layers='X')
    
    # make sure data is on GPU for rapids_singlecell
    if USE_GPU:
        sc.get.anndata_to_GPU(adata)
    
    logging.info(f'Select features with arguments: {args}...')
    sc.pp.highly_variable_genes(
        adata,
        **args
    )

    # add HVG info back to adata
    hvg_column_map = {
        'highly_variable': False,
        'means': 0,
        'dispersions': 0,
        'dispersions_norm': 0,
        'highly_variable_nbatches': 0,
        'highly_variable_intersection': False,
    }
    for column, default_value in hvg_column_map.items():
        if column not in adata.var.columns:
            continue
        dtype = adata.var[column].dtype
        var[column] = default_value
        var[column] = var[column].astype(dtype)
        var.loc[adata.var_names, column] = adata.var[column]

logging.info(f'Write to {output_file}...')
files_to_keep = ['uns', 'var']
# if subset_to_hvg:
#     logging.info('Subset to highly variable genes...')
#     adata = adata[:, adata.var['highly_variable']].copy()
#     files_to_keep.extend(['X', 'layers', 'varm', 'varp'])
#     logging.info(adata.__str__())
# else:
adata = ad.AnnData(var=var, uns=adata.uns)
write_zarr_linked(
    adata,
    in_dir=input_file,
    out_dir=output_file,
    files_to_keep=files_to_keep
)