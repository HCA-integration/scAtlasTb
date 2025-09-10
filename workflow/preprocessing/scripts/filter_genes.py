"""
Highly variable gene selection with Scanpy
"""
from pathlib import Path
import logging
logging.basicConfig(level=logging.INFO)
import anndata as ad
from dask import array as da
from dask import config as da_config
da_config.set(num_workers=snakemake.threads)

from utils.io import read_anndata, write_zarr_linked
from utils.processing import _filter_genes

input_file = snakemake.input[0]
output_file = snakemake.output[0]
dask = snakemake.params.get('dask', True)

logging.info(f'Read {input_file}...')
kwargs = dict(
    X='X',
    obs='obs',
    var='var',
    uns='uns',
    backed=dask,
    dask=dask,
    stride=10_000,
)
# TODO: read with chunks (None, 'auto') for optimal gene-wise operations

adata = read_anndata(input_file, **kwargs)
logging.info(adata.__str__())
var = adata.var.copy()

logging.info('Determine nonzero genes...')
adata.var['nonzero_genes'] = _filter_genes(adata, min_cells=1)

write_zarr_linked(
    adata,
    in_dir=input_file,
    out_dir=output_file,
    files_to_keep=['var/nonzero_genes'],
)