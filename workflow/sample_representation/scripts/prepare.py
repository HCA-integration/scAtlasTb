import numpy as np
from matplotlib import pyplot as plt
import anndata
import pandas as pd
import scanpy as sc
import warnings
warnings.simplefilter("ignore", UserWarning)
import logging
logging.basicConfig(level=logging.INFO)

import patpy as pr

from utils.io import read_anndata
from utils.processing import get_pseudobulks


sc.set_figure_params(dpi=100, frameon=False)
input_zarr = snakemake.input.zarr
output_zarr = snakemake.output.zarr
sample_key = snakemake.params.get('sample_key')
cell_type_key = snakemake.params.get('cell_type_key')
aggregate = snakemake.params.get('aggregate')
layer = snakemake.params.get('norm_counts')
min_cells_per_sample = snakemake.params.get('min_cells_per_sample')
min_cells_per_cell_type = snakemake.params.get('min_cells_per_cell_type')

logging.info(f'Read "{input_zarr}"...')
n_obs = read_anndata(input_zarr, X=layer, dask=True, backed=True, verbose=False).n_obs
dask = n_obs > 1e6
adata = read_anndata(
    input_zarr,
    X=layer,
    obs='obs',
    var='var',
    backed=dask,
    dask=dask,
    stride=int(n_obs / 5),
)

# # filter small samples and cell types
# adata = pr.pp.filter_small_samples(
#     adata,
#     sample_key,
#     sample_size_threshold=min_cells_per_sample,
# )
# if cell_type_key is not None:
#     adata = pr.pp.filter_small_cell_types(
#         adata,
#         sample_key,
#         cell_type_key,
#         cluster_size_threshold=min_cells_per_cell_type,
#     )

logging.info(f'Pseudobulk by "{sample_key}"...')
sample_key = [x.strip() for x in sample_key.split(',')]
adata.obs['group'] = adata.obs[sample_key].astype(str).agg('-'.join, axis=1)
adata_bulk = get_pseudobulks(adata, group_key='group', agg=aggregate)
sc.pp.normalize_total(adata_bulk)
sc.pp.log1p(adata_bulk)

logging.info(f'Write "{output_zarr}"...')
logging.info(adata_bulk.__str__())
adata_bulk.write_zarr(output_zarr)
