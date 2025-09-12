from pathlib import Path
import scanpy as sc
import warnings
warnings.simplefilter("ignore", UserWarning)
import logging
logging.basicConfig(level=logging.INFO)
# import patpy as pr

from utils.io import read_anndata, write_zarr, write_zarr_linked, ALL_SLOTS
from utils.misc import dask_compute
from utils.processing import get_pseudobulks


sc.set_figure_params(dpi=100, frameon=False)
input_zarr = snakemake.input.zarr
output_zarr = snakemake.output.zarr
output_bulk_zarr = snakemake.output.bulks

sample_key = snakemake.params.get('sample_key')
cell_type_key = snakemake.params.get('cell_type_key')
aggregate = snakemake.params.get('aggregate')
layer = snakemake.params.get('layer')
min_cells_per_sample = snakemake.params.get('min_cells_per_sample')
min_cells_per_cell_type = snakemake.params.get('min_cells_per_cell_type')

logging.info(f'Read "{input_zarr}"...')
n_obs = read_anndata(input_zarr, X=layer, dask=True, backed=True, verbose=False).n_obs
dask = n_obs > 5e5

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
adata_bulk = get_pseudobulks(
    adata,
    group_key='group',
    agg=aggregate,
)

# preprocess pseudobulks
sc.pp.normalize_total(adata_bulk)
sc.pp.log1p(adata_bulk)
logging.info(adata_bulk.__str__())

logging.info(f'Write "{output_bulk_zarr}"...')
write_zarr(adata_bulk, output_bulk_zarr)

logging.info(f'Write "{output_zarr}"...')

if Path(input_zarr).suffix == '.h5ad':
    obs = adata.obs
    del adata
    logging.info('Read full h5ad file...')
    adata = read_anndata(input_zarr, dask=True, backed=True, verbose=False)
    adata.obs = obs

write_zarr_linked(
    adata,
    in_dir=input_zarr,
    out_dir=output_zarr,
    files_to_keep=['obs'],
)