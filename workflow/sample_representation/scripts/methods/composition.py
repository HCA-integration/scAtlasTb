import logging
import warnings

import patpy as pr
import pandas as pd
import scanpy as sc

from utils.io import read_anndata, write_zarr_linked
from utils.misc import dask_compute

warnings.simplefilter("ignore", UserWarning)
logging.basicConfig(level=logging.INFO)

sc.set_figure_params(dpi=100, frameon=False)
input_file = snakemake.input.zarr
prepare_file = snakemake.input.prepare
output_file = snakemake.output.zarr
sample_key = snakemake.params.get('sample_key')
cell_type_key = snakemake.params.get('cell_type_key')
var_mask = snakemake.params.get('var_mask')

logging.info(f'Read "{input_file}"...')
n_obs = read_anndata(input_file, obs='obs').n_obs
dask = n_obs > 1e6
adata = read_anndata(
    input_file,
    obs='obs',
    backed=dask,
    dask=dask,
    stride=int(n_obs / 5),
)

# parse sample key
sample_columns = [x.strip() for x in sample_key.split(',')]
adata.obs['group'] = adata.obs[sample_columns].agg('-'.join, axis=1)
sample_key = 'group'

# subset HVGs
if var_mask is not None:
    adata.var = read_anndata(input_file, var='var').var
    adata = adata[:, adata.var[var_mask]].copy()
dask_compute(adata)

logging.info(f'Calculating composition representation for "{cell_type_key}"')
representation_method = pr.tl.CellGroupComposition(
    sample_key=sample_key,
    cell_group_key=cell_type_key,
)
representation_method.prepare_anndata(adata)
distances = representation_method.calculate_distance_matrix(force=True, dist='euclidean')

# create new AnnData object for patient representations
adata = sc.AnnData(
    obs=pd.DataFrame(index=representation_method.samples),
    obsm={
        'X_pca': sc.pp.pca(distances),
        'distances': distances,
    },
)
samples = read_anndata(prepare_file, obs='obs').obs_names
adata = adata[samples].copy()

# compute kNN graph
sc.pp.neighbors(adata, use_rep='distances', metric='precomputed', transformer='sklearn')

logging.info(f'Write "{output_file}"...')
logging.info(adata.__str__())
write_zarr_linked(
    adata,
    in_dir=prepare_file,
    out_dir=output_file,
    files_to_keep=['obsm', 'obsp', 'uns']
)
