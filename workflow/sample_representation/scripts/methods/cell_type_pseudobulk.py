import logging
import warnings

import patpy as pr
import numpy as np
import pandas as pd
import scanpy as sc

from utils.io import read_anndata, write_zarr_linked
from utils.misc import dask_compute

warnings.simplefilter("ignore", UserWarning)
logging.basicConfig(level=logging.INFO)

sc.set_figure_params(dpi=100, frameon=False)
input_file = snakemake.input.zarr
bulk_file = snakemake.input.bulks
output_file = snakemake.output.zarr

cell_type_key = snakemake.params.get('cell_type_key')
use_rep = snakemake.params.get('use_rep')
var_mask = snakemake.params.get('var_mask')

logging.info(f'Read "{input_file}"...')
n_obs = read_anndata(input_file, obs='obs').n_obs
dask = n_obs > 2e6
adata = read_anndata(
    input_file,
    X=use_rep,
    obs='obs',
    backed=dask,
    dask=dask,
    stride=int(n_obs / 5),
)

# subset HVGs
if var_mask is not None:
    adata.var = read_anndata(input_file, var='var').var
    adata = adata[:, adata.var[var_mask]].copy()
dask_compute(adata)

logging.info(f'Calculating Cell Type Pseudobulk representation for "{cell_type_key}", using cell features from "{use_rep}"')
representation_method = pr.tl.GroupedPseudobulk(
    sample_key='group',
    cell_group_key=cell_type_key,
    layer='X',
)
representation_method.prepare_anndata(adata)  # Assuming that small cell types and samples are already filtered out

# compute distances
distances = representation_method.calculate_distance_matrix(
    force=True,
    aggregate="mean",
    dist="euclidean"
)

# create new AnnData object for patient representations
adata = sc.AnnData(
    obs=pd.DataFrame(index=representation_method.samples),
    obsm={'distances': distances} | {
        cell_type: representation_method.sample_representation[i]
        for i, cell_type in enumerate(representation_method.cell_groups)
    }
)
adata.obsm['X_emb'] = np.hstack(representation_method.sample_representation)
adata.obsm['X_pca'] = sc.pp.pca(adata.obsm['X_emb'])

samples = read_anndata(bulk_file, obs='obs').obs_names
adata = adata[samples].copy()

# compute kNN graph
sc.pp.neighbors(adata, use_rep='distances', metric='precomputed', transformer='sklearn')
# sc.pp.neighbors(adata, use_rep='X_emb', key_added='X_emb')

logging.info(f'Write "{output_file}"...')
logging.info(adata.__str__())
write_zarr_linked(
    adata,
    in_dir=bulk_file,
    out_dir=output_file,
    files_to_keep=['obsm', 'obsp', 'uns'],
)
