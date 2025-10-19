import numpy as np
import logging
logging.basicConfig(level=logging.INFO)

from utils.io import read_anndata, write_zarr_linked
from subset_functions import SUBSET_MAP

input_file = snakemake.input[0]
output_file = snakemake.output[0]
strategy = snakemake.params.get('strategy')
n_cell_max = snakemake.params.get('n_cells')
n_cell_max = np.iinfo(int).max if n_cell_max is None else int(n_cell_max)
sample_key = snakemake.params.get('sample_key')
per_sample = snakemake.params.get('per_sample')
random_state = snakemake.params.get('seed', 42)

files_to_keep = ['uns']

logging.info(f'Read {input_file}...')
adata = read_anndata(
    input_file,
    backed=True,
    dask=True,
)
logging.info(f'Shape before filtering: {adata.shape}')

try:
    adata.obs[sample_key]
except KeyError:
    logging.info(adata.__str__())
    raise AssertionError(f'sample key "{sample_key}" not in adata')

logging.info(adata.obs[sample_key].value_counts())

kwargs = dict(
    n_cell_max=n_cell_max,
    sample_key=sample_key,
    n_cells_per_sample=per_sample,
    random_state=random_state,
)

logging.info(f'Apply subset strategy: {strategy} with {kwargs}...')
subset_func = SUBSET_MAP[strategy]
subset_mask = subset_func(adata, **kwargs)

# subset_key = f'subset_{strategy}'
# adata.obs[subset_key] = False
# adata.obs.loc[subset_mask, subset_key] = True

n_cells = subset_mask.sum()
n_samples = len(adata.obs.loc[subset_mask, sample_key].unique())

logging.info(f'Subset to {n_cells} cells, {n_samples} samples...')
adata = adata[subset_mask].copy()
adata.uns['subset'] = {
    'strategy': strategy,
    'n_cells': n_cell_max,
    'sample_key': sample_key,
    'per_sample': per_sample,
    'seed': random_state,
}
logging.info(adata.__str__())

# save
logging.info('Write...')
write_zarr_linked(
    adata,
    in_dir=input_file,
    out_dir=output_file,
    files_to_keep=files_to_keep,
    compute=True,
    subset_mask=(subset_mask, adata.var_names.isin(adata.var_names)),
)
