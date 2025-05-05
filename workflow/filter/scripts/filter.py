"""
Data filtering
"""
from pathlib import Path
from pprint import pformat
import numpy as np
import pandas as pd
import scanpy as sc
from dask import array as da
import logging
logging.basicConfig(level=logging.INFO)

from utils.io import read_anndata, write_zarr_linked, get_file_reader
from utils.misc import dask_compute


input_file = snakemake.input[0]
output_file = snakemake.output[0]
params = dict(snakemake.params)

backed = params.get('backed', True)
dask = params.get('dask', True)
subset = params.get('subset', False)

kwargs = {'obs': 'obs', 'var': 'var'}
adata = read_anndata(input_file, **kwargs)
logging.info(adata.__str__())

mask = pd.Series(np.full(adata.n_obs, True, dtype=bool), index=adata.obs_names)
subset_mask = None

ex_filters = params.get('remove_by_column', {})
logging.info(pformat(ex_filters))
for column, values in ex_filters.items():
    logging.info(f'remove cells matching {len(values)} value(s) from column="{column}"...')
    values = [str(v) for v in values]
    mask &= ~adata.obs[column].astype(str).isin(values)

for query in params.get('remove_by_query', []):
    logging.info(f'remove by query="{query}"...')
    mask &= adata.obs.eval(query)

logging.info('Add filtered column...')
adata.obs['filtered'] = mask
value_counts = adata.obs['filtered'].value_counts()
logging.info(value_counts)

if subset and False in value_counts.index:
    kwargs |= {
        'X': 'X',
        'layers': 'layers',
        'raw': 'raw',
        'obsm': 'obsm',
        'obsp': 'obsp',
    }
    # filter out slots that aren't present in the input
    func, _ = get_file_reader(input_file)
    store = func(input_file, 'r')
    kwargs = {k: v for k, v in kwargs.items() if k in store.keys()}
    
    logging.info('Read all slots for subsetting...')
    obs = adata.obs # save updated obs
    adata = read_anndata(
        input_file,
        backed=backed,
        dask=dask,
        **{k: v for k, v in kwargs.items() if k != 'obs'},
    )
    adata.obs = obs # updated obs
    
    logging.info('Subset data by filters...')
    adata = dask_compute(adata[adata.obs['filtered']].copy())
    logging.info(adata.__str__())
    files_to_keep = kwargs.keys()
else:
    var_mask = np.full(adata.n_vars, True, dtype=bool)
    subset_mask = (mask.values, var_mask)
    files_to_keep = []

logging.info(f'Write to {output_file}...')
write_zarr_linked(
    adata,
    input_file,
    output_file,
    files_to_keep=files_to_keep,
    subset_mask=subset_mask,
)