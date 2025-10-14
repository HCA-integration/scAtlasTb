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

from utils.io import read_anndata, write_zarr, write_zarr_linked, get_file_reader
from utils.misc import dask_compute


input_file = snakemake.input[0]
output_file = snakemake.output[0]
params = dict(snakemake.params)

backed = params.get('backed', True)
dask = params.get('dask', True)
subset = params.get('subset', True)

adata = read_anndata(input_file, dask=dask, backed=backed)
logging.info(adata.__str__())

mask = pd.Series(np.full(adata.n_obs, True, dtype=bool), index=adata.obs_names)

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

# update subset flag
subset = subset and False in value_counts.index

if not subset:
    # don't subset, just keep filtering annotation
    logging.info(f'Write to {output_file}...')
    write_zarr_linked(
        adata,
        input_file,
        output_file,
        files_to_keep=['obs'],
    )
else:
    logging.info('Subset data by filters...')
    adata = adata[adata.obs['filtered']].copy()
    logging.info(adata.__str__())
    
    var_mask = np.full(adata.n_vars, True, dtype=bool)

    logging.info(f'Write to {output_file}...')
    write_zarr_linked(
        adata,
        input_file,
        output_file,
        files_to_keep=['obs'],
        compute=True, # for h5ad files, for zarr all slots other than obs will be dropped before writing
        subset_mask=(mask.values, var_mask),
    )