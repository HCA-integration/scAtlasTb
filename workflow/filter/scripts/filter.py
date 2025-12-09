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

subset = params.get('subset', True)
write_copy = params.get('write_copy', False) or input_file.endswith('.h5ad')

adata = read_anndata(input_file, dask=True, backed=True)
logging.info(adata.__str__())

mask = pd.Series(np.full(adata.n_obs, True, dtype=bool), index=adata.obs_names)

ex_filters = params.get('remove_by_column', {})
logging.info(pformat(ex_filters))
for column, values in ex_filters.items():
    logging.info(f'remove cells matching {len(values)} value(s) from column="{column}"...')
    values = [str(v) for v in values]
    mask &= ~adata.obs[column].astype(str).isin(values)
    logging.info(f'{mask.sum()} cells remaining')

for query in params.get('keep_by_query', []):
    logging.info(f'keep by query="{query}"...')
    mask &= adata.obs.eval(query)
    logging.info(f'{mask.sum()} cells remaining')

for query in params.get('remove_by_query', []):
    logging.info(f'remove by query="{query}"...')
    mask &= ~adata.obs.eval(query)
    logging.info(f'{mask.sum()} cells remaining')

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

    logging.info(f'Write to {output_file}...')
    if write_copy:
        write_zarr(
            adata,
            output_file,
            compute=True,
        )
    else:
        write_zarr_linked(
            adata,
            input_file,
            output_file,
            files_to_keep=['obs'],
            subset_mask=(mask.values, None),
        )