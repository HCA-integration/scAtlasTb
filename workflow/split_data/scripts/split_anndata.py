import faulthandler
faulthandler.enable()
from pathlib import Path
import numpy as np
import warnings
warnings.filterwarnings("ignore")
import logging
logging.basicConfig(level=logging.INFO)
import anndata as ad
from pprint import pformat
from scipy.sparse import csr_matrix, coo_matrix
import sparse
from dask import array as da
from dask import config as da_config
da_config.set(num_workers=snakemake.threads)

from utils.accessors import adata_to_memory
from utils.annotate import add_wildcards
from utils.io import read_anndata, write_zarr, write_zarr_linked, ALL_SLOTS
import time

input_file = snakemake.input[0]
output_dir = snakemake.output[0]
split_key = snakemake.wildcards.key
values = snakemake.params.get('values', [])
dask = snakemake.params.get('dask', False)
write_copy = snakemake.params.get('write_copy', False) or input_file.endswith('.h5ad')
slots = snakemake.params.get('slots', {})

if not slots:
    slots = {s: s for s in ALL_SLOTS}

exclude_slots = [
    slot for slot in ALL_SLOTS
    if slot not in slots
]

assert 'obs' in slots, "obs slot must be read to split anndata"
assert 'obs' not in exclude_slots, "obs slot cannot be excluded when splitting anndata"

out_dir = Path(output_dir)
if not out_dir.exists():
    out_dir.mkdir()

# minimal slots to be read
if not write_copy and all(k == v for k, v in slots.items()):
    # only obs needed if everything else gets linked directly
    slots = dict(obs='obs')
    # TODO: also optimise case where slot key and value does not match fully

logging.info(f'Read anndata file {input_file}...')
adata = read_anndata(
    input_file,
    **slots,
    backed=True,
    dask=True,
)
logging.info(adata.__str__())

# convert split_key column to string
adata.obs[split_key] = adata.obs[split_key].astype(str)
logging.info(adata.obs[split_key].value_counts())

file_value_map = {
    s.replace(' ', '_').replace('/', '_'): s
    for s in adata.obs[split_key].astype(str).unique()
}
logging.info(f'file_value_map: {pformat(file_value_map)}')
logging.info(f'splits: {pformat(values)}')

for i, value in enumerate(values):
    split = file_value_map.get(value, value)
    out_file = out_dir / f"value~{value}.zarr"
    
    # split anndata
    logging.info(f'Split by {split_key}={split}')
    obs_mask = adata.obs[split_key] == split
    adata_sub = adata[obs_mask]
    logging.info(adata_sub.__str__())

    # add wildcards
    add_wildcards(adata_sub, {'key': split_key, 'value': split} , 'split_data')
        
    if write_copy:
        logging.info(f'Write to {out_file}...')
        start_time = time.time()
        write_zarr(adata_sub.copy(), out_file, compute=not dask)
        logging.info(f'Took {time.time() - start_time:.2f} seconds.')
    else:
        logging.info(f'Write to {out_file} with subset mask...')
        write_zarr_linked(
            adata_sub,
            in_dir=input_file,
            out_dir=out_file,
            subset_mask=(obs_mask, None),
            files_to_keep=['uns']+exclude_slots
        )
    del adata_sub
    logging.info(f'Finished {i+1} out of {len(values)}.')

logging.info(f'Finished splitting data by {split_key}.')
# touch done file
Path(snakemake.output.done).touch()