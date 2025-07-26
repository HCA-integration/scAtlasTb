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
from utils.io import read_anndata, write_zarr, write_zarr_linked
from utils.misc import dask_compute

input_file = snakemake.input[0]
output_dir = snakemake.output[0]
split_key = snakemake.wildcards.key
values = snakemake.params.get('values', [])
backed = snakemake.params.get('backed', True)
dask = snakemake.params.get('dask', True)
write_copy = snakemake.params.get('write_copy', False)
exclude_slots = snakemake.params.get('exclude_slots', [])

out_dir = Path(output_dir)
if not out_dir.exists():
    out_dir.mkdir()

write_copy = write_copy or input_file.endswith('.h5ad')
# kwargs = dict(obs='obs', var='var', uns='uns')
kwargs = dict(
    backed=backed,
    dask=dask,
    exclude_slots=exclude_slots,
)

logging.info(f'Read anndata file {input_file}...')
adata = read_anndata(input_file, **kwargs)
logging.info(adata.__str__())

# convert split_key column to string
adata.obs[split_key] = adata.obs[split_key].astype(str)
logging.info(adata.obs[split_key].value_counts())

file_value_map = {
    s.replace(' ', '_').replace('/', '_'): s
    for s in adata.obs[split_key].astype(str).unique()
}
logging.info(f'file_value_map: {pformat(file_value_map)}')
# split_files = list(file_value_map.keys())
# split_files = set(split_files + values)
split_files = values
logging.info(f'splits: {split_files}')

for i, split_file in enumerate(split_files):
    split = file_value_map.get(split_file, split_file)
    out_file = out_dir / f"value~{split_file}.zarr"
    
    # split anndata
    logging.info(f'Split by {split_key}={split}')
    adata_sub = adata[adata.obs[split_key] == split]
    logging.info(adata_sub.__str__())

    # add wildcards
    add_wildcards(adata_sub, {'key': split_key, 'value': split} , 'split_data')
    
    if adata_sub.n_obs == 0:
        write_copy = True
        adata_sub = ad.AnnData(
            X=np.zeros(adata_sub.shape),
            obs=adata_sub.obs,
            obsm=adata_sub.obsm,
            obsp=adata_sub.obsp,
            var=adata_sub.var,
            varm=adata_sub.varm,
            varp=adata_sub.varp,
            layers={
                layer: np.zeros(adata_sub.shape)
                for layer in adata_sub.layers
            },
            uns=adata_sub.uns,
        )
    
    if write_copy:
        adata_sub = dask_compute(adata_sub.copy())
        logging.info(f'Write to {out_file}...')
        write_zarr(adata_sub, out_file)
    else:
        logging.info(f'Write to {out_file} with subset mask...')
        obs_mask = adata.obs_names.isin(adata_sub.obs_names)
        var_mask = adata.var_names.isin(adata_sub.var_names)
        write_zarr_linked(
            adata_sub,
            in_dir=input_file,
            out_dir=out_file,
            subset_mask=(obs_mask, var_mask),
            files_to_keep=['uns']+exclude_slots
        )
    del adata_sub
    logging.info(f'Finished {i+1} out of {len(split_files)}.')

logging.info(f'Finished splitting data by {split_key}.')
# touch done file
Path(snakemake.output.done).touch()