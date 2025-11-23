import logging
logging.basicConfig(level=logging.INFO)
import scanpy as sc
import anndata as ad
import pandas as pd
from pprint import pformat
from tqdm import tqdm
from scipy.sparse import csr_matrix
import re

from utils.io import read_anndata, get_file_reader, write_zarr_linked, link_file
from utils.misc import dask_compute
from collect_utils import get_same_columns, merge_df


dataset = snakemake.wildcards.dataset
files = snakemake.input
output_file = snakemake.output.zarr
sep = snakemake.params.get('sep', '--')
obs_index_col_map = snakemake.params.get('obs_index_col', {})
same_slots = snakemake.params.get('same_slots', [])
merge_slots = snakemake.params.get('merge_slots', [])
skip_slots = snakemake.params.get('skip_slots', [])
kwargs = {k: k for k in same_slots+merge_slots if k not in skip_slots}

LARGE_SLOTS = ['layers', 'obsm', 'obsp']

# parse obs index
if isinstance(obs_index_col_map, str):
    obs_index_col_map = {file_id: obs_index_col_map for file_id in files.keys()}
assert isinstance(obs_index_col_map, dict), 'obs_index_col_map must be a dict'

# parse regex keys
obs_index_col_map = {
    file_id: next(
        (v for k, v in obs_index_col_map.items() if re.fullmatch(k, file_id)), 
        None  # Default if no match is found
    )
    for file_id in files.keys()
}

if len(files) == 1:
    logging.info('Single file, write unmodified...')
    adata = read_anndata(files[0], **kwargs, dask=True, backed=True)
    print(adata, flush=True)
    if files[0].endswith('.zarr'):
        write_zarr_linked(adata, in_dir=files[0], out_dir=output_file)
    else:
        adata.write_zarr(output_file)
    exit(0)


def read_file(file, **kwargs):
    if file.endswith('.zarr'):
        kwargs = {k: v for k, v in kwargs.items() if k not in LARGE_SLOTS+['X', 'raw']}
        
        adata = read_anndata(file, verbose=False, **kwargs)
        func, _ = get_file_reader(file)
        
        default_values = [
            csr_matrix(adata.shape), # for layers
            csr_matrix((adata.n_obs, 0)), # for obsm
            csr_matrix((adata.n_obs, adata.n_obs)), # for obsp
        ]
        for slot_name, default_value in zip(LARGE_SLOTS, default_values):
            slot_keys = func(file, 'r').get(slot_name, {}).keys()
            setattr(
                adata,
                f'_{slot_name}',
                {key: default_value for key in slot_keys}
            )
        
        return adata
    return read_anndata(file, verbose=False, **kwargs)


adatas = {
    file_id: read_file(file, **kwargs, dask=True, backed=True)
    for file_id, file in tqdm(
        files.items(),
        desc='Read files',
        miniters=1,
        total=len(files),
    )
}
logging.info(f'File ids: \n{pformat(list(adatas.keys()))}')

# Ensure all adatas have the same obs_names and var_names, and match their order to the first adata
adata_iter = iter(adatas.items())
first_key, first_adata = next(adata_iter)
first_obs_names = first_adata.obs_names
first_var_names = first_adata.var_names
missing_slots = {slot: slot for slot in LARGE_SLOTS if slot in merge_slots}
write_copy_for = {slot: [] for slot in missing_slots}

for file_id, _ad in adata_iter:
    # get file path
    file_name = files[file_id]
    
    # Check obs_names
    if not _ad.obs_names.equals(first_obs_names):
        # read slots that need to be reordered
        _ad_tmp = read_anndata(
            file_name,
            **missing_slots,
            dask=True,
            backed=True,
            dask_slots=['layers', 'obsm', 'obsp'],
            verbose=False,
        )
        for slot_name in missing_slots:
            setattr(_ad, f'_{slot_name}', getattr(_ad_tmp, f'_{slot_name}'))
        
        # set write status for affected slots
        write_copy_for |= {slot: write_copy_for[slot] + [file_id] for slot in missing_slots}
        
        # reorder obs_names
        if set(_ad.obs_names) == set(first_obs_names):
            logging.info(f"Reordering obs_names for file_id={file_id} to match the first adata.")
            adatas[file_id] = _ad[first_obs_names, :]
        else:
            logging.error(f'file_id{first_key}:\n{adatas[file_id].obs_names}')
            logging.error(f'file_id={file_id}:\n{first_obs_names}')
            raise ValueError(f"obs_names do not match between first adata '{first_key}' and file_id='{file_id}'")
    
    # Check var_names
    if not _ad.var_names.equals(first_var_names):
        # read slots that need to be reordered
        if 'layers' in merge_slots:
            _ad.layers = read_anndata(
                file_name,
                layers='layers',
                dask=True,
                backed=True,
            ).layers

        # set write status for affected slots
        write_copy_for['layers'] = write_copy_for['layers'] + [file_id]

        # reorder var_names
        if set(_ad.var_names) == set(first_var_names):
            logging.info(f"Reordering var_names for file_id={file_id} to match the first adata.")
            adatas[file_id] = adatas[file_id][:, first_var_names]
        else:
            raise ValueError(f"var_names do not match between first adata and file_id={file_id}")
    
    if adatas[file_id].is_view:
        adatas[file_id] = adatas[file_id].copy()

if 'obs' in merge_slots:
    for file_id, _ad in adatas.items():
        obs_index_column = obs_index_col_map.get(file_id)
        if obs_index_column is not None:
            logging.info(f'Set obs index "{obs_index_column}" for file_id={file_id}...')
            assert obs_index_column in _ad.obs.columns, \
                f'Index column "{obs_index_column}" not found for {file_id}\n{_ad.obs}'
            adatas[file_id].obs = _ad.obs.set_index(obs_index_column)
    logging.info('Determine which columns are the same...')
    same_obs_columns = get_same_columns(adatas, n_threads=snakemake.threads)
    logging.info(f'Same columns:\n{pformat(same_obs_columns)}')
else:
    same_obs_columns = []

# TODO: link slots with new slot names for merge_slots

slots = dict()
slot_link_map = dict()
in_dir_map = dict()
file_to_link = None

for file_id, _ad in adatas.items():

    file_name = files[file_id]
    # intialize file_to_link
    if not file_to_link and file_name.endswith('.zarr'):
        file_to_link = file_name

    for slot_name in merge_slots:
        logging.info(f'Merge slot "{slot_name}" for file_id={file_id}...')
        slot = _ad.__dict__.get(f'_{slot_name}')
        update_slot_link_map = dict()
        to_link = file_name.endswith('.zarr') and file_id not in write_copy_for.get(slot_name, [])
        
        if isinstance(slot, pd.DataFrame):
            slots[slot_name] = merge_df(
                df_current=slot,
                file_id=file_id,
                df_previous=slots.get(slot_name),
                same_columns=same_obs_columns,
                sep=sep,
            )
        elif slot_name == 'X':
            if to_link:
                update_slot_link_map[f'layers/{slot_name}{sep}{file_id}'] = f'{slot_name}'
            else:
                new_slot = {f'{slot_name}{sep}{file_id}': slot}
                slots['layers'] = slots.get('layers', {}) | new_slot
        elif hasattr(slot, 'items'):
            if to_link:
                update_slot_link_map = {
                    f'{slot_name}/{key}{sep}{file_id}': f'{slot_name}/{key}'
                    for key in slot.keys()
                }
            else:
                new_slot = {
                    f'{key}{sep}{file_id}': value
                    for key, value in slot.items()
                }
                slots[slot_name] = slots.get(slot_name, {}) | new_slot
        else:
            raise NotImplementedError(f'Slot "{slot_name}" not supported')
        
        slot_link_map |= update_slot_link_map
        in_dir_map |= {
            key: file_name for key in update_slot_link_map.keys()
        }

# deal with same slots
if file_to_link:
    slot_link_map |= {slot: slot for slot in same_slots}
else:
    for slot_name in same_slots:
        slots[slot_name] = _ad.__dict__.get(f'_{slot_name}')

logging.info('Create AnnData object...')
adata = ad.AnnData(**slots)
print(adata, flush=True)

logging.info('Compute matrix...')
compute_layers = [v for slot in LARGE_SLOTS for v in getattr(adata, f'_{slot}').keys() if slot ]
print(f'Compute layers: {compute_layers}', flush=True)
adata = dask_compute(adata, layers=compute_layers)

logging.info(f'Write to {output_file}...')
write_zarr_linked(
    adata=adata,
    in_dir=file_to_link,
    out_dir=output_file,
    files_to_keep=merge_slots+skip_slots,
    # slot_map=slot_link_map,
    # in_dir_map=in_dir_map,
)

# workaround hack - need to fix bug in write_zarr_linked
print('in_dir_map\n', pformat(in_dir_map), flush=True)
for slot, file in in_dir_map.items():
    link_file(in_file=f'{file}/{slot_link_map[slot]}', out_file=f'{output_file}/{slot}')
