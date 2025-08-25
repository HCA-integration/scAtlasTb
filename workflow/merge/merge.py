import logging
from pathlib import Path
logging.basicConfig(level=logging.INFO)
import gc
import faulthandler
faulthandler.enable()
from scipy.sparse import issparse
import tqdm.dask as tdask
from tqdm import tqdm
import pandas as pd
import yaml
import scanpy as sc
from anndata.experimental import AnnCollection
from anndata import AnnData
from dask import array as da
from dask import config as da_config
# from dask.diagnostics import ProgressBar

da_config.set(
    **{
        'num_workers': snakemake.threads,
        'array.slicing.split_large_chunks': False
    }
)
logging.info(f"Dask using {da_config.get('num_workers')} workers")

from utils.io import read_anndata, link_zarr, write_zarr
from utils.misc import apply_layers, dask_compute, ensure_sparse


def read_adata(
    file,
    file_id=None,
    backed=False,
    dask=False,
    stride=10_000,
    chunks=(-1, -1),
    log=True,
    **slots
):
    if file_id is None:
        file_id = file
    if log:
        logging.info(f'Read {file}...')
    try:
        adata = read_anndata(
            file,
            backed=backed,
            dask=dask,
            stride=stride,
            chunks=chunks,
            **slots,
            verbose=False,
        )
        adata.obs['file_id'] = file_id
        adata = ensure_sparse(adata)
    except Exception as e:
        raise Exception(f'Error reading {file}') from e
    
    return adata


def remove_slots(adata):
    for slot in ['X', 'layers', 'obsm', 'obsp', 'varm', 'varp']:
        if hasattr(adata, slot):
            delattr(adata, slot)


dataset = snakemake.wildcards.dataset
files = snakemake.input
out_file = snakemake.output.zarr

merge_strategy = snakemake.params.get('merge_strategy', 'inner')
keep_all_columns = snakemake.params.get('keep_all_columns', False)
backed = snakemake.params.get('backed', False)
dask = snakemake.params.get('dask', False)
stride = snakemake.params.get('stride', 500_000)
slots = snakemake.params.get('slots', {})
if slots is None:
    slots = {}
    
logging.info(f'Read slots: {slots}')

if len(files) == 1:
    link_zarr(in_dir=files[0], out_dir=out_file)
    exit(0)

# subset to non-empty datasets
def check_cells(file):
    if file.endswith('.zarr'):
        zattr_path = Path(file) / slots.get('X', 'X') / '.zattrs'
        with open(zattr_path, 'r') as f:
            zattrs = yaml.safe_load(f)
        return 'shape' in zattrs
        # return zattrs['shape'][0] > 0
    return read_anndata(file, obs='obs', verbose=False).n_obs > 0
files = {
    file_id: file
    for file_id, file
    in zip(files.keys(), files)
    if check_cells(file)
}

if len(files) == 0:
    logging.info('All adatas are empty, skip concatenation...')
    AnnData().write_zarr(out_file)
    exit(0)

adatas = []

if dask:
    logging.info('Read all files with dask...')
    for file_id, file_path in tqdm(files.items(), desc='Read files', miniters=1):
        _ad = read_adata(
            file_path,
            file_id=file_id,
            backed=backed,
            dask=dask,
            stride=stride,
            **slots,
            log=False,
        )
        # logging.info(f'{file_id} shape: {_ad.shape}')
        
        #with tdask.TqdmCallback(desc='Persist'):
            # _ad = apply_layers(_ad, func=lambda x: x.persist())
        
        adatas.append(_ad)
    
    # concatenate
    adata = sc.concat(adatas, join=merge_strategy)
    print(adata, flush=True)

    for _ad in adatas:
        remove_slots(_ad)
        gc.collect()
    
    #if backed:
    #    with tdask.TqdmCallback(desc='Persist'): # ProgressBar():
    #        adata = apply_layers(adata, func=lambda x: x.rechunk((stride, -1)).persist() if isinstance(x, da.Array) else x)
    #         adata = apply_layers(adata, func=lambda x: x.persist() if isinstance(x, da.Array) else x)
    
    # with tdask.TqdmCallback(desc='Rechunk'):
    #     adata  = apply_layers(adata, func=lambda x: x.rechunk((stride, -1)))
    
elif backed:
    logging.info('Read all files in backed mode...')
    adatas = [
        read_adata(file_path, file_id=file_id, backed=backed, dask=dask)
        for file_id, file_path in files.items()
    ]
    dc = AnnCollection(
        adatas,
        join_obs='outer',
        join_obsm=None,
        join_vars=merge_strategy,
        indices_strict=not backed,
    )
    logging.info(dc.__str__())

    logging.info('Subset AnnDataCollection, returning a View...')
    adata = dc[:].to_adata()
    assert adata.X is not None
    
    for _ad in adatas:
        remove_slots(_ad)
        gc.collect()

else:

    adata = None
    for file_id, file_path in tqdm(files.items()):
        logging.info(f'Read {file_path}...')
        _ad = read_adata(
            file_path,
            file_id=file_id,
            **slots,
            backed=backed,
            dask=dask
        )
        logging.info(f'{file_id} shape: {_ad.shape}')
        
        if adata is None:
            adata = _ad
            continue
        
        # merge adata
        adata = sc.concat([adata, _ad], join=merge_strategy)
        logging.info(f'after merge:\n{adata}')

        remove_slots(_ad)
        adatas.append(_ad)
        gc.collect()


# merge lost annotations
logging.info('Add gene info...')
for _ad in adatas:
    # get intersection of var_names
    var_names = list(set(adata.var_names).intersection(set(_ad.var_names)))
    
    # conver categorical columns to str
    categorical_columns = _ad.var.select_dtypes(include='category').columns
    _ad.var[categorical_columns] = _ad.var[categorical_columns].astype(str)
    
    # add new columns to adata.var
    adata.var.loc[var_names, _ad.var.columns] = _ad.var.loc[var_names, :]

# fix dtypes
adata.var = adata.var.infer_objects()
logging.info(adata.var)

# set new indices
adata.obs[f'obs_names_before_{dataset}'] = adata.obs_names
adata.obs_names = dataset + '-' + adata.obs.reset_index(drop=True).index.astype(str)

# add uns data
adata.uns['dataset'] = dataset
logging.info(adata.__str__())

logging.info(f'Write to {out_file}...')
with tdask.TqdmCallback():
    if isinstance(adata.X, da.Array):
        print(adata.X.dask, flush=True)
    write_zarr(adata, out_file)
