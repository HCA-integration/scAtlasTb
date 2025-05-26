import logging
logging.basicConfig(level=logging.INFO)
import gc
import faulthandler
faulthandler.enable()
import pandas as pd
import scanpy as sc
from anndata.experimental import AnnCollection
from anndata import AnnData

from load_data_utils import SCHEMAS, get_union
from utils.io import read_anndata, link_zarr
from utils.misc import apply_layers


def read_adata(file, keep_columns, backed=False, dask=False):
    logging.info(f'Read {file}...')
    adata = read_anndata(
        file,
        backed=backed,
        dask=dask,
        X='X',
        obs='obs',
        var='var',
        uns='uns',
    )
    if not keep_columns:
        logging.info(f'Keep only mandatory columns for {file}...')
        adata.obs = adata.obs[get_union(SCHEMAS["CELLxGENE_OBS"], SCHEMAS["EXTRA_COLUMNS"])]
    adata.var = adata.var[SCHEMAS['CELLxGENE_VARS']]
    return adata


dataset = snakemake.params.dataset
files = snakemake.input
out_file = snakemake.output.zarr

merge_strategy = snakemake.params.get('merge_strategy', 'inner')
keep_all_columns = snakemake.params.get('keep_all_columns', False)
backed = snakemake.params.get('backed', False)
dask = snakemake.params.get('dask', False)

if len(files) == 1:
    link_zarr(in_dir=files[0], out_dir=out_file)
    exit(0)

adatas = [read_anndata(file, obs='obs', var='var', uns='uns') for file in files]
adatas = [adata for adata in adatas if adata.n_obs > 0]

if len(adatas) == 0:
    logging.info('All adatas are empty, skip concatenation...')
    AnnData().write_zarr(out_file)
    exit(0)

if dask:
    from dask import array as da
    from dask import config as da_config

    da_config.set(
        **{
            'num_workers': snakemake.threads,
            'array.slicing.split_large_chunks': True
        }
    )
    logging.info('Read all files with dask...')
    logging.info(f'n_threads: {snakemake.threads}')
    adatas = [read_adata(file, keep_all_columns, backed, dask) for file in files]
    adatas = [adata for adata in adatas if adata.n_obs > 0]
    
    # concatenate
    adata = sc.concat(adatas, join=merge_strategy)
    
elif backed:
    logging.info('Read all files in backed mode...')
    adatas = [read_adata(file, keep_all_columns, backed, dask) for file in files]
    adatas = [adata for adata in adatas if adata.n_obs > 0]
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

else:
    logging.info(f'Read first file {files[0]}...')
    adata = read_adata(files[0], keep_all_columns, backed=backed, dask=dask)
    logging.info(adata.__str__())

    uns_per_dataset = {
        adata.uns['meta']['dataset']: adata.uns['meta']
    }

    for file in files[1:]:
        logging.info(f'Read {file}...')
        _adata = read_adata(file, keep_all_columns, backed=backed, dask=dask)
        logging.info(_adata.__str__())
        
        # merge adata
        adata = sc.concat([adata, _adata], join=merge_strategy)

        del _adata
        gc.collect()

# merge lost annotations
var_map = adatas[0].var
uns = adatas[0].uns
uns_per_dataset = {
    uns['meta']['dataset']: uns['meta']
}
for _adata in adatas[1:]:
    # uns
    uns_per_dataset[_adata.uns['meta']['dataset']] = _adata.uns['meta']
    
    # genes
    var_map = pd.merge(
        var_map,
        _adata.var,
        how=merge_strategy,
        on=['feature_id'] + SCHEMAS['CELLxGENE_VARS']
    )
    var_map = var_map[~var_map.index.duplicated()]

logging.info('Add gene info...')
adata = adata[:, var_map.index]
adata.var = var_map.loc[adata.var_names]
logging.info(adata.var)
assert 'feature_name' in adata.var

# add uns data
organ = adata.obs['organ'].unique()
assert len(organ) == 1
adata.uns['dataset'] = dataset
adata.uns['organ'] = organ
adata.uns['meta'] = {
    'dataset': dataset,
    'organ': organ,
    'per_dataset': uns_per_dataset,
}
logging.info(adata.__str__())

logging.info(f'Write to {out_file}...')
adata.write_zarr(out_file)