import logging
logging.basicConfig(level=logging.INFO)
from tqdm import tqdm
from tqdm.dask import TqdmCallback
from pathlib import Path
import pandas as pd
from dask import dataframe as dd
from pprint import pformat
from dask import config as da_config

da_config.set(num_workers=snakemake.threads)
logging.info(f"Dask using {da_config.get('num_workers')} workers")


from utils.io import read_anndata, write_zarr_linked
from utils.misc import dask_compute


def convert_dtypes(col):
    # check if can be cast to boolean
    if set(col.dropna().str.lower()) <= {"true", "false"}:
        bool_map = {"true": True, "false": False}
        return col.str.lower().map(bool_map).astype('boolean')
    try:
        return pd.to_numeric(col)
    except Exception:
        return col


def read_table(file, index_name, obs_names, usecols):
    kwargs = dict(
        usecols=usecols,
        comment='#',
        na_values=['nan', 'NaN', 'NA', 'na', ''],
        keep_default_na=True,
        engine='c',
        dtype=str,
        low_memory=False,
        blocksize='200MB',
    )
    
    # define file reader
    if file.endswith('.csv'):
        read_func = dd.read_csv
    elif file.endswith('.tsv'):
        read_func = dd.read_table
    elif file.endswith('.parquet'):
        read_func = pd.read_parquet
        kwargs = {}
    else:
        raise ValueError(f'Unsupported file type: {file}')
    
    df = read_func(file, **kwargs)
    
    # load to memory
    if isinstance(df, dd.DataFrame):
        if index_name in df.columns:
            # subset to entries in adata.obs
            overlapping_idx = df[index_name].isin(obs_names)
            df = df.loc[overlapping_idx, :].copy()
        with TqdmCallback(desc=f'Loading {file}', total=df.shape[0], miniters=1, leave=True):
            df = df.compute()
    
    if index_name in df.columns:
        df = df.set_index(index_name)
    if df.index.name == index_name and index_name is not None:
        # make sure shape matches adata.obs
        index_old = df.index
        df = df.reindex(obs_names)
    
        # check that values aren't all NaN
        if df.empty or df.isna().all().all():
            raise ValueError(
                'No matching data found in mapping file.\n'
                f'Mapping file index:\n{index_old}\nobs_names:\n{obs_names}\n'
                f'Mapping file "{file}" might be mismatched with adata.obs'
            )
    
    return df.apply(convert_dtypes)
    

input_file = snakemake.input.anndata
input_new_cols = snakemake.input.get('new_columns')
input_merge_cols = snakemake.input.get('merge_columns')
output_file = snakemake.output.zarr
file_id = snakemake.wildcards.file_id
rename_columns = snakemake.params.get('rename_columns')
selective_update = snakemake.params.get('selective_update')

logging.info('Read adata...')
kwargs = dict(dask=True, backed=True)
if input_file.endswith(('.zarr', '.zarr/')):
    kwargs |= dict(obs='obs')
adata = read_anndata(input_file, **kwargs)

# merge new columns
if rename_columns:
    logging.info(f'Rename columns:\n{pformat(rename_columns)}')
    for old_col, new_col in tqdm(rename_columns.items(), desc='Renaming'):
        if old_col not in adata.obs.columns:
            adata.obs[old_col] = 'nan'
        else:
            adata.obs[new_col] = adata.obs[old_col]

if input_new_cols:
    index_col = snakemake.params.get('index_col', 'index')
    mapping_order = snakemake.params.get('mapping_order')
    old_index_name = adata.obs.index.name or 'index'
    
    if index_col in adata.obs.columns:
        adata.obs[index_col] = adata.obs[index_col].astype(str)
        adata.obs = adata.obs.reset_index().set_index(index_col)
    else:
        adata.obs_names.name = index_col
    
    label_mapping = read_table(
        file=input_new_cols,
        index_name=index_col,
        obs_names=adata.obs_names,
        usecols=mapping_order
    )    
    logging.info(f'\n{pformat(label_mapping)}')
    
    logging.info(f'Mapping order: {pformat(mapping_order)}')
    label_key = None
    
    for mapping_label in tqdm(mapping_order, desc='Mapping new columns'):
        
        if label_key is None:
            assert mapping_label in adata.obs.columns.tolist() + [index_col], \
                f'"{mapping_label}" not found in adata.obs. ' \
                'Please make sure the first entry in the mapping order is a column in adata.obs.'
            label_key = mapping_label
            continue        
        
        if label_key == index_col:
            adata.obs[mapping_label] = label_mapping[mapping_label]
            # Note: keep label_key as index
            continue
        
        # get unique mapping
        df = label_mapping[[label_key, mapping_label]].drop_duplicates()
        
        # remove trailing whitespaces
        remove_trailing_whitespaces = lambda x: x.strip() if isinstance(x, str) else x
        adata.obs[label_key] = adata.obs[label_key].apply(remove_trailing_whitespaces)
        df[label_key] = df[label_key].apply(remove_trailing_whitespaces)
        df[mapping_label] = df[mapping_label].apply(remove_trailing_whitespaces)
        
        # apply mapping
        map_dict = df.set_index(label_key)[mapping_label].to_dict()
        adata.obs[mapping_label] = pd.Series(adata.obs[label_key].map(map_dict), dtype="category")
        
        # check mapping success
        n_na = adata.obs[mapping_label].isna().sum()
        if n_na > 0:
            logging.info(f'Number of unmapped entries: {n_na}')

    value_counts = adata.obs[
        [x for x in mapping_order if x != index_col]
    ].value_counts(dropna=False, sort=False)
    # sort by unique values per column
    value_counts = value_counts.sort_index(level=list(range(value_counts.index.nlevels)))
    logging.info(f'New columns:\n{pformat(value_counts)}')

    if old_index_name in adata.obs.columns:
        logging.info(f'Reset index "{old_index_name}"...')
        adata.obs = adata.obs.reset_index().set_index(old_index_name)

# merge existing columns
if input_merge_cols:
    sep = snakemake.params.get('merge_sep', '-')
    merge_config_df = pd.read_table(input_merge_cols, comment='#')
    merge_config_df = merge_config_df[merge_config_df['file_id'] == file_id]
    logging.info(f'Merge columns:\n{pformat(merge_config_df)}')
    
    # check if required columns are present
    for col in ['file_id', 'column_name', 'columns']:
        assert col in merge_config_df.columns, \
            f'"{col}" not found in {input_merge_cols}\n{merge_config_df}'    
    assert not merge_config_df.duplicated().any(), \
        f'Duplicated rows in {input_merge_cols} for {file_id}\n{merge_config_df[merge_config_df.duplicated()]}'
    
    for _, row in tqdm(
        merge_config_df.iterrows(),
        total=merge_config_df.shape[0],
        desc=f'Merge values with sep="{sep}"',
    ):
        col_name = row['column_name']
        cols = row['columns'].split(',')
        adata.obs[col_name] = adata.obs[cols].apply(
            lambda x: sep.join(x.values.astype(str)), axis=1
        )

if selective_update:
    base_column = selective_update.get('base_column')
    new_column = selective_update.get('new_column', base_column)
    update_map = selective_update.get('update_map', {})
    assert isinstance(update_map, dict)
        
    adata.obs[new_column] = adata.obs[base_column].astype(str)
    
    for col, _dict in tqdm(update_map.items(), desc='Selective update'):
        assert col in adata.obs.columns, f'"{col}" not found in adata.obs'
        
        # parse to dictionary if not already
        if not isinstance(_dict, dict):
            df = pd.read_table(_dict, comment='#')
            _dict = df[[col, base_column]].drop_duplicates().set_index(col)[base_column].to_dict()
        
        labels = list(_dict.keys())
        query = f'{col}.isin(@labels)'
        
        # map selective values
        s = adata.obs[col].astype(str).map(_dict)
        adata.obs.loc[~s.isna(), new_column] = s.dropna()
        
        # check values before and after mapping
        before = adata.obs.query(query)[[col, base_column]].value_counts(dropna=False, sort=False)
        after = adata.obs.query(query)[[col, new_column]].value_counts(dropna=False, sort=False)
        # fix this
        if before.to_dict() != after.to_dict():
            logging.info(f'Before mapping:\n{pformat(before)}')
            logging.info(f'After mapping:\n{pformat(after)}')


logging.info(f'Write to {output_file}...')
dask_compute(adata)
write_zarr_linked(
    adata,
    in_dir=input_file,
    out_dir=output_file,
    files_to_keep=['obs']
)