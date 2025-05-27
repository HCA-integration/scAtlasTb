import logging
logging.basicConfig(level=logging.INFO)
from tqdm import tqdm
from tqdm.dask import TqdmCallback
from pathlib import Path
import pandas as pd
from pprint import pformat
import anndata

from utils.io import read_anndata, write_zarr_linked, link_zarr_partial

input_file = snakemake.input.anndata
input_new_cols = snakemake.input.get('new_columns')
input_merge_cols = snakemake.input.get('merge_columns')
output_file = snakemake.output.zarr
file_id = snakemake.wildcards.file_id
rename_columns = snakemake.params.get('rename_columns')
selective_update = snakemake.params.get('selective_update')

logging.info('Read adata...')
if input_file.endswith('.zarr'):
    adata = read_anndata(input_file, obs='obs')
else:
    adata = anndata.read_h5ad(input_file)  # read complete file to convert to zarr
if adata.obs_names.name in (None, ''):
    adata.obs_names.name = 'index'

# merge new columns
if rename_columns is not None:
    logging.info(f'Rename columns...')
    # adata.obs = adata.obs.rename(columns=rename_columns)
    
    for old_col, new_col in rename_columns.items():
        if old_col not in adata.obs.columns:
            adata.obs[old_col] = 'nan'
        else:
            adata.obs[new_col] = adata.obs[old_col]

if input_new_cols is not None:
    mapping_order = snakemake.params.get('mapping_order')
    logging.info(f'Mapping order:\n{mapping_order}')
    label_mapping = pd.read_table(input_new_cols, comment='#')
    
    label_key = None
    for mapping_label in mapping_order:
        if label_key is None:
            if mapping_label == adata.obs_names.name:
                adata.obs.reset_index(inplace=True)
                label_key = mapping_label
            assert mapping_label in adata.obs.columns, f'"{mapping_label}" not found in adata.obs.columns. ' \
                'Please make sure the first entry in the mapping order is a column in adata.obs.'
            label_key = mapping_label
            continue
        
        logging.info(f'mapping "{label_key}" to "{mapping_label}"...')
        
        # get unique mapping
        df = label_mapping[[label_key, mapping_label]].drop_duplicates()
        
        # remove trailing whitespaces
        remove_trailing_whitespaces = lambda x: x.strip() if isinstance(x, str) else x
        adata.obs[label_key] = adata.obs[label_key].apply(remove_trailing_whitespaces)
        df[label_key] = df[label_key].apply(remove_trailing_whitespaces)
        df[mapping_label] = df[mapping_label].apply(remove_trailing_whitespaces)
        
        # apply mapping
        map_dict = df.set_index(label_key)[mapping_label].to_dict()
        logging.info(f'Mapping:\n{pformat(map_dict)}')
        adata.obs[mapping_label] = pd.Series(adata.obs[label_key].map(map_dict), dtype="category")
        
        logging.info(f'Number of unmapped entries: {adata.obs[mapping_label].isna().sum()}')
        logging.info(pformat(adata.obs[[mapping_label, label_key]].value_counts(dropna=False, sort=False)))

        # set current mapping label as new label key
        # label_key = mapping_label

# merge existing columns
if input_merge_cols is not None:
    sep = snakemake.params.get('merge_sep', '-')
    logging.info(f'Merge existing columns with sep="{sep}"...')
    
    merge_config_df = pd.read_table(input_merge_cols, comment='#')
    for col in ['file_id', 'column_name', 'columns']:
        assert col in merge_config_df.columns, f'"{col}" not found in {input_merge_cols}\n{merge_config_df}'
    
    merge_config_df = merge_config_df[merge_config_df['file_id'] == file_id]
    assert not merge_config_df.duplicated().any(), f'Duplicated rows in {input_merge_cols} for {file_id}\n{merge_config_df.duplicated()}'
    
    for _, row in merge_config_df.iterrows():
        col_name = row['column_name']
        cols = row['columns'].split(',')
        logging.info(f'merge {col_name} from {cols}...')
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
write_zarr_linked(
    adata,
    in_dir=input_file,
    out_dir=output_file,
    files_to_keep=['obs']
)