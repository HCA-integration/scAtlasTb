import warnings
import glob
import yaml
from pprint import pformat
import zarr
import logging
from itertools import product
logging.basicConfig(level=logging.INFO)

from utils.misc import check_sparse_equal
from utils.io import read_anndata
from PreprocessingConfig import HVG_PARAMS as HVG_KEYS
def _build_suffixes(param_dict):
    """Return list of sorted key=value suffix strings for all combinations."""
    if not isinstance(param_dict, dict) or len(param_dict) == 0:
        return []
    filtered = {k: v for k, v in param_dict.items() if k in HVG_KEYS}
    if not filtered:
        return []
    
    # Generate all combinations of parameter values
    keys = sorted(filtered.keys())
    values = [filtered[k] if isinstance(filtered[k], list) else [filtered[k]] for k in keys]
    
    suffixes = []
    for combo_values in product(*values):
        combo = dict(zip(keys, combo_values))
        suffix = '--' + '--'.join(f"{k}={combo[k]}" for k in keys)
        suffixes.append(suffix)
    return suffixes


def expected_hvg_columns(hvg_config):
    """Expected 'highly_variable--...' columns from HVG config."""
    return ['highly_variable' + s for s in _build_suffixes(hvg_config)]


def expected_extra_hvg_columns(extra_config):
    """Expected 'extra_hvgs--...' columns from extra_hvgs overwrite_args.
    If no overwrite_args provided, expect plain 'extra_hvgs'.
    """
    if not isinstance(extra_config, dict):
        return []
    overwrite = extra_config.get('overwrite_args', {})
    suffixes = _build_suffixes(overwrite)
    if suffixes:
        return ['extra_hvgs' + s for s in suffixes]
    return ['extra_hvgs']

config = yaml.safe_load(open('test/config.yaml', 'r'))
config_no_gpu = yaml.safe_load(open('test/config_no_gpu.yaml', 'r'))
config['DATASETS'] |= config_no_gpu['DATASETS']

single_outputs = glob.glob('test/out/preprocessing/dataset~*/file_id~*.zarr') \
    + glob.glob('test/out/preprocessing/dataset~no_gpu/file_id~*.zarr')
assert len(single_outputs) > 0, 'No output files found'

for file in single_outputs:
    matching_datasets = [x for x in config['DATASETS'] if x in file]
    if not matching_datasets or matching_datasets[0] == 'empty':
        logging.info(f'No matching dataset found for {file}, skipping...')
        continue
    
    dataset = matching_datasets[0]
    logging.info(f'Check dataset "{dataset}", file {file}')

    preprocessing_config = config['DATASETS'][dataset]['preprocessing']

    z = zarr.open(file)
    
    if 'X' in z:
        X = read_anndata(file, X='X', dask=True, backed=True, verbose=False).X    

    try:
        if 'normcounts' in preprocessing_config['assemble']:
            assert 'raw' in z, list(z)
            adata_raw = read_anndata(file, raw='raw', dask=True, backed=True, verbose=False)
            raw = adata_raw.raw
            counts = read_anndata(file, X='layers/counts', dask=True, backed=True, verbose=False).X
            assert counts.min() == 0, f'counts.min() = {counts.min()}'
            assert counts.max() > 100, f'counts.max() = {counts.max()}'
            assert not isinstance(raw, type(None))
            if counts.shape[1] == raw.n_vars:
                raw_range = (raw.X.min(), raw.X.max())
                counts_range = (counts.min(), counts.max())
                message = f"matrices aren't equal\nlayers/counts range: {counts_range} vs. raw.X range: {raw_range}"
                assert check_sparse_equal(counts, raw.X), message
                assert any(X.data != raw.X.data)
            
            # check if values are log-transformed
            assert X.min() == 0, f'X.min() = {X.min()}'
            assert X.max() < 20, f'X.max() = {X.max()}'
            
            assert 'normcounts' in z['layers']
            adata_norm = read_anndata(file, X='layers/normcounts', dask=True, backed=True, verbose=False)
            normcounts = adata_norm.X
            assert check_sparse_equal(X, normcounts)

        if 'highly_variable_genes' in preprocessing_config['assemble']:
            var = read_anndata(file, var='var', verbose=False).var
            assert 'highly_variable' in var.columns
            # Build and assert expected HVG variant columns from config
            for col in expected_hvg_columns(preprocessing_config.get('highly_variable_genes')):
                assert col in var.columns, (
                    f'Expected HVG column "{col}" not found. '
                    f'Available HVG columns: {[c for c in var.columns if c.startswith("highly_variable")]}'
                )
                logging.info(f'Found expected HVG column: {col}')
            if preprocessing_config['highly_variable_genes'] is False:
                assert var['highly_variable'].sum() == var.shape[0], \
                    f'Expected all genes to be marked highly_variable, got {var["highly_variable"].sum()} / {var.shape[0]}'

        if 'extra_hvgs' in preprocessing_config['assemble']:
            var = read_anndata(file, var='var', verbose=False).var
            for col in expected_extra_hvg_columns(preprocessing_config.get('extra_hvgs')):
                assert col in var.columns, (
                    f'Expected extra_hvgs column "{col}" not found. '
                    f'Available columns: {[c for c in var.columns if c.startswith("extra_hvgs")]}'
                )
                logging.info(f'Found expected extra_hvgs column: {col}')

        if 'pca' in preprocessing_config['assemble']:
            assert 'X_pca' in z['obsm']

        if 'neighbors' in preprocessing_config['assemble']:
            assert 'neighbors' in z['uns']
            assert 'distances' in z['obsp']
            assert 'connectivities' in z['obsp']
            
        assert 'wildcards' in z['uns']
    except AssertionError as e:
        logging.error(pformat(preprocessing_config))
        logging.error(f'Assertion failed for dataset "{dataset}", file {file}')
        raise e
