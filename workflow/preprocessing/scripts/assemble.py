"""
Assemble multiple preprocessing outputs into single file
"""
from pathlib import Path
from pprint import pformat, pprint
import numpy as np
import logging
logging.basicConfig(level=logging.INFO)

from utils.io import read_anndata, write_zarr_linked
from utils.annotate import add_wildcards


def assemble_adata(file, file_type, adata, backed=True):

    def deep_update(d, u):
        """Update metadata from uns
        adapted from:
        https://stackoverflow.com/questions/3232943/update-value-of-a-nested-dictionary-of-varying-depth#3233356
        """
        for k, v in u.items():
            d[k] = deep_update(d.get(k, {}), v) if isinstance(v, dict) else v
        return d

    if adata.n_obs == 0 or adata.n_vars == 0:
        return adata

    if file_type == 'normalize':
        logging.info('add normalised counts')
        adata_pp = read_anndata(file, X='X', backed=backed)
        adata.layers['normcounts'] = adata_pp[:, adata.var_names].X
        adata.X = adata.layers['normcounts']
    elif file_type == 'highly_variable_genes':
        logging.info('add highly variable genes')
        adata_pp = read_anndata(file, var='var')
        # subset if adata has been subsetted by HVG
        if any(adata.var_names != adata_pp.var_names):
            adata = adata[:, adata_pp.var_names]
        hvg_column_map = {
            'highly_variable': False,
            'means': 0,
            'dispersions': 0,
            'dispersions_norm': 0,
            'highly_variable_nbatches': 0,
            'highly_variable_intersection': False,
        }
        hvg_columns = [
            column
            for column in hvg_column_map
            if column in adata_pp.var.columns
        ]
        adata.var[hvg_columns] = adata_pp.var[hvg_columns]
        adata.var = adata.var.fillna(hvg_column_map)
    elif file_type == 'extra_hvgs':
        var = read_anndata(file, var='var').var
        hvg_col = [col for col in var.columns if col.startswith('extra_hvgs')][0]
        adata.var[hvg_col] = var[hvg_col]
    elif file_type == 'pca':
        logging.info('add PCA')
        adata_pp = read_anndata(file, obsm='obsm', uns='uns', varm='varm')
        adata.obsm['X_pca'] = adata_pp.obsm['X_pca']
        adata.uns['pca'] = adata_pp.uns['pca']
        adata.varm = adata_pp.varm
    elif file_type == 'neighbors':
        logging.info('add neighbors')
        adata_pp = read_anndata(file, obsp='obsp', uns='uns')
        adata.uns['neighbors'] = adata_pp.uns['neighbors']
        adata.obsp['distances'] = adata_pp.obsp['distances']
        adata.obsp['connectivities'] = adata_pp.obsp['connectivities']
    elif file_type == 'umap':
        logging.info('add UMAP')
        adata_pp = read_anndata(file, obsm='obsm')
        adata.obsm['X_umap'] = adata_pp.obsm['X_umap']
    else:
        ValueError(f'Unknown file type {file_type}')
    
    adata.uns = deep_update(adata.uns, read_anndata(file, uns='uns').uns)
    return adata


def assemble_zarr(file, file_type, slot_map, in_dir_map):
    update_slot_map = {}

    if file_type == 'normalize':
        logging.info('add normalised counts')
        update_slot_map = {
            'X': 'X',
            'layers': 'layers',
            'raw': 'raw',
            'uns/log1p': 'uns/log1p',
            'uns/preprocessing/log-transformed': 'uns/preprocessing/log-transformed',
            'uns/preprocessing/normalization': 'uns/preprocessing/normalization',
        }
    elif file_type == 'highly_variable_genes':
        logging.info('add highly variable genes')
        update_slot_map = {
            'var/highly_variable': 'var/highly_variable',
            'uns/preprocessing/highly_variable_genes': 'uns/preprocessing/highly_variable_genes',
        }
    elif file_type == 'extra_hvgs':
        logging.info('add highly variable genes')
        var = read_anndata(file, var='var').var
        hvg_col = [col for col in var.columns if col.startswith('extra_hvgs')][0]
        update_slot_map = {
            f'var/{hvg_col}': f'var/{hvg_col}',
            f'uns/preprocessing/{hvg_col}': f'uns/preprocessing/{hvg_col}',
        }
    elif file_type == 'pca':
        logging.info('add PCA')
        update_slot_map = {
            'obsm/X_pca': 'obsm/X_pca',
            'uns/pca': 'uns/pca',
            # 'varm/PCs': 'varm/PCs',
        }
    elif file_type == 'neighbors':
        logging.info('add neighbors')
        update_slot_map = {
            'obsp/connectivities': 'obsp/connectivities',
            'obsp/distances': 'obsp/distances',
            'uns/neighbors': 'uns/neighbors',
        }
    elif file_type == 'umap':
        logging.info('add UMAP')
        update_slot_map = {
            'obsm/X_umap': 'obsm/X_umap'
        } 
    else:
        ValueError(f'Unknown file type {file_type}')
    slot_map |= update_slot_map
    in_dir_map |= {slot: file for slot in update_slot_map.values()}
    return slot_map, in_dir_map


output_file = Path(snakemake.output[0])

adata = None
slot_map = {}
in_dir_map = {}
backed=True

file_map = []
for file_type, file in snakemake.input.items():
    if isinstance(file, str):
        file_map.append((file_type, file))
    else:
        file_map.extend((file_type, f) for f in file)

for file_type, file in file_map:
    if adata is None: # read first file
        logging.info(f'Read first file {file}...')
        adata = read_anndata(file, obs='obs', var='var')
        # from scipy import sparse
        # if adata.X is None:
        #     adata.X = sparse.csr_matrix(adata.shape, dtype=np.int8)
        if adata.n_obs == 0:
            logging.info('No data, write empty file...')
            adata.write_zarr(output_file)
            exit(0)

    if file.endswith('.h5ad'):
        adata = assemble_adata(
            file=file,
            file_type=file_type,
            adata=adata,
            backed=backed
        )
    elif file.endswith('.zarr'):
        # if file_type in ['counts']: #, 'highly_variable_genes']:
        #     adata = assemble_adata(file, file_type, adata, backed=backed)
        slot_map, in_dir_map = assemble_zarr(
            file=file,
            file_type=file_type,
            slot_map=slot_map,
            in_dir_map=in_dir_map,
        )
    else:
        ValueError(f'Unknown file type {file}')


logging.info(f'Write to {output_file}...')
adata.uns = read_anndata(file, uns='uns').uns
add_wildcards(adata, snakemake.wildcards, 'preprocessing')
write_zarr_linked(
    adata=adata,
    in_dir=file_map[0][1],
    out_dir=output_file,
    files_to_keep=['uns/wildcards'],
    slot_map=slot_map,
    in_dir_map=in_dir_map,
)