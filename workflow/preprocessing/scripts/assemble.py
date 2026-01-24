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


def assemble_adata(file, file_type, adata, is_default=True):

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
        adata_pp = read_anndata(file, X='X', dask=True, backed=True, verbose=False)
        adata.layers['normcounts'] = adata_pp[:, adata.var_names].X
        adata.X = adata.layers['normcounts']
    
    elif file_type == 'highly_variable_genes':
        logging.info('add highly variable genes')
        adata_pp = read_anndata(file, var='var', verbose=False)
        # subset if adata has been subsetted by HVG
        if any(adata.var_names != adata_pp.var_names):
            adata = adata[:, adata_pp.var_names]
        
        # Find the specific HVG column (should be highly_variable or highly_variable--)
        hvg_col = next((col for col in var.columns if col.startswith('highly_variable--')), 'highly_variable')
        adata.var[hvg_col] = adata_pp.var[hvg_col]
        if is_default:
            adata.var['highly_variable'] = adata_pp.var[hvg_col]

    elif file_type == 'extra_hvgs':
        var = read_anndata(file, var='var', verbose=False).var
        hvg_col = [col for col in var.columns if col.startswith('extra_hvgs')][0]
        adata.var[hvg_col] = var[hvg_col]
    
    elif file_type == 'pca':
        logging.info('add PCA')
        adata_pp = read_anndata(file, obsm='obsm', uns='uns', varm='varm', verbose=False)
        adata.obsm['X_pca'] = adata_pp.obsm['X_pca']
        adata.uns['pca'] = adata_pp.uns['pca']
        adata.varm = adata_pp.varm
    
    elif file_type == 'neighbors':
        logging.info('add neighbors')
        adata_pp = read_anndata(file, obsp='obsp', uns='uns', verbose=False)
        adata.uns['neighbors'] = adata_pp.uns['neighbors']
        adata.obsp['distances'] = adata_pp.obsp['distances']
        adata.obsp['connectivities'] = adata_pp.obsp['connectivities']
    
    elif file_type == 'umap':
        logging.info('add UMAP')
        adata_pp = read_anndata(file, obsm='obsm', verbose=False)
        adata.obsm['X_umap'] = adata_pp.obsm['X_umap']
    
    else:
        ValueError(f'Unknown file type {file_type}')
    
    adata.uns = deep_update(adata.uns, read_anndata(file, uns='uns', verbose=False).uns)
    return adata


def assemble_zarr(file, file_type, slot_map, in_dir_map, is_default=True):
    update_slot_map = {}

    if file_type == 'normalize':
        logging.info('add normalised counts')
        update_slot_map |= {
            'X': 'X',
            'layers': 'layers',
            'raw': 'raw',
            'uns/log1p': 'uns/log1p',
            'uns/preprocessing/log-transformed': 'uns/preprocessing/log-transformed',
            'uns/preprocessing/normalization': 'uns/preprocessing/normalization',
        }
    
    elif file_type == 'highly_variable_genes':
        logging.info('add highly variable genes')
        var = read_anndata(file, var='var', verbose=False).var

        # map specific hvg_column
        hvg_col = next((col for col in var.columns if col.startswith('highly_variable--')), 'highly_variable')
        update_slot_map |= {
            f'var/{hvg_col}': f'var/{hvg_col}',
            f'uns/preprocessing/{hvg_col}': f'uns/preprocessing/highly_variable_genes',
        }
        
        # additionally pick default for 'highly_variable' slot
        if is_default:
            update_slot_map |= {
                'var/highly_variable': f'var/{hvg_col}',
                'uns/preprocessing/highly_variable': f'uns/preprocessing/highly_variable_genes',
            }
    
    elif file_type == 'extra_hvgs':
        logging.info('add highly variable genes')
        var = read_anndata(file, var='var', verbose=False).var
        
        hvg_col = [col for col in var.columns if col.startswith('extra_hvgs')][0]
        update_slot_map |= {
            f'var/{hvg_col}': f'var/{hvg_col}',
            f'uns/preprocessing/{hvg_col}': f'uns/preprocessing/{hvg_col}',
        }

        # additionaly pick default for 'extra_hvgs' slot
        if is_default:
            update_slot_map |= {
                'var/extra_hvgs': f'var/{hvg_col}',
                'uns/preprocessing/extra_hvgs': f'uns/preprocessing/{hvg_col}',
            }
    
    elif file_type == 'pca':
        logging.info('add PCA')
        update_slot_map |= {
            'obsm/X_pca': 'obsm/X_pca',
            'uns/pca': 'uns/pca',
            'varm/PCs': 'varm/PCs',
        }
    
    elif file_type == 'neighbors':
        logging.info('add neighbors')
        update_slot_map |= {
            'obsp/connectivities': 'obsp/connectivities',
            'obsp/distances': 'obsp/distances',
            'uns/neighbors': 'uns/neighbors',
        }
    
    elif file_type == 'umap':
        logging.info('add UMAP')
        update_slot_map |= {
            'obsm/X_umap': 'obsm/X_umap'
        } 
    
    else:
        ValueError(f'Unknown file type {file_type}')
    
    slot_map |= update_slot_map
    in_dir_map |= {slot: file for slot in update_slot_map.values()}
    
    return slot_map, in_dir_map


input_files = snakemake.input
output_file = Path(snakemake.output.zarr)

adata = None
slot_map = {}
in_dir_map = {}

# Build list of (file_type, file_path) tuples from snakemake input
file_map = []  # all files
default_files = {}  # default per file_type (first in list or only file)

for file_type, file in input_files.items():
    if isinstance(file, str):
        file_map.append((file_type, file))
        default_files[file_type] = file
    else:
        # add all files
        file_map.extend((file_type, f) for f in file)
        # only first file
        default_files[file_type] = file[0]


# Mapping of file extensions to processing functions
def process_file(file_type, file, adata, slot_map, in_dir_map, is_default=True):
    """Process file based on its extension."""
    
    if file.endswith('.h5ad'):
        return assemble_adata(
            file=file,
            file_type=file_type,
            adata=adata,
            is_default=is_default
        ), slot_map, in_dir_map
    elif file.endswith('.zarr') or file.endswith('.zarr/'):
        slot_map, in_dir_map = assemble_zarr(
            file=file,
            file_type=file_type,
            slot_map=slot_map,
            in_dir_map=in_dir_map,
            is_default=is_default
        )
        return adata, slot_map, in_dir_map
    else:
        raise ValueError(f'Unknown file format: {file}')


for file_type, file in file_map:
    if adata is None: # read first file
        logging.info(f'Read first file {file}...')
        adata = read_anndata(file, obs='obs', var='var', verbose=False)
        # from scipy import sparse
        # if adata.X is None:
        #     adata.X = sparse.csr_matrix(adata.shape, dtype=np.int8)
        if adata.n_obs == 0:
            logging.info('No data, write empty file...')
            write_zarr_linked(adata, in_dir=file, out_dir=output_file)
            exit(0)

    adata, slot_map, in_dir_map = process_file(
        file_type,
        file,
        adata,
        slot_map,
        in_dir_map,
        is_default=(file == default_files.get(file_type)),
    )

logging.info(f'Write to {output_file}...')
adata.uns = read_anndata(file, uns='uns', verbose=False).uns
add_wildcards(adata, snakemake.wildcards, 'preprocessing')
write_zarr_linked(
    adata=adata,
    in_dir=file_map[0][1],
    out_dir=output_file,
    files_to_keep=['uns/wildcards'],
    slot_map=slot_map,
    in_dir_map=in_dir_map,
)
