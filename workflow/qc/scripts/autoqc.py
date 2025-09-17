import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import sctk
import matplotlib.pyplot as plt
import logging
logging.basicConfig(level=logging.INFO)
import dask
dask.config.set(scheduler="single-threaded")

from utils.io import read_anndata, write_zarr_linked

input_file = snakemake.input[0]
output_file = snakemake.output[0]
layer = snakemake.params['layer']
metrics_params_file = snakemake.input.get('metrics_params')
gaussian_kwargs = snakemake.params.get('gaussian_kwargs', {})
QC_FLAGS = [
    'n_counts',
    'n_genes',
    'percent_mito',
    'n_counts_mito',
    'percent_ribo',
    'n_counts_ribo',
    'percent_hb',
    'n_counts_hb',
]

files_to_keep = ['obs', 'uns']

adata = read_anndata(
    input_file,
    X=layer,
    obs='obs',
    var='var',
    dask=True,
    backed=True,
)

if adata.n_obs == 0:
    logging.info(f'Write empty zarr file to {output_file}...')
    columns = adata.obs.columns.tolist() + QC_FLAGS
    adata.obs = pd.DataFrame(columns=list(dict.fromkeys(columns)))
    write_zarr_linked(adata, input_file, output_file, files_to_keep=files_to_keep)
    exit(0)

print('Calculate QC stats...')
if 'feature_name' in adata.var.columns:
    adata.var_names = adata.var['feature_name'].astype(str)

logging.info('Calculate QC metrics...')
sctk.calculate_qc(adata)

logging.info('Determine parameters for scAutoQC...')
metrics_params = sctk.default_metric_params_df
if metrics_params_file:
    user_params = pd.read_table(metrics_params_file, index_col=0)
    # update default parameters with user-provided parameters
    metrics_params.update(user_params)

logging.info(f'\n{metrics_params}')

logging.info('Calculate cell-wise QC...')
sctk.cellwise_qc(
    adata,
    metrics=metrics_params,
    **gaussian_kwargs,
)

adata.uns['scautoqc_ranges'] = adata.uns['scautoqc_ranges'].astype('float32')
logging.info(f"\n{adata.uns['scautoqc_ranges']}")

# sctk.generate_qc_clusters(adata, metrics=["log1p_n_counts", "log1p_n_genes", "percent_mito"])
# adata.obs['qc_cell'] = np.where(adata.obs['consensus_passed_qc'], 'pass', 'fail')

logging.info(f'Write zarr file to {output_file}...')
write_zarr_linked(
    adata,
    input_file,
    output_file,
    files_to_keep=files_to_keep,
)
