import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import sctk
import logging
logging.basicConfig(level=logging.INFO)

from utils.io import read_anndata, write_zarr_linked

input_file = snakemake.input[0]
output_file = snakemake.output[0]
gauss_threshold = snakemake.params['gauss_threshold']
layer = snakemake.params['layer']

files_to_keep = ['obs', 'uns']

adata = read_anndata(
    input_file,
    X=layer,
    obs='obs',
    var='var'
)

if adata.n_obs == 0:
    logging.info(f'Write empty zarr file to {output_file}...')
    adata.obs = pd.DataFrame(
        columns=adata.obs.columns.tolist() \
            +["n_counts", "n_genes", "percent_mito", "percent_ribo", "percent_hb"]
    )
    write_zarr_linked(adata, input_file, output_file, files_to_keep=files_to_keep)
    exit(0)

print('Calculate QC stats...')
if 'feature_name' in adata.var.columns:
    var_names = adata.var['feature_name'].astype(str)
else:
    var_names = adata.var_names

adata.var["mito"] = var_names.str.startswith("MT-")
# sc.pp.calculate_qc_metrics(adata, qc_vars=["mito"], inplace=True)

logging.info('Calculate QC metrics...')
sctk.calculate_qc(adata)
metrics = sctk.default_metric_params_df.loc[["n_counts", "n_genes", "percent_mito", "percent_ribo", "percent_hb"], :]

logging.info('Calculate cell-wise QC...')
print(metrics)
sctk.cellwise_qc(adata, metrics=metrics, threshold=gauss_threshold)
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