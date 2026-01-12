from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad
from pprint import pformat
from tqdm import tqdm
import traceback
import logging
logging.basicConfig(level=logging.INFO)

from utils.io import read_anndata
from utils.processing import sc as rsc
import scanpy as sc
from utils.processing import USE_GPU
from utils.misc import dask_compute

input_zarr = snakemake.input.zarr
output_tsv = snakemake.output.tsv
batches_txt = snakemake.input.batch
layer = snakemake.params.get('layer', 'X')
batch_key = snakemake.params.get('batch_key')

kwargs = dict(
    batch_key=None,
    sim_doublet_ratio=2.0,
    expected_doublet_rate=0.05,
    stdev_doublet_rate=0.02,
    synthetic_doublet_umi_subsampling=1.0,
    knn_dist_metric='euclidean',
    normalize_variance=True,
    log_transform=False,
    mean_center=True,
    use_approx_neighbors=True,
    get_doublet_neighbor_parents=False,
    n_neighbors=None,
    threshold=None,
    verbose=True,
    copy=False,
    random_state=0,
)


def run_method(
    adata,
    batch_key,
    batch,
    min_cells=100,
    use_gpu=USE_GPU,
    **kwargs
):
    if batch_key in adata.obs.columns:
        adata = adata[adata.obs[batch_key].astype(str) == batch, :].copy()

    if adata.n_obs < min_cells:
        columns = ['scrublet_score', 'scrublet_prediction']
        return pd.DataFrame(index=adata.obs.index, columns=columns, dtype=float).fillna(0)

    adata = dask_compute(adata, verbose=False)
    
    kwargs |= dict(
        n_prin_comps=np.min([30, adata.n_obs-1, adata.n_vars-1])
    )

    if use_gpu:
        rsc.get.anndata_to_GPU(adata)

    try: 
        rsc.pp.scrublet(adata, **kwargs)
    except Exception as e:    
        if not use_gpu:
            raise e

        traceback.print_exc()
        print(adata, flush=True)

        logging.info('Retry on CPU...')
        rsc.get.anndata_to_CPU(adata)
        sc.pp.scrublet(adata, **kwargs)

    df = adata.obs[['doublet_score', 'predicted_doublet']].rename(
        columns={
            'doublet_score': 'scrublet_score',
            'predicted_doublet': 'scrublet_prediction',
        }
    )

    return df


logging.info(f'Read {input_zarr}...')
adata = read_anndata(
    input_zarr,
    X=layer,
    obs='obs',
    backed=True,
    dask=True,
)

with open(batches_txt, 'r') as f:
    batches = [line.strip().split('\t')[0] for line in f if line.strip()]

if batch_key is not None:
    adata = adata[adata.obs[batch_key].isin(batches), :].copy()
    logging.info(f'Subset to {adata.n_obs} cells with specified batches.')

results = (
    run_method(adata, batch_key, batch)
    for batch in tqdm(batches, desc='Run scrublet', miniters=1)
)
df = pd.concat(results)

df.to_csv(output_tsv, sep='\t')
