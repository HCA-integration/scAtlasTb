from pathlib import Path
import numpy as np
import pandas as pd
import doubletdetection
import anndata as ad
from pprint import pformat
from tqdm import tqdm
import logging
logging.basicConfig(level=logging.INFO)

from utils.io import read_anndata
from utils.misc import dask_compute

input_zarr = snakemake.input.zarr
output_tsv = snakemake.output.tsv
batches_txt = snakemake.input.batch
layer = snakemake.params.get('layer', 'X')
batch_key = snakemake.params.get('batch_key')
threads = snakemake.threads


def run_method(adata, batch_key, batch, min_cells=100, threads=threads):
    if batch_key in adata.obs.columns:
        adata = adata[adata.obs[batch_key].astype(str) == batch, :].copy()

    if adata.n_obs < min_cells:
        columns = ['doubletdetection_score', 'doubletdetection_prediction']
        return pd.DataFrame(index=adata.obs.index, columns=columns, dtype=float).fillna(0)

    adata = dask_compute(adata, verbose=False)

    clf = doubletdetection.BoostClassifier(
        n_iters=10,
        n_top_var_genes=4000,
        n_components=np.min([adata.n_obs, adata.n_vars, 30]),
        clustering_algorithm="leiden",
        clustering_kwargs=dict(flavor='igraph', n_iterations=2),
        n_jobs=threads,
    )
    labels = clf.fit(adata.X).predict(p_thresh=1e-16, voter_thresh=0.5)
    scores = clf.doublet_score()

    adata.obs['doubletdetection_score'] = scores
    adata.obs['doubletdetection_prediction'] = labels
    return adata.obs[['doubletdetection_score', 'doubletdetection_prediction']]


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
    adata = dask_compute(adata, layers='X')
    logging.info(f'Subset to {adata.n_obs} cells with specified batches.')

logging.info('Run doubletdetection...')
results = (
    run_method(adata, batch_key, batch)
    for batch in tqdm(batches, desc='Run doubletdetection', miniters=1)
)
df = pd.concat(results)

df.to_csv(output_tsv, sep='\t')
