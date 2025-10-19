from pathlib import Path
import pandas as pd

from utils.io import read_anndata, link_zarr, write_zarr_linked, ALL_SLOTS


input_anndata = snakemake.input[0]
input_scrublet = snakemake.input.get('scrublet')
input_doubletdetection = snakemake.input.get('doubletdetection')
output_zarr = snakemake.output.zarr
layer = snakemake.params.get('layer', 'X')

if input_anndata.endswith('.h5ad'):
    kwargs = {x: x for x in ALL_SLOTS} | dict(X=layer)
else:
    kwargs = dict(obs='obs')

adata = read_anndata(input_anndata, **kwargs)

if adata.n_obs == 0:
    write_zarr(adata, output_zarr)
    exit(0)

if input_scrublet:
    scrub_scores = pd.concat([pd.read_table(f, index_col=0) for f in input_scrublet])
    scrub_scores.index = scrub_scores.index.astype(str)
    scrub_scores['scrublet_prediction'] = scrub_scores['scrublet_prediction'].astype(str)
    print(scrub_scores)
    adata.obs = adata.obs.merge(scrub_scores, left_index=True, right_index=True, how='left')

if input_doubletdetection:
    doub_scores = pd.concat([pd.read_table(f, index_col=0) for f in input_doubletdetection])
    doub_scores.index = doub_scores.index.astype(str)
    print(doub_scores)
    adata.obs = adata.obs.merge(doub_scores, left_index=True, right_index=True, how='left')

print(adata.obs)

write_zarr_linked(
    adata,
    input_anndata,
    output_zarr,
    files_to_keep=['obs'],
    slot_map={'X': layer},
)