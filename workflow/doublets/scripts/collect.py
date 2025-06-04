from pathlib import Path
import pandas as pd

from utils.io import read_anndata, link_zarr


input_anndata = snakemake.input[0]
input_scrublet = snakemake.input.scrublet
input_doubletdetection = snakemake.input.doubletdetection
output_zarr = snakemake.output.zarr


# read AnnData
kwargs = dict() if input_zarr.endswith('.h5ad') else dict(obs='obs')
adata = read_anndata(input_zarr, **kwargs)

if adata.n_obs == 0:
    adata.write_zarr(output_zarr)
    exit(0)

scrub_scores = pd.concat([pd.read_table(f, index_col=0) for f in input_scrublet])
scrub_scores.index = scrub_scores.index.astype(str)
scrub_scores['scrublet_prediction'] = scrub_scores['scrublet_prediction'].astype(str)

doub_scores = pd.concat([pd.read_table(f, index_col=0) for f in input_doubletdetection])
doub_scores.index = doub_scores.index.astype(str)

print(scrub_scores)
print(doub_scores)

adata.obs = adata.obs.merge(scrub_scores, left_index=True, right_index=True, how='left')
adata.obs = adata.obs.merge(doub_scores, left_index=True, right_index=True, how='left')
print(adata.obs)

adata.write_zarr(output_zarr)

if input_anndata.endswith('.zarr'):
    input_files = [f.name for f in Path(input_anndata).iterdir()]
    files_to_keep = [f for f in input_files if f not in ['obs']]
    link_zarr(
        in_dir=input_anndata,
        out_dir=output_zarr,
        file_names=files_to_keep,
        overwrite=True,
    )