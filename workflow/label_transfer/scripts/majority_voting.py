import logging
logging.basicConfig(level=logging.INFO)
from pathlib import Path
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from pprint import pformat

from utils.io import read_anndata, write_zarr_linked


input_file = snakemake.input.zarr
output_file = snakemake.output.zarr
reference_key = snakemake.params.get('reference_key')
query_key = snakemake.params.get('query_key')
crosstab_kwargs = snakemake.params.get('crosstab_kwargs', {})

logging.info('Read adata...')
adata = read_anndata(input_file, obs='obs')

logging.info('Compute majority reference for...')
logging.info(f'reference_key: {reference_key}')
logging.info(f'query_key: {query_key}')

adata.obs[reference_key] = adata.obs[reference_key].astype(str).replace('nan', float('nan'))

map_majority = pd.crosstab(
    adata.obs[reference_key],
    adata.obs[query_key],
    **crosstab_kwargs
).idxmax(axis=0)

adata.obs['majority_reference'] = pd.Categorical(
    adata.obs[query_key].map(map_majority),
    categories=map_majority.dropna().unique(),
)

logging.info(f'Write to {output_file}...')
write_zarr_linked(
    adata=adata,
    in_dir=input_file,
    out_dir=output_file,
    files_to_keep=['obs']
)