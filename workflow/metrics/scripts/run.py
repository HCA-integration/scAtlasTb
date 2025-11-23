import numpy as np
import pandas as pd
import logging
logger = logging.getLogger('Run metric')
logger.setLevel(logging.INFO)
try:
    from sklearnex import patch_sklearn
    patch_sklearn()
    logging.getLogger("sklearnex").setLevel(logging.WARNING)
except ImportError:
    logger.info('no hardware acceleration for sklearn')

from metrics import metric_map
from metrics.utils import write_metrics
from utils.io import read_anndata


def replace_last(source_string, replace_what, replace_with, cluster_key='leiden'):
    """
    Helper function for parsing clustering columns
    adapted from
    https://stackoverflow.com/questions/3675318/how-to-replace-some-characters-from-the-end-of-a-string/3675423#3675423
    """
    if not source_string.startswith(cluster_key):
        return source_string
    if not source_string.endswith(replace_what):
        return source_string + '_orig'
    head, _sep, tail = source_string.rpartition(replace_what)
    return head + replace_with + tail


input_file = snakemake.input[0]
output_file = snakemake.output[0]
wildcards = snakemake.wildcards
params = snakemake.params
threads = snakemake.threads

dataset = wildcards.dataset
file_id = wildcards.file_id
batch_key = wildcards.batch
label_key = wildcards.label
metric = wildcards.metric

metric_type = params.get('metric_type')
assert metric_type in ['batch_correction', 'bio_conservation'], f'Unknown metric_type: {metric_type}'
allowed_output_types = params.get('output_types')
input_type = params.get('input_type')
comparison = params.get('comparison', False)
cluster_key = params.get('cluster_key', 'leiden')
use_covariate = params.get('use_covariate', False)
use_gene_set = params.get('use_gene_set', False)
covariates = params.get('covariates', [])
gene_sets = params.get('gene_sets', {})

metric_function = metric_map[metric]

uns = read_anndata(input_file, uns='uns', verbose=False).uns
output_type = uns.get('output_type', 'full') # Same as in prepare.py

if output_type not in allowed_output_types:
    logging.info(
        f'Skip metric={metric} for data output type={output_type}\n'
        f'allowed output types={allowed_output_types}'
    )
    write_metrics(
        metric_names=[metric],
        scores=[np.nan],
        output_types=[output_type],
        metric=metric,
        metric_type=metric_type,
        batch=batch_key,
        label=label_key,
        file_id=file_id,
        dataset=dataset,
        filename=output_file,
        **uns.get('wildcards', {}),
        )
    exit(0)

kwargs = {'obs': 'obs', 'uns': 'uns'}
if 'knn' in input_type:
    kwargs |= {'obsp': 'obsp'}
if 'embed' in input_type:
    kwargs |= {'obsm': 'obsm'}
if 'full' in input_type:
    kwargs |= {'X': 'X', 'var': 'var', 'dask': True, 'backed': True}

logger.info(f'Read {input_file} of input_type {input_type}...')
adata = read_anndata(input_file, **kwargs)
if 'feature_name' in adata.var.columns:
    adata.var_names = adata.var['feature_name'].astype(str)
print(adata, flush=True)

adata_raw = None
if comparison:
    logging.info(f'Read raw data for comparison from {input_file} ...')
    adata_raw = read_anndata(
        input_file,
        X='raw/X',
        obs='raw/obs',
        obsm='raw/obsm',
        var='raw/var',
        varm='raw/varm',
        uns='raw/uns',
        dask=True,
        backed=True,
    )
    if 'feature_name' in adata_raw.var.columns:
        adata_raw.var_names = adata_raw.var['feature_name'].astype(str)
    adata_raw.obs[batch_key] = adata_raw.obs[batch_key].astype(str)


# set default covariates
if use_covariate and len(covariates) == 0:
    covariates = [label_key, batch_key]

gene_score_columns = []
if use_gene_set:
    gene_score_columns = [x for x in adata.obs.columns if x.startswith('gene_score:')]

# subset obs columns for metrics
cluster_columns = [col for col in adata.obs.columns if col.startswith(cluster_key) and col.endswith('_1')]
columns = [batch_key, label_key] + covariates + cluster_columns + gene_score_columns
columns = list(set(columns))  # make unique
adata.obs = adata.obs[columns].copy()
adata.obs.rename(columns=lambda x: replace_last(x, '_1', ''), inplace=True)

logger.info(f'Run metric {metric} for {output_type}...')
adata.obs[batch_key] = adata.obs[batch_key].astype(str)

# TODO: deal with bootstrapping
scores = metric_function(
    adata,
    output_type,
    batch_key=batch_key,
    label_key=label_key,
    adata_raw=adata_raw,
    cluster_key=cluster_key,
    covariate=covariates,
    gene_set=gene_sets,
    n_threads=threads,
)

# unpack scores
if isinstance(scores, tuple):
    scores, metric_names = scores
elif isinstance(scores, (int, float, np.number)):
    scores = [scores]
    metric_names = [metric]
elif isinstance(scores, list):
    metric_names = [metric] * len(scores)
else:
    raise ValueError(f'Unknown scores type {type(scores)}')

write_metrics(
    metric_names=metric_names,
    scores=scores,
    output_types=[output_type]*len(scores),
    metric=metric,
    metric_type=metric_type,
    batch=batch_key,
    label=label_key,
    file_id=file_id,
    dataset=dataset,
    filename=output_file,
    **adata.uns.get('wildcards', {}),
)