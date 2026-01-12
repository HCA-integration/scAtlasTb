import pandas as pd
from scipy import sparse
import numpy as np
import scib
from anndata import AnnData
import yaml
from tqdm import tqdm
from joblib import Parallel, delayed

import logging
import gc
logging.basicConfig(level=logging.INFO)

from utils.io import read_anndata


input_file = snakemake.input[0]
setup_file = snakemake.input.setup
output_file = snakemake.output.tsv
covariate = snakemake.wildcards.covariate
sample_key = snakemake.params.get('sample_key')
n_threads = np.max([snakemake.threads, 1])

print('Read anndata file...', flush=True)
adata = read_anndata(
    input_file,
    X='obsm/X_pca',
    obs='obs',
    uns='uns',
)
assert 'pca' in adata.uns, f'.uns["pca"] is missing, please make sure PCA is computed and PC loadings are saved in adata.uns as provided by scanpy'
adata.uns = {'pca': adata.uns['pca']}

# make sure the PCA embedding is an array
if not isinstance(adata.X, np.ndarray):
    adata.X = adata.X.toarray()
adata.obsm['X_pca'] = adata.X
del adata.X

adata = adata[adata.obs[covariate].notna()].copy()
n_covariate = adata.obs[covariate].nunique()

# set default sample key
if sample_key is None or sample_key == 'None':
    print('Using index as sample key...', flush=True)
    sample_key = 'index'
    adata.obs[sample_key] = adata.obs.index
sample_keys = [x.strip() for x in sample_key.split(',')]
if len(sample_keys) > 1:
    sample_key = '--'.join(sample_keys)
    adata.obs[sample_key] = adata.obs[sample_keys].astype(str).agg('-'.join, axis=1).astype('category')

# minimise adata.obs
select_columns = list({sample_key, covariate})
adata.obs = adata.obs[select_columns].copy()

# read setup file
with open(setup_file, 'r') as f:
    setup = yaml.safe_load(f)
n_permute = setup['n_permute']


def is_categorical_series(s, max_unique=10):
    return (
        pd.api.types.is_categorical_dtype(s) or
        pd.api.types.is_object_dtype(s) or
        (pd.api.types.is_numeric_dtype(s) and s.nunique() <= max_unique)
    )

if not is_categorical_series(adata.obs[covariate]):
    adata.obs[covariate] = adata.obs[covariate].astype(np.float32)
else:
    nonunique_map = (
        adata.obs.groupby(sample_key, observed=False)[covariate]
        .unique()
        .apply(sorted)
        .loc[lambda x: x.str.len() > 1]
    )
    assert nonunique_map.shape[0] == 0, \
        f'Each sample (defined by {sample_key}) must have exactly one value for covariate {covariate}, '  \
        f'but found values:\n{nonunique_map}'

    value_counts = adata.obs[[sample_key, covariate]].drop_duplicates().value_counts(covariate)
    logging.info(value_counts[value_counts > 1])

    if value_counts.max() == 1:
        logging.info('Sample key is the same as covariate, skipping permutation...')
        n_permute = 0


logging.info(f'Calculating PCR scores for {n_permute} permutations (using {n_threads} threads)')

X_pca = adata.obsm['X_pca'].astype(np.float32)
pca_var = adata.uns['pca']['variance']
obs_cov = adata.obs[[sample_key, covariate]].copy()

sample_codes, uniques = pd.factorize(obs_cov[sample_key], sort=False)
cov_per_sample = obs_cov \
    .groupby(sample_key, observed=True)[covariate] \
    .first() \
    .reindex(uniques) \
    .to_numpy()

def run_pcr(i, seed, **kwargs):
    if i > 0:
        rng = np.random.default_rng([i, seed])
        covariate_array = rng.permutation(cov_per_sample)[sample_codes]
    else:
        covariate_array = obs_cov[covariate].values

    return f"{covariate}-{i}", i > 0, scib.me.pc_regression(
        X_pca,
        pca_var=pca_var,
        covariate=covariate_array,
        **kwargs,
    )
    

def chunked_parallel_jobs(jobs, chunk_size, n_threads):
    results = []
    for i in tqdm(range(0, len(jobs), chunk_size), desc='PCR chunks'):
        chunk = jobs[i:i+chunk_size]
        chunk_results = Parallel(
            n_jobs=n_threads,
            require='sharedmem',
        )(chunk)
        results.extend(chunk_results)
    return results


all_jobs = [
    delayed(run_pcr)(
        i=i,
        seed=42,
        verbose=False,
        linreg_method='numpy',
    ) for i in range(n_permute + 1)
]

chunk_size = max(1, n_threads)
pcr_scores = chunked_parallel_jobs(all_jobs, chunk_size, n_threads)

# Set permuted score when covariate is the same as the group variable
if n_permute == 0:
    # permutations wouldn't change values in this case, impute same value as covariate score
    perm_score = (f'{covariate}-1', True, pcr_scores[0][2])
    pcr_scores.append(perm_score)

df = pd.DataFrame.from_records(
    pcr_scores,
    columns=['covariate', 'permuted', 'pcr'],
)

# calculate summary stats
df['covariate'] = df['covariate'].str.split('-', expand=True)[0].astype('category')
df['n_covariates'] = n_covariate
df['perm_mean'] = df.loc[df['permuted'], 'pcr'].mean()
df['perm_std'] = df.loc[df['permuted'], 'pcr'].std()
df['z_score'] = (df['pcr'] - df['perm_mean']) / df['perm_std']
df['signif'] = df['z_score'] > 1.5
if n_permute == 0:
    df['p-val'] = np.nan
else:
    # df['p-val'] = df.loc[df['permuted'], 'signif'].sum() / n_permute
    stat = np.abs(df['pcr'] - df['perm_mean'])
    null = stat.loc[df['permuted']].values
    obs = stat.loc[~df['permuted']].values[0]
    df['p-val'] = (np.sum(null >= obs) + 1) / (n_permute + 1)
    df['signif'] = df['p-val'] <= 0.05

print(df, flush=True)
df.to_csv(output_file, sep='\t', index=False)
