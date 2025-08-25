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


input_file = snakemake.input.anndata
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
adata.uns = {k: v for k, v in adata.uns.items() if k == 'pca'}
assert 'pca' in adata.uns, f'.uns["pca"] is missing, please make sure PCA is computed and PC loadings are saved in adata.uns as provided by scanpy'

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
# n_permute = setup['n_permute']
# n_permute = min(snakemake.params.get('n_permute', 0), n_permute)
n_permute = snakemake.params.get('n_permute', 0)

if adata.obs[[sample_key, covariate]].value_counts().max() == 1:
    logging.info('Sample key is the same as covariate, skipping permutation...', flush=True)
    n_permute = 0

logging.info(f'Calculating PCR scores for {n_permute} permutations (using {n_threads} threads)')


def run_pcr(adata, covariate, i, sample_key, sample_ids, cov_values, **kwargs):
    is_permuted = i > 0
    if is_permuted:
        # Permute covariate per sample
        permuted_values = np.random.permutation(cov_values)
        cov_map = dict(zip(sample_ids, permuted_values))
        covariate_array = adata.obs[sample_key].map(cov_map)
    else:
        # Use original covariate
        covariate_array = adata.obs[covariate]

    result = scib.me.pc_regression(
        adata.obsm['X_pca'],
        pca_var=adata.uns['pca']['variance'],
        covariate=covariate_array,
        **kwargs
    )
    
    return f'{covariate}-{i}', is_permuted, result


def chunked_parallel_jobs(jobs, chunk_size, n_threads):
    results = []
    for i in tqdm(range(0, len(jobs), chunk_size), desc='PCR chunks'):
        chunk = jobs[i:i+chunk_size]
        chunk_results = Parallel(
            n_jobs=n_threads,
            require='sharedmem',
        )(chunk)
        results.extend(chunk_results)
        # Explicit cleanup
        del chunk_results
        del chunk
        gc.collect()
    return results


# aggregate covariate per sample for permutation
cov_per_sample = adata.obs.groupby(sample_key, observed=True).agg({covariate: 'first'})

all_jobs = [
    delayed(run_pcr)(
        adata=adata,
        covariate=covariate,
        i=i,
        sample_key=sample_key,
        sample_ids=cov_per_sample.index.values,
        cov_values=cov_per_sample[covariate].values,
        verbose=False,
        linreg_method='numpy',
    ) for i in range(n_permute + 1)
]

chunk_size = max(1, n_threads)
pcr_scores = chunked_parallel_jobs(all_jobs, chunk_size, n_threads)

# Set permuted score when covariate is the same as the group variable
if n_permute == 0 :
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
df['p-val'] = df.loc[df['permuted'], 'signif'].sum() / n_permute

print(df, flush=True)
df.to_csv(output_file, sep='\t', index=False)
