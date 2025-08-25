import pandas as pd
from scipy import sparse
import numpy as np
import scib
from anndata import AnnData
import yaml
from tqdm import tqdm
# from concurrent.futures import ProcessPoolExecutor, as_completed
# from joblib import Parallel, delayed
# try:
#     from sklearnex import patch_sklearn
#     patch_sklearn()
# except ImportError:
#     print('no hardware acceleration for sklearn', flush=True)

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
if not isinstance(adata.X, np.ndarray):
    # make sure the PCA embedding is an array
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
    adata.obs[sample_key] = adata.obs[sample_keys].astype(str).agg('-'.join, axis=1)

# read setup file
with open(setup_file, 'r') as f:
    setup = yaml.safe_load(f)
# n_permute = setup['n_permute']
# n_permute = min(snakemake.params.get('n_permute', 0), n_permute)
n_permute = snakemake.params.get('n_permute', 0)
print(f'n_permute: {n_permute}', flush=True)

# Permute covariates
perm_covariates = []

if covariate == sample_key:
    print('Sample key is the same as covariate, skipping permutation...', flush=True)
    n_permute = 0

series_dict = {}
for i in range(n_permute):
    # aggregate and permutate covariate
    cov_per_sample = adata.obs.groupby(sample_key, observed=True).agg({covariate: 'first'})
    cov_map = dict(zip(cov_per_sample.index, cov_per_sample[covariate].sample(frac=1)))
    
    # apply permutation
    covariate_perm = f'{covariate}-{i}'
    series_dict[covariate_perm] = adata.obs[sample_key].map(cov_map)
    perm_covariates.append(covariate_perm)

# concat new columns to adata to avoid fragmentation
adata.obs = pd.concat(
    [
        adata.obs[[sample_key, covariate]],
        pd.DataFrame(series_dict)
    ],
    axis=1
).copy()


pcr_scores = []
for covariate in tqdm([covariate]+perm_covariates):
    is_permuted = covariate in perm_covariates
    score = scib.me.pcr(
        adata,
        covariate=covariate,
        recompute_pca=False,
        verbose=False,
        linreg_method='numpy',
    )
    pcr_scores.append((covariate, is_permuted, score))

# pcr_scores = tqdm(
#     Parallel(
#         n_jobs=n_threads,
#         require='sharedmem',
#         return_as='generator'
#     )(
#         delayed(compute_pcr)(adata, covariate, perm_covariates)
#         for covariate in [covariate]+perm_covariates
#     ),
#     total=1+len(perm_covariates),
# )
# pcr_scores = list(pcr_scores)

# Set permuted score when covariate is the same as the group variable
if covariate == sample_key:
    # permutations wouln't change values in this case, impute same value as covariate score
    perm_score = (f'{covariate}-0', True, pcr_scores[0][2])
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

print(df.to_string(), flush=True)
df.to_csv(output_file, sep='\t', index=False)
