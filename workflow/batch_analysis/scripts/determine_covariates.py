import logging
logging.basicConfig(level=logging.INFO)
from pprint import pformat
from pandas.api.types import is_numeric_dtype
from pathlib import Path
from utils.io import read_anndata

input_file = snakemake.input[0]
output_dir = snakemake.output[0]
sample_key = 'group'
Path(output_dir).mkdir(parents=True)

covariates = snakemake.params.get('covariates', [])
perm_covariates = snakemake.params.get('permute_covariates')
if perm_covariates is None:
    perm_covariates = covariates
n_perms = snakemake.params.get('n_permute')
na_strings = snakemake.params.get('na_strings', [])

logging.info(f'Read {input_file}...')
obs = read_anndata(input_file, obs='obs', verbose=False).obs
if sample_key not in obs.columns:
    raise KeyError(
        f'Missing required "{sample_key}" column in obs. '
        'Run sample_representation/prepare first to create it.'
    )

def covariate_valid(obs, covariate, na_strings):
    if covariate not in obs.columns:
        return False

    obs[covariate] = obs.loc[
        ~obs[covariate].isin(na_strings) & obs[covariate].notna(),
        covariate
    ]
    
    if is_numeric_dtype(obs[covariate]):
        return obs[covariate][obs[covariate].notna()].nunique() > 1
    
    nonunique_map = (
        obs.groupby(sample_key, observed=False)[covariate]
        .unique()
        .apply(sorted)
        .loc[lambda x: x.str.len() > 1]
    )
    if nonunique_map.shape[0] > 0:
        logging.warning(f'Skipping covariate "{covariate}" because only one value per sample "{sample_key}" allowed')
        logging.debug(f'\n{nonunique_map.head()}')
        return False

    return obs[covariate][obs[covariate].notna()].nunique() >= 2

covariates = [
    covariate for covariate in covariates
    if covariate_valid(obs, covariate, na_strings)
]
perm_covariates = [
    covariate for covariate in perm_covariates
    if covariate_valid(obs, covariate, na_strings) and covariate != sample_key
]

logging.info(f'covariates:\n {pformat(covariates)}')
logging.info(f'perm_covariates:\n {pformat(perm_covariates)}')
logging.info(f'n_permutations:\n {pformat(n_perms)}')

for covariate in covariates + perm_covariates:
    covariate_file = f'{output_dir}/{covariate}.yaml'
    with open(covariate_file, 'w') as f:
        n_permute = n_perms if covariate in perm_covariates else 0
        is_numeric = bool(is_numeric_dtype(obs[covariate]))
        na_strings_yaml = str(na_strings).replace("'", '"')
        f.write(f'n_permute: {n_permute}\nis_numeric: {str(is_numeric).lower()}\nna_strings: {na_strings_yaml}\n')
