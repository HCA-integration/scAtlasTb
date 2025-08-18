from pathlib import Path
import logging
logging.basicConfig(level=logging.INFO)
import numpy as np
import anndata as ad
from tqdm import tqdm
from scipy import sparse

from metrics.bootstrap import get_bootstrap_adata
from utils.io import read_anndata, write_zarr_linked
from utils.processing import sc
from utils.misc import dask_compute
from utils.accessors import parse_gene_names


def bins_by_quantiles(matrix: [np.array, np.matrix], n_quantiles: int):
    assert n_quantiles > 1
    
    def bin_values(x):
        bins = np.linspace(0, 1, n_quantiles + 1)
        quantiles = np.quantile(x, bins)
        binned_idx = np.digitize(x, quantiles) - 1
        return bins[binned_idx]
    
    return np.apply_along_axis(bin_values, axis=0, arr=matrix)


input_file = snakemake.input[0]
output_file = snakemake.output[0]

params = snakemake.params
unintegrated_layer = params.get('unintegrated_layer', 'X')
raw_counts_layer = params.get('raw_counts_layer', unintegrated_layer)
gene_sets = params['gene_sets']
n_random_permutations = params.get('n_permutations', 100)
n_quantiles = params.get('n_quantiles', 2)
var_key = 'metrics_features'

files_to_keep = ['obs', 'obsm']   # TODO: only overwrite specific obsm slots


logger.info(f'Read {input_file} ...')
adata = read_anndata(
    input_file,
    X=unintegrated_layer,
    obs='obs',
    obsm='obsm',
    var='var',
    dask=True,
    backed=True,
)

# adata.layers['raw_counts'] = read_anndata(
#     input_file,
#     X=raw_counts_layer,
#     dask=True,
#     backed=True,
# ).X
 
if 'feature_name' in adata.var.columns:
    adata.var['feature_id'] = adata.var_names
    adata.var_names = adata.var['feature_name'].astype(str).values

all_obs_names = adata.obs_names.copy()
all_var_names = adata.var_names.copy()

# filter all gene sets to genes in adata
for set_name, gene_list in gene_sets.items():
    gene_sets[set_name] = parse_gene_names(adata, gene_list)

# get all genes from gene sets
genes = list(set().union(*gene_sets.values()))

if len(genes) == 0:
    logging.warning('No genes of interest found in dataset')
    write_zarr_linked(
        adata,
        in_dir=input_file,
        out_dir=output_file,
    )
    exit(0)

# deal with duplicate gene names
if adata.var_names.duplicated().sum() > 0:
    adata = adata[:, ~adata.var_names.duplicated()]

# subset adata if dataset too large TODO: move to prepare script?
n_subset = int(4e6)
if adata.n_obs > n_subset:
    adata = get_bootstrap_adata(adata, size=n_subset)

# subset to HVGs (used for control genes) + genes of interest
adata.var.loc[adata.var_names.isin(genes), var_key] = True
if not adata.var[var_key].all():
    adata = adata[:, adata.var[var_key]].copy()

adata = dask_compute(adata, layers='X')

for set_name, gene_list in tqdm(gene_sets.items(), desc='Compute Gene scores', miniters=1):
    n_genes = len(gene_list)
    if n_genes == 0:
        logging.info(f'Gene set with {n_genes} genes is too small, skip')
        continue
    
    logging.info(f'Gene score for gene set with {len(gene_list)} genes...')
    sc.tl.score_genes(
        adata,
        gene_list=gene_list,
        score_name=f'gene_score:{set_name}'
    )

    logging.info('Add random gene scores...')
    random_key = f'random_gene_scores:{n_genes}'
    if random_key in adata.obsm.keys():
        logging.info(f'Random gene scores for {n_genes} genes already computed, skip')
        continue
    
    scores = []
    for _ in range(n_random_permutations):
        sc.tl.score_genes(
            adata,
            gene_list=np.random.choice(
                adata.var_names,
                size=n_genes,
                replace=False
            )
        )
        scores.append(adata.obs['score'].values)
    del adata.obs['score']
    adata.obsm[random_key] = sparse.csr_matrix(np.vstack(scores).T)

# logging.info('Bin expression...')
adata = adata[:, genes].copy()
# adata = dask_compute(adata[:, genes].copy(), layers='raw_counts')
# adata.X = bins_by_quantiles(
#     adata.layers['raw_counts'].toarray(),
#     n_quantiles=n_quantiles,
# )

logging.info(f'Write to {output_file}...')
logging.info(adata.__str__())
obs_mask = all_obs_names.isin(adata.obs_names)
var_mask = all_var_names.isin(adata.var_names)
write_zarr_linked(
    adata,
    in_dir=input_file,
    out_dir=output_file,
    subset_mask=(obs_mask, var_mask),
    files_to_keep=files_to_keep,
)