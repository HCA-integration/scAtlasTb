import numpy as np
import scanpy as sc
import anndata
from pathlib import Path
from matplotlib import pyplot as plt
from pprint import pformat
import logging
logging.basicConfig(level=logging.INFO)

from utils.io import read_anndata
from utils.misc import dask_compute, ensure_dense


sc.set_figure_params(frameon=False)
plt.rcParams['figure.figsize'] = 30, 25
input_file = snakemake.input[0]
output_png = snakemake.output.dotplot
markers = snakemake.params.markers

Path(output_png) .mkdir(parents=True, exist_ok=True)

wildcards = snakemake.wildcards
group_key = wildcards.group
file_id = wildcards.file_id
title = f'Marker genes for file_id={file_id} group={group_key}'

layer = snakemake.params.layer
args = snakemake.params.get('args', {})
n_groups_per_split = args.pop('n_groups_per_split', 10)

logging.info(f'Reading {input_file}...')
adata = read_anndata(
    input_file,
    X=layer,
    obs='obs',
    var='var',
    uns='uns',
    dask=True,
    backed=True,
)

if 'feature_name' in adata.var.columns:
    adata.var_names = adata.var['feature_name']

genes_to_plot = []
gene_sets = {'markers': {}}
if isinstance(markers, list):
    gene_sets['markers']['generic'] = markers
    genes_to_plot.append(markers)
elif isinstance(markers, dict):
    for k, v in markers.items():
        if isinstance(v, list):
            gene_sets['markers'][k] = v  # collect in generic markers set
            genes_to_plot.extend(v)
        elif isinstance(v, dict):
            gene_sets[k] = v
            all_genes = [element for sublist in v.values() for element in sublist]
            genes_to_plot.extend(all_genes)
        else:
            raise ValueError(f'Invalid markers format\n{k}: {v}')
else:
    raise ValueError(f'Invalid markers format\n{markers}')

if not gene_sets['markers']:
    del gene_sets['markers']

# subset to genes for plotting
genes_to_plot = list(set(genes_to_plot))
logging.info(f'{len(genes_to_plot)} genes requested...')
gene_mask = adata.var_names.isin(genes_to_plot)
logging.info(f'{gene_mask.sum()} genes found overlapping...')
adata = adata[:, gene_mask].copy()

logging.info('Load genes to memory...')
adata = dask_compute(adata)
adata = ensure_dense(adata)
logging.info(adata.__str__())

if min(adata.shape) == 0:
    logging.info('No data, exit...')
    print(adata, flush=True)
    exit()

for set_name, gene_set in gene_sets.items():
    logging.info(f'Processing gene set "{set_name}"...')
    # match marker genes and var_names
    print(
        pformat({k: len(v) for k, v in gene_set.items()}),
        flush=True
    )
    gene_set = {
        k: adata.var_names[adata.var_names.isin(v)].to_list()
        for k, v in gene_set.items()
    }
    
    # remove groups without matching genes
    gene_set = {k: v for k, v in gene_set.items() if len(v) > 0}
    print(
        f'matching genes:\n',
        pformat({k: len(v) for k, v in gene_set.items()}),
        flush=True
    )
    
    if not gene_set:
        logging.info(f'No genes in marker gene set {set_name}, skipping...\n{gene_set}')
        continue

    logging.info('Dotplot...')
    dp = sc.pl.dotplot(
        adata,
        gene_set,
        groupby=group_key,
        use_raw=False,
        standard_scale='var',
        # title=title,
        show=False,
        swap_axes=False,
        return_fig=True,
    )
    dp.add_totals().style(dot_edge_color='black', dot_edge_lw=0.5)
    dp.savefig(f'{output_png}/{set_name}.png', dpi=100, bbox_inches='tight')
