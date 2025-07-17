import logging
logging.basicConfig(level=logging.INFO)
from pathlib import Path
import traceback
import warnings
warnings.filterwarnings("ignore")

from matplotlib import pyplot as plt
import scanpy as sc
import numpy as np
import pandas as pd
from pandas.api.types import is_numeric_dtype, is_string_dtype, is_categorical_dtype
from pprint import pformat
from tqdm import tqdm
from joblib import Parallel, delayed

from utils.io import read_anndata, get_file_reader
from utils.misc import ensure_dense, remove_outliers, dask_compute


input_file = snakemake.input[0]
output_dir = Path(snakemake.output.plots)
output_dir.mkdir(exist_ok=True)

params = dict(snakemake.params.items())
basis = params['basis']
wildcards_string = '\n'.join([f'{k}: {v}' for k, v in snakemake.wildcards.items()])
logging.info(f'Wildcard string: {wildcards_string}')

logging.info(f'Read file {input_file}...')
kwargs = dict(
    obs='obs',
    obsm='obsm',
    var='var',
)

# check if .X exists
read_func, _ = get_file_reader(input_file)
if 'X' in read_func(input_file, 'r'):
    kwargs |= dict(X='X', dask=True, backed=True)

logging.info(f'Read {input_file}...')
adata = read_anndata(input_file, **kwargs)
assert basis in adata.obsm.keys(), f'"{basis}" not in adata.obsm'
n_cells = adata.n_obs # save original number of cells
obs_columns = adata.obs.columns.tolist()
ensure_dense(adata, basis)

if adata.n_obs == 0:
    logging.info('No cells, skip...')
    exit()

if 'feature_name' in adata.var.columns:
    adata.var_names = adata.var['feature_name'].astype(str)

# set minimum cells per category
min_cells_per_category = params.pop('min_cells_per_category', 1e-4)
if min_cells_per_category < 1:
    min_cells_per_category *= n_cells
print('Remove categories with fewer than', min_cells_per_category, 'cells')

# parse colors
colors = params.get('color', [None])

gene_colors = sorted([col for col in colors if col not in obs_columns])
pattern = '|'.join(gene_colors)
gene_colors = adata.var_names[adata.var_names.astype(str).str.contains(pat=pattern)].tolist()

if 'color' in params:
    logging.info(f'Configured colors:\n{pformat(colors)}')
    colors = colors if isinstance(colors, list) else [colors]
    colors = [color for color in colors if color in obs_columns]
    # filter colors with too few or too many categories
    logging.info(f'Colors after filtering:\n{pformat(colors)}')
    
    for color in colors:
        if color not in obs_columns:
            continue
        column = adata.obs[color]

        if is_categorical_dtype(column) or is_string_dtype(column):
            # parse data types
            column = column.astype(object) \
                .replace(['NaN', 'None', '', 'nan', 'unknown'], float('nan')) \
                .astype('category')

            # remove categories with too few cells
            value_counts = column.value_counts()
            categories_to_remove = value_counts[value_counts <= min_cells_per_category].index
            column = column.cat.remove_categories(categories_to_remove)
            
            # update column
            adata.obs[color] = column
            # adata.obs[color] = column.codes if len(column.categories) > 102 else column

    if len(colors) == 0:
        logging.info('No valid colors, skip...')
        colors = [None]

    del params['color']


logging.info('Remove outliers...')
outlier_factor = params.pop('outlier_factor', 0)
adata = remove_outliers(adata, 'max', factor=outlier_factor, rep=basis)
adata = remove_outliers(adata, 'min', factor=outlier_factor, rep=basis)

# match gene patterns
n_gene_colors = len(gene_colors)
if n_gene_colors > 0:
    logging.info(f'Subset to {n_gene_colors} requested genes...')
    adata = dask_compute(adata[:, adata.var_names.isin(gene_colors)].copy(), layers='X')
    print(adata.var, flush=True)

logging.info('Shuffle cells...')
if n_cells > 1e6:
    adata = adata[adata.obs.sample(frac=0.7).index]

if adata.is_view:
    logging.info('Convert view to copy...')
    adata = adata.copy()

# set minimum point size
default_size = 200_000 / adata.n_obs
size = params.get('size', default_size)
if size is None:
    size = default_size
params['size'] = np.min([np.max([size, 0.4, default_size]), 200])


def plot_color(color, title, output_dir, verbose=True):
    palette = None
    if title is None:
        title = str(color)
    colors = color if isinstance(color, list) else [color]

    fig_params = dict(
        frameon=False,
        vector_friendly=True,
        fontsize=9,
        figsize=(6,6),
        dpi=200,
        dpi_save=200,
        format='png',
    )
    
    if len(colors) > 4:
        fig_params = dict(frameon=False)
    sc.set_figure_params(**fig_params)

    # select palette according to obs column
    for color in colors:
        if color in adata.obs.columns:
            color_vec = adata.obs[color]
            if is_categorical_dtype(color_vec):
                if color_vec.nunique() > 102:
                    palette = 'turbo'
                elif color_vec.nunique() > 20:
                    palette = sc.pl.palettes.godsnot_102
            elif is_numeric_dtype(color_vec):
                if color_vec.min() < 0:
                    palette = 'coolwarm'
                else:
                    palette = 'plasma'
    
    try:
        fig = sc.pl.embedding(
            adata,
            color=colors,
            show=False,
            return_fig=True,
            palette=palette,
            **params
        )

        legend = fig.get_axes()[0].get_legend()
        if palette == 'turbo':
            legend.remove()
        elif legend:
            legend_bbox = legend.get_window_extent()
            fig_width, fig_height = fig.get_size_inches()
            fig_width = fig_width + (legend_bbox.width / fig.dpi)
            fig.set_size_inches((fig_width, fig_height))
        
        # set title & adjust figure layout
        suptitle = fig.suptitle(f'{wildcards_string}\nn={n_cells}')
        fig.canvas.draw()
        bbox = suptitle.get_window_extent(renderer=fig.canvas.get_renderer())
        bbox_inches = bbox.transformed(fig.transFigure.inverted())

        fontsize_inches = suptitle.get_size() / fig.dpi
        top_margin = bbox_inches.y1 + fontsize_inches * 1.5
        top_margin = min(top_margin, 1.0)  # clamp to avoid errors

        fig.tight_layout(rect=[0, 0, 1, top_margin])

        if verbose:
            logging.info(f'Plotting color "{title}" successful.')
    
    except Exception as e:
        logging.error(f'Failed to plot {title}: {e}')
        traceback.print_exc()
        plt.plot([])
    
    out_path = output_dir / f'{title}.png'
    
    try:        
        plt.savefig(out_path, bbox_inches='tight')
    except Exception as e:
        logging.error(f'Failed to save plot "{title}" to {out_path}: {e}')
        traceback.print_exc()

logging.info('Parameters:\n' + pformat(params))
# Run plotting in parallel
list(tqdm(
    Parallel(return_as='generator')(delayed(plot_color)(
        color,
        title=color,
        output_dir=output_dir,
    ) for color in set(colors)),
    desc="Plotting colors",
    total=len(colors),
    miniters=1,
))

chunk_size = 24
gene_colors = {
    f'genes_group={idx}': gene_colors[i:i + chunk_size]
    for idx, i in enumerate(range(0, len(gene_colors), chunk_size))
}
if gene_colors:
    Path(output_dir / 'genes').mkdir(exist_ok=True)
    list(tqdm(
        Parallel(return_as='generator')(delayed(plot_color)(
            color,
            title=title,
            output_dir=output_dir / 'genes',
            verbose=False
        ) for title, color in gene_colors.items()),
        desc="Plotting genes",
        total=len(gene_colors),
        miniters=1,
    ))