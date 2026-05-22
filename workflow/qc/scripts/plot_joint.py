from pathlib import Path
import warnings
import pandas as pd
from pandas.api.types import is_numeric_dtype
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec
import seaborn as sns
import numpy as np
from tqdm import tqdm
import traceback
from joblib import Parallel, delayed
import scanpy as sc
from pprint import pformat
import logging
from matplotlib.colors import to_rgba
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

logging.basicConfig(level=logging.INFO)
warnings.filterwarnings(
    'ignore',
    message=r'.*pkg_resources is deprecated as an API.*',
    category=UserWarning,
)
sns.set_theme(style='ticks')
sns.set_context('paper', font_scale=1.2)

from utils.io import read_anndata
from qc_utils import parse_parameters, get_thresholds, plot_qc_joint, plot_density, QC_FLAGS


input_zarr = snakemake.input.zarr
output_joint = Path(snakemake.output.joint)
output_joint.mkdir(parents=True, exist_ok=True)
dpi = snakemake.params.get('dpi', 150)
# whether to generate density plots (optional)
plot_density_enabled = snakemake.params.get('plot_density', True)
threshold_color = snakemake.params.get('threshold_color', 'black')
threshold_color2 = snakemake.params.get('threshold_color2', 'black')
threshold_linestyle = snakemake.params.get('threshold_linestyle', '-')
threshold_linestyle2 = snakemake.params.get('threshold_linestyle2', '--')

logging.info(f'Read {input_zarr}...')
adata = read_anndata(input_zarr, obs='obs', uns='uns')
print(adata, flush=True)
adata.obs_names_make_unique()

# get parameters
file_id = snakemake.wildcards.file_id
threads = snakemake.threads
dataset, hues = parse_parameters(adata, snakemake.params, filter_hues=True)
scautoqc_metrics = snakemake.params.get('scautoqc_metrics', [])
hues = list(dict.fromkeys(hues + ['qc_status']))

logging.info(f'{hues=}')


# if no cells filtered out, save empty plots
if adata.obs.shape[0] == 0:
    logging.info('No data, skip plotting...')
    exit()


def _log1p_base(values, base):
    values = np.asarray(values, dtype=float)
    return np.log1p(values) / np.log(base)


def _thresholds_equal(threshold_a, threshold_b):
    if threshold_a is None or threshold_b is None:
        return False
    a = np.asarray(threshold_a, dtype=float)
    b = np.asarray(threshold_b, dtype=float)
    return a.shape == b.shape and np.allclose(a, b, equal_nan=True)


def create_facet_figure(df, out_file, hue, joint_title, scatter_plot_kwargs, dpi=150):
    # determine plotting parameters
    if hue is None:
        palette = None
        legend = False
    elif is_numeric_dtype(df[hue]):
        palette = 'plasma'
        legend = 'brief'
    else:
        n_unique = df[hue].nunique()
        categories = df[hue].dropna().unique()
        if n_unique > 100:
            palette = 'turbo'
            legend = False
        else:
            alpha = 0.5 if n_unique > 20 else 1.0
            if n_unique <= 10:
                base = sns.color_palette('tab10', n_colors=n_unique)
            elif n_unique <= 20:
                base = sns.color_palette('tab20', n_colors=n_unique)
            else:
                base = sc.pl.palettes.default_102
            palette = {cat: to_rgba(base[i], alpha=alpha) for i, cat in enumerate(categories)}
            legend = 'full'

    use_marginal_hue = hue is not None and hue in df.columns and not is_numeric_dtype(df[hue])
    if use_marginal_hue and df[hue].nunique() > 100:
        use_marginal_hue = False

    plot_kwargs = dict(scatter_plot_kwargs)
    plot_kwargs.update(
        palette=palette,
        legend=False,
        marginal_kwargs=dict(
            palette=palette,
            legend=False,
            # stat='density',
        ),
    )

    n_rows = len(coordinates)
    panel_size = 5.0
    fig = plt.figure(
        figsize=(2.4 * panel_size, n_rows * panel_size),
        constrained_layout=True,
    )
    outer = GridSpec(
        n_rows,
        2,
        figure=fig,
    )

    for row_idx, (x, y, log_x, log_y) in enumerate(coordinates):
        common_kwargs = dict(
            df=df,
            x=x,
            y=y,
            hue=hue,
            marginal_hue=hue if use_marginal_hue else None,
            x_threshold=updated_thresholds.get(x),
            y_threshold=updated_thresholds.get(y),
            x_threshold2=auto_thresholds.get(x),
            y_threshold2=auto_thresholds.get(y),
            threshold_color=threshold_color,
            threshold_color2=threshold_color2,
            threshold_linestyle=threshold_linestyle,
            threshold_linestyle2=threshold_linestyle2,
            title='',
            fig=fig,
        ) | plot_kwargs

        g_linear = plot_qc_joint(
            **common_kwargs,
            log_x=1,
            log_y=1,
            subplot_spec=outer[row_idx, 0],
        )
        plot_qc_joint(
            **common_kwargs,
            log_x=log_x,
            log_y=log_y,
            subplot_spec=outer[row_idx, 1],
            # sharey=g_linear.ax_joint if log_y <= 1 else None,
        )

    def _make_legend_section(title, handles, labels):
        return (
            [Line2D([], [], linestyle='none')] + handles,
            [title] + labels,
        )

    threshold_handles = [Line2D([0], [0], color=threshold_color, lw=2, linestyle=threshold_linestyle)]
    threshold_labels  = ['Applied threshold']

    has_autoqc = any(
        not _thresholds_equal(updated_thresholds.get(x), auto_thresholds.get(x))
        or not _thresholds_equal(updated_thresholds.get(y), auto_thresholds.get(y))
        for x, y, _, _ in coordinates
    )
    if has_autoqc:
        threshold_handles.append(Line2D([0], [0], color=threshold_color2, lw=2, linestyle=threshold_linestyle2))
        threshold_labels.append('sctk scAutoQC')

    hue_handles = [
        Line2D([0], [0], marker='o', linestyle='', color=to_rgba(palette[cat], alpha=1.0), markersize=5)
        for cat in palette
    ] if isinstance(palette, dict) else []
    hue_labels = list(palette.keys()) if isinstance(palette, dict) else []

    all_handles, all_labels = [], []

    h, l = _make_legend_section('Thresholds', threshold_handles, threshold_labels)
    all_handles += h
    all_labels += l

    if hue_handles:
        # Explicit spacer row between threshold and hue sections.
        all_handles += [Line2D([], [], linestyle='none')]
        all_labels += ['']
        h, l = _make_legend_section(hue, hue_handles, hue_labels)
        all_handles += h
        all_labels += l

    fig.legend(
        all_handles,
        all_labels,
        loc='outside right center',
        frameon=False,
        fontsize=10,
    )

    fig.suptitle(joint_title, fontsize=16)
    fig.savefig(out_file, dpi=dpi)
    plt.close(fig)


def call_plot(hue, scatter_plot_kwargs):
    # Handle None hue case
    if hue is None:
        plot_path = output_joint / 'main.svg'
    else:
        plot_path = output_joint / f'hue={hue}.svg'

    create_facet_figure(
        adata.obs,
        out_file=plot_path,
        hue=hue,
        joint_title=f'Joint QC for: {dataset}',
        scatter_plot_kwargs=scatter_plot_kwargs,
        dpi=dpi,
    )


updated_thresholds = get_thresholds(
    threshold_keys=scautoqc_metrics,
    autoqc_thresholds=adata.uns.get('scautoqc_ranges'),
    user_thresholds=snakemake.params.get('thresholds'),
)
auto_thresholds = get_thresholds(
    threshold_keys=scautoqc_metrics,
    autoqc_thresholds=adata.uns.get('scautoqc_ranges'),
    user_thresholds=None,
)
logging.info(f'\nupdated_thresholds={pformat(updated_thresholds)}\nauto_thresholds={pformat(auto_thresholds)}')

scatter_plot_kwargs = dict(
    s=8,
    alpha=.8,
    linewidth=0,
    rasterized=True,
)

coordinates = [
    ('n_counts', 'n_genes', 10, 2),
    ('n_genes', 'percent_mito', 2, 1),
    ('n_genes', 'percent_ribo', 2, 1),
    ('n_genes', 'percent_hb', 2, 1),
    ('n_genes', 'scrublet_score', 2, 1),
    ('n_counts', 'scrublet_score', 10, 1),
]
# filter to configured metrics
coordinates = [
    c for c in coordinates if
    all(x in adata.obs.columns for x in c[:2])
]

# reduce obs to required columns only
required_columns = [
    col for col in
    {*hues, *[c for coords in coordinates for c in coords[:2]]}
    if col in adata.obs.columns
]
adata.obs = adata.obs[required_columns].copy()


for x, y, _, _ in coordinates:
    logging.info(f'Joint QC density plot for "{x}" vs "{y}"...')

    if plot_density_enabled:
        logging.info('Plotting density...')
        plot_density(
            adata.obs,
            x,
            y,
            updated_thresholds,
            auto_thresholds=auto_thresholds,
            datashader_kwargs=dict(sigma=5, outlier_boost=2),
        )
        plt.tight_layout()
        plt.savefig(output_joint / f'{x}_vs_{y}_density.svg', bbox_inches='tight', dpi=dpi)
        plt.close('all')

batch_size = max(1, len(hues) // threads)
results = Parallel(n_jobs=threads, batch_size=batch_size, return_as='generator')(
    delayed(call_plot)(
        hue=hue,
        scatter_plot_kwargs=scatter_plot_kwargs,
    ) for hue in hues
)
# Consume generator to trigger execution and display progress
for _ in tqdm(
    results,
    total=len(hues),
    desc=f'Plotting all coordinates faceted by hue using {threads=}, {batch_size=}'
):
    pass
