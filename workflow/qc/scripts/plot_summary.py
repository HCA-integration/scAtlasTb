from pathlib import Path
import re
import warnings
import logging
from tqdm import tqdm

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D

from qc_utils import QC_FLAGS


warnings.simplefilter("ignore", UserWarning)
logging.basicConfig(level=logging.INFO)
sns.set_theme(style='white', context='paper')


def infer_file_id(path, fallback):
    match = re.search(r'file_id~([^/]+)', str(path))
    if match:
        return match.group(1)
    return fallback


def parse_file_id_parts(file_id):
    file_id = str(file_id)

    # Keep only the part after the last ':' before the first '--'.
    prefix = file_id.split('--', 1)[0]
    file_name = prefix.rsplit(':', 1)[-1]

    if file_name == '':
        file_name = file_id

    split_data_value = ''
    for token in file_id.split('--')[1:]:
        if token.startswith('split_data_value='):
            split_data_value = token.split('=', 1)[1] or ''
            break

    return file_name, split_data_value


def safe_name(name):
    return re.sub(r'[^A-Za-z0-9_.-]+', '_', str(name)).strip('_')


def _extract_bounds(df, study, lineage, lineage_key, thresh_type, metrics):
    """Helper to cleanly extract (min, max) threshold tuples for a framework."""
    match = df[(df["study"] == study) & (df[lineage_key] == lineage) & (df["threshold_type"] == thresh_type)]
    if match.empty:
        return {m: (None, None) for m in metrics}

    row = match.iloc[0]
    return {m: (row.get(f"{m}_min"), row.get(f"{m}_max")) for m in metrics}


def ridge_plot(
    df,
    row_group,
    col_group,
    metric,
    out_file=None,
    dpi=100,
    thresholds_df=None,
    threshold_color='black',
    hspace=-0.5,
    linewidth=1.5,
):
    THRESHOLD_TYPES = [
        ('sctk_autoqc', '--', 'sctk scAutoQC'),
        ('updated',     '-',  'Applied threshold'),
        ('alternative', ':',  'Alternative threshold'),
    ]

    # ── Grid setup ───────────────────────────────────────────────────────────
    sns.set_theme(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})

    grid = sns.FacetGrid(
        df,
        row=row_group,
        col=col_group,
        hue=row_group,
        aspect=6,
        height=0.7,
        palette=sns.cubehelix_palette(df[row_group].nunique(), rot=-.25, light=.7),
        margin_titles=False,
    )

    # ── KDE layers ───────────────────────────────────────────────────────────
    kde_kwargs = dict(bw_adjust=1, clip_on=False)
    grid.map(sns.kdeplot, metric, fill=True, alpha=0.9, **kde_kwargs)
    grid.map(sns.kdeplot, metric, color='white', linewidth=linewidth, **kde_kwargs)
    grid.refline(y=0, linewidth=linewidth, linestyle="-", color=None, clip_on=False)

    # ── Layout ───────────────────────────────────────────────────────────────
    grid.set(yticks=[], ylabel='')
    grid.set_xlabels(metric)
    grid.despine(bottom=True, left=True)
    grid.figure.subplots_adjust(left=0.18, right=0.98, top=0.98, hspace=hspace, wspace=0.12)

    # ── Column titles (workaround for seaborn bug) ────────────────────────────
    for ax in grid.axes.flat:
        ax.set_title("")
    if df[col_group].nunique() > 1:
        for ax, title in zip(grid.axes[0], grid.col_names):
            ax.set_title(title.split(" | ")[-1].split(" = ")[-1])

    # ── Row labels ───────────────────────────────────────────────────────────
    for row_idx, row_name in enumerate(grid.row_names):
        grid.axes[row_idx, 0].text(
            -0.08, 0.04, str(row_name),
            transform=grid.axes[row_idx, 0].transAxes,
            ha='right', va='bottom', fontsize=12
        )

    if thresholds_df is None or thresholds_df.empty:
        return grid

    # ── Threshold lines ───────────────────────────────────────────────────────
    base_metric = metric.split('_', 1)[1] if metric.startswith('log1p_') else metric
    log_transform = metric.startswith('log1p_')
    study_first = (row_group == 'file_name' and col_group == 'split_data_value')
    seen_types = set()

    # draw threshold lines in a separate overlay axis to ensure they are visible above the KDE layers
    overlay = grid.figure.add_axes([0, 0, 1, 1], facecolor='none', zorder=10)
    overlay.set_xlim(0, 1)
    overlay.set_ylim(0, 1)
    overlay.axis('off')

    to_display = grid.figure.transFigure.inverted().transform

    def to_fig(ax, x_data, y_axes):
        x_fig, _ = to_display(ax.transData.transform((x_data, 0)))
        _, y_fig = to_display(ax.transAxes.transform((0, y_axes)))
        return x_fig, y_fig

    def draw_threshold(ax, study_name, lineage_name, thresh_type, linestyle):
        bounds = _extract_bounds(thresholds_df, study_name, lineage_name, 'lineage', thresh_type, [base_metric])
        vmin, vmax = bounds.get(base_metric, (None, None))
        drawn = False
        for v in (vmin, vmax):
            if pd.notna(v):
                v = np.log1p(float(v)) if log_transform else float(v)
                xlim = ax.get_xlim()
                # Only draw threshold if it falls within the current plot x-limits
                if v < xlim[0] or v > xlim[1]:
                    continue
                x, y_bot = to_fig(ax, v, 0)
                _, y_top = to_fig(ax, v, 1 - abs(hspace))
                overlay.vlines(x, y_bot, y_top, colors=threshold_color, linestyles=linestyle, linewidth=1)
                drawn = True
        return drawn

    for i, r_name in enumerate(grid.row_names):
        for j, c_name in enumerate(grid.col_names):
            study, lineage = (str(r_name), str(c_name)) if study_first else (str(c_name), str(r_name))
            for thresh_type, linestyle, _ in THRESHOLD_TYPES:
                if draw_threshold(grid.axes[i, j], study, lineage, thresh_type, linestyle):
                    seen_types.add(thresh_type)

    # ── Threshold legend ──────────────────────────────────────────────────────
    legend_handles = [
        Line2D([], [], color=threshold_color, linestyle=ls, lw=linewidth, label=label)
        for thresh_type, ls, label in THRESHOLD_TYPES
        if thresh_type in seen_types
    ]
    if legend_handles:
        grid.figure.legend(
            handles=legend_handles,
            labels=[h.get_label() for h in legend_handles],
            loc='center left',
            bbox_to_anchor=(1.02, 0.5),
            bbox_transform=grid.figure.transFigure,
            frameon=False,
            fontsize=10,
        )
    
    if out_file:
        grid.savefig(out_file, dpi=dpi, bbox_inches='tight')
        plt.close(grid.figure)


def plot_heatmap(data_df, title, filename, out_dir, annot=False, dpi=200):
    """Write a heatmap SVG; kept near `ridge_plot` for shared presentation defaults."""
    plot_df = data_df.copy()
    transposed = plot_df.shape[0] < plot_df.shape[1]
    if transposed:
        plot_df = plot_df.T

    figsize = (3, 5) if not transposed else (5, 3)
    fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
    sns.heatmap(plot_df, cmap='Blues', annot=annot, ax=ax)
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position('top')
    ax.set_xlabel('')
    ax.set_ylabel('')
    plt.xticks(rotation=90)
    fig.suptitle(title)
    out_file = out_dir / filename
    fig.savefig(out_file, dpi=dpi, bbox_inches='tight')
    plt.close(fig)
    return out_file.name


# removed wrapper: heatmap creation inlined where needed below


input_files = list(snakemake.input.qc_metrics)
threshold_files = list(snakemake.input.tsv)
output_dir = Path(snakemake.output.plots)
output_dir.mkdir(parents=True, exist_ok=True)
dpi = snakemake.params.get('dpi', 150)
qc_metric_columns = snakemake.params.get('scautoqc_metrics', QC_FLAGS)
qc_metric_columns += [
    f'log1p_{col}' for col in qc_metric_columns
    if col in ['n_counts', 'n_genes']
]
threshold_color = snakemake.params.get('threshold_color', '#E66101')
# Hardcoded linestyles for thresholds
threshold_linestyle = '-'
threshold_linestyle2 = '--'
output_ridge = output_dir / 'ridge_plots'
output_ridge.mkdir(parents=True, exist_ok=True)
written_files = []


frames = []
if len(input_files) == 0:
    logging.info('No qc_metrics input files found. Skipping ridge plots.')

for idx, path in tqdm(
    enumerate(input_files),
    total=len(input_files),
    desc='Reading qc_metrics files',
):
    try:
        df = pd.read_parquet(path)
    except Exception as exc:
        logging.warning(f'Could not read qc metrics parquet: {path} ({exc})')
        continue

    if df.shape[0] == 0:
        continue
    file_id = infer_file_id(path, fallback=f'file_{idx + 1}')
    file_name, split_data_value = parse_file_id_parts(file_id)
    df = df.copy()
    df['file_id'] = file_id
    df['file_name'] = file_name
    df['split_data_value'] = split_data_value
    frames.append(df)

qc_df = pd.concat(frames, ignore_index=True) if len(frames) > 0 else None
if qc_df is None:
    logging.info('All qc_metrics files are empty or unreadable. Skipping ridge plots.')

threshold_frames = []
for path in threshold_files:
    try:
        thresholds = pd.read_table(path)
    except Exception as exc:
        logging.warning(f'Could not read threshold TSV: {path} ({exc})')
        continue

    if thresholds.shape[0] == 0:
        continue
    thresholds = thresholds.copy()
    thresholds['file_id'] = thresholds['file_id'].astype(str)
    parsed = thresholds['file_id'].apply(parse_file_id_parts)
    thresholds['study'] = parsed.apply(lambda x: x[0])
    thresholds['lineage'] = parsed.apply(lambda x: x[1])
    threshold_frames.append(thresholds)

# If threshold files were provided, build thresholds DataFrame now
if len(threshold_frames) > 0:
    thresholds_df = pd.concat(threshold_frames, ignore_index=True)
else:
    thresholds_df = None

if qc_df is not None:
    metric_columns = [
        col for col in qc_metric_columns
        if col in qc_df.columns and pd.api.types.is_numeric_dtype(qc_df[col])
    ]

    if len(metric_columns) == 0:
        logging.info('No numeric QC metric columns found. Skipping ridge plots.')
    else:
        logging.info(f'Found {len(metric_columns)} numeric QC metric columns for plotting: {metric_columns}')
        for metric in tqdm(metric_columns, desc='Creating ridge plots'):
            metric_df = qc_df[['file_id', 'file_name', 'split_data_value', metric]].copy()
            metric_df[metric] = pd.to_numeric(metric_df[metric], errors='coerce')
            metric_df = metric_df[np.isfinite(metric_df[metric])]

            if metric_df.shape[0] < 2 or metric_df['file_name'].nunique() == 0:
                continue
            if metric_df[metric].nunique() < 2:
                continue

            n_file_names = metric_df['file_name'].nunique()
            n_split_values = metric_df['split_data_value'].nunique()
            if n_file_names < n_split_values:
                row_group = 'split_data_value'
                col_group = 'file_name'
            else:
                row_group = 'file_name'
                col_group = 'split_data_value'

            out_file = output_ridge / f'{safe_name(metric)}.svg'
            ridge_plot(
                metric_df,
                row_group=row_group,
                col_group=col_group,
                metric=metric,
                out_file=out_file,
                dpi=dpi,
                thresholds_df=thresholds_df,
                threshold_color=threshold_color,
            )
            written_files.append(out_file.name)

if len(written_files) > 0:
    logging.info(f'Created {len(written_files)} ridge plots in {output_ridge}')
elif qc_df is not None:
    logging.info('No ridge plots were created after metric filtering.')


logging.info('Creating removed-cell heatmaps...')
if thresholds_df is not None:
    thresholds = thresholds_df.query('threshold_type == "updated"')
    if thresholds.shape[0] == 0:
        logging.info('No updated thresholds found. Skipping removed-cell heatmaps.')
    else:
        thresholds = thresholds.drop_duplicates(subset=['study', 'lineage'])
        removed_df = thresholds.pivot(index=['study'], columns='lineage', values='n_removed')
        removed_frac_df = thresholds.pivot(index=['study'], columns='lineage', values='removed_frac')

        plot_heatmap(
            removed_df,
            'Cells removed from semi-automated QC',
            'n_removed_heatmap.svg',
            output_dir,
            annot=False,
            dpi=dpi
        )
        
        plot_heatmap(
            removed_frac_df,
            'Fraction removed from semi-automated QC',
            'removed_frac_heatmap.svg',
            output_dir,
            annot=True,
            dpi=dpi
        )
        logging.info(f'Created removed-cell heatmaps in {output_dir}')
else:
    logging.info('No thresholds TSV rows found, skipping removed-cell heatmaps.')
