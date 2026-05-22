from pathlib import Path
import re
import warnings
import logging
from tqdm import tqdm

import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt


warnings.simplefilter("ignore", UserWarning)
logging.basicConfig(level=logging.INFO)
sns.set_theme(style='white', context='paper')


def infer_file_id(path, fallback):
    match = re.search(r'file_id~([^/]+)', str(path))
    if match:
        return match.group(1)
    return fallback


def safe_name(name):
    return re.sub(r'[^A-Za-z0-9_.-]+', '_', str(name)).strip('_')


def ridge_plot(df, group, metric, out_file=None, dpi=100):
    sns.set_theme(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})

    def add_label(x, color, label, **kwargs):
        ax = plt.gca()
        ax.text(
            0,
            0.1,
            label,
            transform=ax.transAxes,
            fontsize=10,
            color=color,
            ha="right",
            va="center",
            clip_on=False,
        )
    
    grid = sns.FacetGrid(
        df,
        row=group,
        hue=group,
        aspect=15,
        height=0.65,
        palette=sns.cubehelix_palette(df[group].nunique(), rot=-.25, light=.7),
    )

    grid.map(sns.kdeplot, metric, fill=True, alpha=1, linewidth=1, bw_adjust=0.8, clip_on=False)
    grid.map(sns.kdeplot, metric, color='white', linewidth=1.5, bw_adjust=0.8, clip_on=False)
    grid.refline(y=0, linewidth=1.5, linestyle="-", color=None, clip_on=False)

    grid.map(add_label, metric)

    # Remove axes details that don't play well with overlap
    grid.set_titles('')
    grid.set(yticks=[], ylabel='')
    grid.despine(bottom=True, left=True)

    grid.set_xlabels(metric)

    grid.figure.subplots_adjust(hspace=-0.5)
    grid.figure.suptitle(f'QC metric ridge plot: {metric}', y=1.02)
    
    if out_file:
        grid.savefig(out_file, dpi=dpi, bbox_inches='tight')
        plt.close(grid.figure)


input_files = list(snakemake.input.qc_metrics)
output_dir = Path(snakemake.output.plots)
output_dir.mkdir(parents=True, exist_ok=True)
dpi = snakemake.params.get('dpi', 200)

frames = []
for idx, path in tqdm(
    enumerate(input_files),
    total=len(input_files),
    desc='Reading qc_metrics files',
):
    df = pd.read_parquet(path)
    if df.shape[0] == 0:
        continue
    file_id = infer_file_id(path, fallback=f'file_{idx + 1}')
    df = df.copy()
    df['file_id'] = file_id
    print(file_id, path)
    frames.append(df)

if len(frames) == 0:
    logging.info('All qc_metrics files are empty, skipping summary ridge plots.')
    (output_dir / 'empty_qc_metrics.txt').write_text('All qc_metrics files are empty.\n')
    exit(0)

qc_df = pd.concat(frames, ignore_index=True)

metric_columns = [
    col
    for col in qc_df.columns
    if col != 'file_id' and pd.api.types.is_numeric_dtype(qc_df[col])
    and col != 'cell_passed_qc'
]

assert len(metric_columns) > 0, 'No numeric QC metric columns found in qc_metrics files.'

# ridge plots
output_ridge = output_dir / 'ridge_plots'
output_ridge.mkdir(parents=True, exist_ok=True)

written_files = []
for metric in tqdm(metric_columns, desc='Creating ridge plots'):
    metric_df = qc_df[['file_id', metric]].copy()
    metric_df[metric] = pd.to_numeric(metric_df[metric], errors='coerce')
    metric_df = metric_df[np.isfinite(metric_df[metric])]

    if metric_df.shape[0] < 2 or metric_df['file_id'].nunique() == 0:
        continue
    if metric_df[metric].nunique() < 2:
        continue

    out_file = output_ridge / f'{safe_name(metric)}.svg'
    ridge_plot(
        metric_df,
        group='file_id',
        metric=metric,
        out_file=out_file, dpi=dpi
    )
    written_files.append(out_file.name)

logging.info(f'Created {len(written_files)} ridge plots in {output_ridge}')
