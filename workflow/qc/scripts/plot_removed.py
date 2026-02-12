from pathlib import Path
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
from pprint import pformat
import logging
from tqdm import tqdm
import traceback
from joblib import Parallel, delayed

logging.basicConfig(level=logging.INFO)

from utils.io import read_anndata
from qc_utils import parse_parameters, QC_FLAGS


def get_fraction_removed(df, group_key, key='qc_status'):
    grouped_counts = df.groupby([group_key, key], observed=False).size().reset_index(name='Counts')
    grouped_frac = grouped_counts.pivot(index=group_key, columns=key, values='Counts')
    grouped_frac['fraction_removed'] = (grouped_frac['failed'] / grouped_frac.sum(axis=1)).round(2)
    return grouped_frac


input_zarr = snakemake.input.zarr
output_plots = snakemake.output.plots
threads = snakemake.threads

output_plots = Path(output_plots)
output_plots.mkdir(parents=True, exist_ok=True)

logging.info(f'Read {input_zarr}...')
adata = read_anndata(input_zarr, obs='obs', uns='uns')

# If no cells filtered out, save empty plots
if adata.obs.shape[0] == 0:
    logging.info('Empty data, skipping plots...')
    exit(0)

# get parameters
file_id = snakemake.wildcards.file_id
dataset, hues = parse_parameters(adata, snakemake.params, filter_hues=True)
scautoqc_metrics = snakemake.params.get('scautoqc_metrics', QC_FLAGS)

logging.info(f'{hues=}')

logging.info('Plot removed cells...')
plt.figure(figsize=(4, 5))
plt.grid(False)
sns.countplot(
    x='qc_status',
    data=adata.obs,
    hue='qc_status',
)
ax = plt.gca()
for pos in ['right', 'top']: 
    ax.spines[pos].set_visible(False)
for container in ax.containers:
    ax.bar_label(container)
plt.xlabel('Cell QC Status')
plt.ylabel('Count')
plt.title(f'Counts of cells QC\'d\n{dataset}')
plt.tight_layout()
plt.savefig(output_plots / 'cells_passed_all.png', bbox_inches='tight')
plt.close()


def plot_composition(df, group_key, plot_dir):
    n_hues = df[group_key].nunique()

    # Get counts by group and QC status
    counts = df.groupby([group_key, 'qc_status'], observed=False).size().unstack(fill_value=0)
    
    # Sort by fraction failed (descending) for consistent ordering
    fraction_failed = (counts.get('failed', 0) / counts.sum(axis=1)).sort_values(ascending=False)
    order = fraction_failed.index
    counts_ordered = counts.loc[order]
    
    # Calculate proportions and totals
    proportions = counts_ordered.div(counts_ordered.sum(axis=1), axis=0) * 100
    total_counts = counts_ordered.sum(axis=1)
    
    # Create JointGrid-like figure with space for legend
    from matplotlib.gridspec import GridSpec
    # Calculate dynamic height based on number of categories with better scaling
    # Use 0.3 inches per category with a minimum of 6 inches
    fig_height = max(6, n_hues * 0.3)
    fig = plt.figure(figsize=(11, fig_height))
    gs = GridSpec(
        nrows=1,
        ncols=3,
        width_ratios=[0.85, 0.1, 0.05],
        wspace=0.1,
        hspace=0,
        top=0.95,
        bottom=0.02
    )
    
    ax_main = fig.add_subplot(gs[0, 0])
    ax_margin = fig.add_subplot(gs[0, 1], sharey=ax_main)
    
    # Main plot: Relative abundance stacked barplot
    bottom = pd.Series(0, index=proportions.index)
    bars_per_group = {idx: [] for idx in proportions.index}
    
    for status in proportions.columns:
        widths = proportions[status]
        bars = ax_main.barh(
            proportions.index,
            widths,
            left=bottom,
            edgecolor='white',
            linewidth=0.5,
            label=status,
            alpha=0.9,
        )
        bottom += widths
        for idx, bar in zip(proportions.index, bars):
            bars_per_group[idx].append((bar, status))
    
    ax_main.set_xlabel('% Cells')
    ax_main.set_ylabel(group_key)
    ax_main.set_xlim([0, 100])
    
    # Add detailed labels for main plot
    fontsize = 8
    fontweight = 'bold'
    # Ensure the figure is drawn before querying the renderer/text extents
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    for idx, bar_status_pairs in bars_per_group.items():
        y_pos = bar_status_pairs[0][0].get_y() + bar_status_pairs[0][0].get_height() / 2

        for bar, status in bar_status_pairs:
            count_value = counts_ordered.loc[idx, status]
            label_text = str(int(count_value))

            temp_text = ax_main.text(0, 0, label_text, fontsize=fontsize, fontweight=fontweight)
            text_width_px = temp_text.get_window_extent(renderer=renderer).width
            temp_text.remove()

            left_px = ax_main.transData.transform((bar.get_x(), 0))[0]
            right_px = ax_main.transData.transform((bar.get_x() + bar.get_width(), 0))[0]
            bar_width_px = right_px - left_px

            # Place label inside bar if it fits, otherwise outside
            if bar_width_px >= text_width_px + 3:
                x_pos = bar.get_x() + bar.get_width() / 2
                color = 'white'
                ha = 'center'
            else:
                x_pos = bar.get_x() + bar.get_width() + proportions.sum(1).max() / 100
                color = bar.get_facecolor()
                ha = 'left'
            
            ax_main.text(
                x_pos,
                y_pos,
                label_text,
                va='center',
                ha=ha,
                fontsize=fontsize,
                fontweight=fontweight,
                color=color,
            )
    
    # Margin plot: Total cell counts as histogram with annotations
    ax_margin.barh(
        total_counts.index,
        total_counts.values,
        edgecolor='white',
        linewidth=0.5,
        alpha=0.9,
        color='#555555',  # Dark gray
    )
    
    ax_margin.tick_params(axis='y', labelleft=False, left=False)
    ax_margin.tick_params(axis='x', labelbottom=True, labelsize=8, rotation=270)
    ax_margin.set_xlabel('# Cells', fontsize=9, labelpad=10)
    ax_margin.locator_params(axis='x', nbins=5)
    ax_margin.ticklabel_format(axis='x', style='plain')
    
    # Format axes
    for pos in ['right', 'top']:
        ax_main.spines[pos].set_visible(False)
        ax_margin.spines[pos].set_visible(False)
    
    # Add legend at figure level in the reserved space
    handles, labels = ax_main.get_legend_handles_labels()
    fig.legend(handles, labels, title='QC Status', loc='center left', bbox_to_anchor=(0.88, 0.5), frameon=False, fontsize=9)
    
    fig.suptitle(f'Cells that passed QC\n{dataset=}, {file_id=}', fontsize=12, y=0.98)
    fig.tight_layout()
    fig.savefig(plot_dir / f'by={group_key}.png', bbox_inches='tight', dpi=100)
    plt.close()


logging.info('Plot compositions...')

def safe_call_plot(*args, **kwargs):
    try:
        plot_composition(*args, **kwargs)
        return None
    except Exception as e:
        logging.error(f"Error in plot job: {e}")
        traceback.print_exc()
        return e

# run in parallel using joblib; results is a list of return values (None or Exception)
results = Parallel(n_jobs=threads)(
    delayed(safe_call_plot)(
        adata.obs,
        group_key=g,
        plot_dir=output_plots
    ) for g in tqdm(hues)
)

# log any errors
errors = [(g, r) for g, r in zip(hues, results) if r is not None]
if errors:
    logging.error(f"{len(errors)} group(s) failed during plotting.")
    for g, e in errors:
        logging.error(f"Group {g} failed: {e}")

logging.info('Plot violin plots per QC metric...')
n_cols = len(scautoqc_metrics)
fig, axes = plt.subplots(1, n_cols, figsize=(4 * n_cols, 6))
plt.grid(False)

for i, qc_metric in enumerate(scautoqc_metrics):
    sns.violinplot(
        data=adata.obs,
        x='qc_status',
        y=qc_metric,
        hue='qc_status',
        fill=False,
        inner='quartile',
        legend=False,
        ax=axes[i]
    )
    axes[i].set_xlabel('QC status')
    axes[i].set_ylabel(qc_metric)
    axes[i].set_title(f'{qc_metric} distribution')
    for pos in ['right', 'top']: 
        axes[i].spines[pos].set_visible(False) 

plt.suptitle(f'Cells that passed QC\n{dataset}', fontsize=12)
plt.tight_layout()
plt.savefig(output_plots / 'per_metric_violin.png', bbox_inches='tight')
