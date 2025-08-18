from pathlib import Path
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
from pprint import pformat
import logging
from tqdm import tqdm
import traceback
import concurrent.futures

logging.basicConfig(level=logging.INFO)

from utils.io import read_anndata
from qc_utils import parse_parameters


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
dataset, groups = parse_parameters(adata, snakemake.params, filter_hues=True)
threshold_keys = ['n_counts', 'n_genes', 'percent_mito'] 

logging.info('Plot removed cells...')
plt.figure(figsize=(4, 5))
plt.grid(False)
sns.countplot(
    x='qc_status',
    data=adata.obs,
    hue='qc_status',
    palette='muted', # 'Set2'
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


def plot_composition(adata, group_key, plot_dir):
    n_groups = adata.obs[group_key].nunique()
    if n_groups > 100:
        logging.info(f'Group {group_key} has too many unique values, skipping...')
        return
    
    grouped_frac = get_fraction_removed(adata.obs, group_key=group_key, key='qc_status')
    order = grouped_frac.sort_values('fraction_removed', ascending=False).index
    adata.obs[group_key] = pd.Categorical(adata.obs[group_key], categories=order)

    f, (ax1, ax2) = plt.subplots(1, 2, sharey=True, figsize=(14, 5 * (1 + n_groups/50)))

    counts = adata.obs.groupby([group_key, 'qc_status'], observed=False).size().unstack(fill_value=0)
    bottom = pd.Series(0, index=counts.index)
    bars_per_group = {idx: [] for idx in counts.index}

    # stacked barplot for QC status
    for status in counts.columns:
        widths = counts[status]
        bars = ax1.barh(
            counts.index,
            widths,
            left=bottom,
            edgecolor='white',
            linewidth=0,
            label=status,
            alpha=0.9,
        )
        bottom += widths
        for idx, bar in zip(counts.index, bars):
            bars_per_group[idx].append(bar)

        ax1.legend(
            title="QC Status",
            loc="best",
            frameon=False,
        )

    sns.barplot(
        data=grouped_frac,
        x='fraction_removed',
        y=group_key,
        order=order,
        ax=ax2,
        alpha=0.9,
    )
    if n_groups < 50:
        for container in ax2.containers:
            ax2.bar_label(container)
    
    for pos in ['right', 'top']: 
        ax1.spines[pos].set_visible(False)
        ax2.spines[pos].set_visible(False)
    f.suptitle(f'Cells that passed QC\n{dataset}')
    f.tight_layout()
    f.canvas.draw()

    # Add labels per bar for ax1 plot
    fontsize = 9
    fontweight = 'bold'
    renderer = ax1.figure.canvas.get_renderer()
    for idx, bars in bars_per_group.items(): # iterate groups
        bbox = ax1.get_window_extent()
        x_min, x_max = ax1.get_xlim()
        offset = counts.sum(1).max() / 100
        
        total_width = sum(bar.get_width() for bar in bars)
        outside_x = total_width + offset
        y_pos = bars[0].get_y() + bars[0].get_height() / 2

        for bar in bars: # iterate bar in stack
            label_text = str(int(bar.get_width()))

            # Get text width in pixels
            temp_text = ax1.text(0, 0, label_text, fontsize=fontsize, fontweight=fontweight)
            text_width_px = temp_text.get_window_extent(renderer=renderer).width
            temp_text.remove()

            # Convert bar edges to pixel coordinates
            left_px = ax1.transData.transform((bar.get_x(), 0))[0]
            right_px = ax1.transData.transform((bar.get_x() + bar.get_width(), 0))[0]
            bar_width_px = right_px - left_px

            if bar_width_px >= text_width_px + 3:
                # Inside white label at right edge of the bar
                ax1.text(
                    bar.get_x() + bar.get_width() - offset,
                    y_pos,
                    label_text,
                    va='center',
                    ha='right',
                    fontsize=fontsize,
                    fontweight=fontweight,
                    color='white',
                )
            else:
                # Outside label in bar's own color at end of full stack
                ax1.text(
                    outside_x + offset,
                    y_pos,
                    label_text,
                    va='center',
                    ha='left',
                    fontsize=fontsize,
                    fontweight=fontweight,
                    color=bar.get_facecolor(),
                )
                # add spacing for next bar in group
                outside_x += (text_width_px * (x_max - x_min) / bbox.width) + offset
    
    
    f.savefig(plot_dir / f'by={group_key}.png', bbox_inches='tight')
    plt.close()


logging.info('Plot compositions...')
# for group_key in tqdm(groups):
#     plot_composition(adata, group_key, plot_dir=output_plots)

with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
    futures = [
        executor.submit(
            plot_composition,
            adata,
            group_key=group_key,
            plot_dir=output_plots
        ) for group_key in groups
    ]
    for future in tqdm(concurrent.futures.as_completed(futures), total=len(futures)):
        try:
            future.result()
        except Exception as e:
            logging.error(f"Exception occurred: {e}")
            traceback.print_exc()

logging.info('Plot violin plots per QC metric...')
n_cols = len(threshold_keys)
fig, axes = plt.subplots(1, n_cols, figsize=(4 * n_cols, 6))
plt.grid(False)

for i, qc_metric in enumerate(threshold_keys):
    sns.violinplot(
        data=adata.obs,
        x='qc_status',
        y=qc_metric,
        hue='qc_status',
        palette='muted', # 'Set2',
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
