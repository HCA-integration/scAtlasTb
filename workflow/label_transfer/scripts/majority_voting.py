import logging
logging.basicConfig(level=logging.INFO)
from pathlib import Path
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from pprint import pformat

from utils.io import read_anndata, write_zarr_linked
from _utils import get_majority_reference, get_majority_consensus


input_file = snakemake.input.zarr
output_file = snakemake.output.zarr
output_plots = Path(snakemake.output.plots)
output_plots.mkdir(parents=True, exist_ok=True)
majority_reference = snakemake.params.get('majority_reference')
majority_consensus = snakemake.params.get('majority_consensus')

logging.info(f'majority_reference: {majority_reference}')
logging.info(f'majority_consensus: {majority_consensus}')

logging.info('Read adata...')
adata = read_anndata(input_file, obs='obs')

if majority_reference is not None:
    logging.info(f'Compute majority reference for:\n{pformat(majority_reference)}')
    
    reference_key = majority_reference['reference_key']
    query_key = majority_reference['query_key']
    
    adata.obs[reference_key] = adata.obs[reference_key].astype(str).replace('nan', float('nan'))
    
    adata.obs['majority_reference'] = get_majority_reference(
        adata.obs,
        reference_key=reference_key,
        query_key=query_key,
        **majority_reference.get('crosstab_kwargs', {})
    )

if majority_consensus is not None and majority_consensus.get('columns'):
    logging.info(f'Compute majority consensus for:\n{pformat(majority_consensus)}')
    
    new_key = 'majority_consensus'
    maj_df = get_majority_consensus(
        adata.obs,
        column_patterns=majority_consensus['columns'],
        new_key=new_key,
    )
    adata.obs[maj_df.columns] = maj_df
    
    # get majority consensus agreement stasts
    counts = maj_df.value_counts(
        subset=[new_key, f'{new_key}_agreement'],
        sort=False,
        dropna=False
    ).reset_index(name='count')
    counts[f'total_cells_{new_key}_label'] = counts.groupby(new_key, observed=True)['count'].transform('sum')
    counts[f'frac_cells_in_{new_key}_label'] = counts['count'] / counts[f'total_cells_{new_key}_label']
    print(counts, flush=True)
    counts.to_csv(output_plots / 'majority_consensus_agreement.tsv', sep='\t', index=False)
    
    # plot majority consensus stats
    plot_df = pd.crosstab(
        maj_df[new_key],
        maj_df[f'{new_key}_low_agreement'],
        dropna=False,
        normalize='index'
    )
    
    if True not in plot_df.columns:
        plot_df[True] = 0
    elif False not in plot_df.columns:
        plot_df[False] = 0
    
    ax = plot_df.sort_values(True).plot(
        kind='barh',
        stacked=True,
        color=sns.color_palette("colorblind").as_hex(),
        figsize=(5, 5 + 0.1 * plot_df.shape[0]),
    )
    ax.legend(title='Low agreement', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.title(f'Fraction of cells with low agreement in {new_key} labels')
    plt.grid(False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.savefig(output_plots / 'majority_consensus_frac.png', bbox_inches='tight')


logging.info(f'Write to {output_file}...')
write_zarr_linked(
    adata=adata,
    in_dir=input_file,
    out_dir=output_file,
    files_to_keep=['obs']
)