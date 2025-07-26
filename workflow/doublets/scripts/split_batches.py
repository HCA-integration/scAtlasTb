import faulthandler
faulthandler.enable()
import logging
logging.basicConfig(level=logging.INFO)
from pathlib import Path
import pandas as pd

from utils.io import read_anndata, link_zarr

input_zarr = snakemake.input.zarr
output_dir = snakemake.output[0]
batch_key = snakemake.params.get('batch_key')
max_cells_per_group = snakemake.params.get('max_cells_per_group', 100_000)

def group_batches(adata, batch_key, max_cells_per_group=100_000):
    value_counts = adata.obs[batch_key].value_counts(dropna=False)
    logging.info(f'Value counts: {value_counts}')

    # Group batches such that each group has <= max_cells_per_group cells
    groups = []
    current_group = []
    current_count = 0

    for batch, count in value_counts.items():
        if current_count + count > max_cells_per_group and current_group:
            groups.append(current_group)
            current_group = [batch]
            current_count = count
        else:
            current_group.append(batch)
            current_count += count

    if current_group:
        groups.append(current_group)

    return groups



adata = read_anndata(input_zarr, obs='obs')

output_dir = Path(output_dir)
output_dir.mkdir(exist_ok=True, parents=True)

if batch_key not in adata.obs.columns:
    logging.info(f'No batch key found in obs columns, setting dummy batch.')
    Path(output_dir / 'no_batch.txt').write_text('no_batch')
else:
    groups = group_batches(adata, batch_key, max_cells_per_group=max_cells_per_group)
    value_counts = adata.obs[batch_key].value_counts(dropna=False)
    for i, group in enumerate(groups):
        group_file = output_dir / f"group_{i + 1}.txt"
        logging.info(f'Writing group {i + 1} to {group_file}: {group}')
        Path(group_file).write_text('\n'.join(group) + '\n')
