import re
import pandas as pd
import warnings

warnings.filterwarnings('ignore')

input_metrics = snakemake.input.metrics
input_benchmark = snakemake.input.benchmark
out_tsv = snakemake.output.tsv
extra_columns = snakemake.output.extra_columns
expanded_wildcards = snakemake.params['wildcards']
# expanded_wildcards = pd.DataFrame(expanded_wildcards)
max_len = int(snakemake.params.get('max_file_name_len', 100))

metrics_df = pd.concat(
    [pd.read_table(file) for file in input_metrics],
    ignore_index=True,
).reset_index(drop=True)
print(metrics_df, flush=True)

# set extra columns
columns_to_ignore = ['score', 'metric_type', 'metric_name']
ex_columns = set(
    [
        col for col in metrics_df.columns
        if col not in expanded_wildcards.columns.tolist()+columns_to_ignore
    ]
)

benchmark_df = pd.concat(
    [pd.read_table(file) for file in input_benchmark],
    ignore_index=True,
).reset_index(drop=True)

print(benchmark_df, flush=True)

benchmark_df = pd.concat([expanded_wildcards, benchmark_df], axis=1)
metrics_df = metrics_df.merge(benchmark_df).drop_duplicates()

# parse any extra entries in file_id
expanded_file_ids = metrics_df['file_id'].str.split('--', expand=True)


def check_existing(df, row, key):
    return key in df.columns \
        and isinstance(df.loc[row, key], str) \
        and df.loc[row, key] != ''


for row, _dict in expanded_file_ids.to_dict('index').items():
    # skip file_id if overwrite_file_id is True
    if not metrics_df.loc[row, 'overwrite_file_id']:
        key = 'file_name'
        metrics_df.loc[row, key] = metrics_df.loc[row, 'file_id']
        ex_columns.add(key)
        continue
    
    for _, value in _dict.items():
        if value is None:
            continue
        
        if '=' in value: # split key, value
            key, value = value.split('=', maxsplit=1)
        elif ':' in value:
            key, value = 'file_name', value.split(':')[-1]
        else:
            key = 'file_name'

        if check_existing(metrics_df, row, key):
            if key == 'file_name' and metrics_df.loc[row, key] != value:
                value = metrics_df.loc[row, key] + '--' + value
            else:
                # don't overwrite existing values
                continue
        
        # add extra metadata
        metrics_df.loc[row, key] = value
        ex_columns.add(key)

if 'file_name' in metrics_df.columns:
    # truncate file_name if too long
    metrics_df['file_name'] = metrics_df['file_name'].apply(
        lambda x: x[:max_len] + '...' if isinstance(x, str) and len(x) > max_len else x
    )

# rename metric names
metrics_df['metric'] = metrics_df['metric_name']
del metrics_df['metric_name']

ex_columns = sorted(ex_columns)
print(metrics_df[ex_columns].drop_duplicates(), flush=True)

# save extra columns
with open(extra_columns, 'w') as f:
    for col in ex_columns:
        f.write(f'{col}\n')

# save metrics
del metrics_df['overwrite_file_id']
metrics_df.to_csv(out_tsv, sep='\t', index=False)

print(
    metrics_df[['metric', 'output_type', 'metric_type', 'score', 's', 'h:m:s']],
    flush=True
)
