import pandas as pd
import warnings

warnings.filterwarnings('ignore')

# load snakemake variables
input_metrics = snakemake.input.metrics
input_benchmark = snakemake.input.benchmark
out_tsv = snakemake.output.tsv
extra_columns_file = snakemake.output.extra_columns
expanded_wildcards = snakemake.params['wildcards']
max_len = int(snakemake.params.get('max_file_name_len', 100))

# efficient data loading
metrics_df = pd.concat([pd.read_table(f) for f in input_metrics], ignore_index=True)
benchmark_df = pd.concat([pd.read_table(f) for f in input_benchmark], ignore_index=True)

# merge wildcards and potential benchmarks
# Assuming expanded_wildcards aligns row-wise with input_benchmark
benchmark_df = pd.concat([expanded_wildcards, benchmark_df], axis=1)

# merge metrics and benchmarks using left join ensures not losing metrics if bench missing
metrics_df = metrics_df.merge(benchmark_df, how='left').drop_duplicates()

# vectorized metada extraction (better than for loop for bigger data bases)
mask_overwrite = metrics_df.get('overwrite_file_id', False) == True

# case a: simple file_id
metrics_df.loc[~mask_overwrite, 'file_name'] = metrics_df.loc[~mask_overwrite, 'file_id']

# case b: complex file_id parsing
if mask_overwrite.any():
    # split the file_id strings by '--' and expand into a long-format
    # 
    extracted = (
        metrics_df.loc[mask_overwrite, 'file_id']
        .str.split('--', expand=True)
        .stack()
        .reset_index(level=1, drop=True)
    )

    # key value pairs
    kv_mask = extracted.str.contains('=')
    if kv_mask.any():
        kv_split = extracted[kv_mask].str.split('=', n=1, expand=True)
        for key, group in kv_split.groupby(0):
            # Map values back to original indices
            metrics_df.loc[group.index, key] = group[1]

    # colon paths or raw strings
    simple_mask = ~kv_mask
    if simple_mask.any():
        
        names = extracted[simple_mask].apply(lambda x: x.split(':')[-1] if ':' in x else x)
        
        # combine multiple file names
        combined_names = names.groupby(level=0).agg('--'.join)
        
        # if already exists append it
        existing_names = metrics_df.loc[combined_names.index, 'file_name']
        metrics_df.loc[combined_names.index, 'file_name'] = (
            combined_names if existing_names.isna() 
            else existing_names + '--' + combined_names
        )

# final clean up and formatting
if 'file_name' in metrics_df.columns:
    metrics_df['file_name'] = metrics_df['file_name'].apply(
        lambda x: (x[:max_len] + '...') if isinstance(x, str) and len(x) > max_len else x
    )

# rename metric column
metrics_df['metric'] = metrics_df.get('metric_name', '')

# identification of extra columns
standard_cols = set(expanded_wildcards.columns.tolist() + ['score', 'metric_type', 'metric_name', 'metric', 'overwrite_file_id'])
ex_columns = sorted([c for c in metrics_df.columns if c not in standard_cols])

# save extra columns list
with open(extra_columns_file, 'w') as f:
    for col in ex_columns:
        f.write(f'{col}\n')

# final export
metrics_df = metrics_df.drop(columns=['overwrite_file_id'], errors='ignore')
metrics_df.to_csv(out_tsv, sep='\t', index=False)

print("Processing complete. Final shape:", metrics_df.shape)
