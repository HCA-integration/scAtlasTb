import numpy as np
import pandas as pd
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

from utils.io import read_anndata, write_zarr_linked
from qc_utils import get_thresholds, apply_thresholds, parse_autoqc


def thresholds_to_df(df, wildcards, qc_type=None, **kwargs):
    if qc_type is None: 
        if 'user_thresholds' in kwargs:
            qc_type = 'user'
            if 'autoqc_thresholds' in kwargs:
                qc_type = 'updated'
        elif 'autoqc_thresholds' in kwargs:
            qc_type = 'sctk_autoqc'
        else:
            qc_type = str(qc_type)
    
    # get thresholds
    thresholds = get_thresholds(**kwargs, transform=False)
    
    # sort threshold columns
    threshold_keys = sorted(list(thresholds.keys()), reverse=True)
    thresholds = {key: thresholds[key] for key in threshold_keys}
    
    if df.shape[0] == 0:
        thresholds = dict(
            dataset=wildcards.dataset,
            file_id=wildcards.file_id,
            threshold_type=qc_type,
            passed_frac=0,
            removed_frac=0,
            n_passed=0,
            n_removed=0,
            n_total=0,
        ) | thresholds
        return pd.DataFrame(thresholds, index=[0])
    
    # apply thresholds
    df['passed_qc'] = True
    for key in kwargs['threshold_keys']:
        df['passed_qc'] = df['passed_qc'] & df[key].between(thresholds[f'{key}_min'], thresholds[f'{key}_max'])
    passed_qc = pd.Categorical(df['passed_qc'], categories=[True, False]).value_counts()
    
    thresholds = dict(
        dataset=wildcards.dataset,
        file_id=wildcards.file_id,
        threshold_type=qc_type,
        passed_frac=passed_qc[True] / passed_qc.sum(),
        removed_frac=passed_qc[False] / passed_qc.sum(),
        n_passed=passed_qc[True],
        n_removed=passed_qc[False],
        n_total=passed_qc.sum(),
    ) | thresholds
    
    return pd.DataFrame(thresholds, index=[0])


input_file = snakemake.input[0]
output_file = snakemake.output[0]
output_tsv = snakemake.output.tsv
output_qc_stats = snakemake.output.qc_stats

adata = read_anndata(input_file, obs='obs', uns='uns', verbose=False)

threshold_keys = snakemake.params.get('scautoqc_metrics')
user_thresholds = snakemake.params.get('thresholds')
alternative_thresholds = snakemake.params.get('alternative_thresholds')
autoqc_thresholds = adata.uns.get('scautoqc_ranges')
if autoqc_thresholds is None:
    autoqc_thresholds = pd.DataFrame()

# Calculate threshold stats
df = pd.concat(
    [
        # autoqc thresholds
        thresholds_to_df(
            df=adata.obs,
            wildcards=snakemake.wildcards,
            threshold_keys=threshold_keys,
            autoqc_thresholds=autoqc_thresholds,
            init_nan=True,
        ),
        # user thresholds
        thresholds_to_df(
            df=adata.obs,
            wildcards=snakemake.wildcards,
            threshold_keys=threshold_keys,
            user_thresholds=user_thresholds,
            init_nan=True,
        ), 
        # alternative thresholds
        thresholds_to_df(
            df=adata.obs,
            wildcards=snakemake.wildcards,
            qc_type='alternative',
            threshold_keys=threshold_keys,
            user_thresholds=alternative_thresholds,
            init_nan=True,
        ),
        # updated thresholds
        thresholds_to_df(
            df=adata.obs,
            wildcards=snakemake.wildcards,
            threshold_keys=threshold_keys,
            autoqc_thresholds=autoqc_thresholds,
            user_thresholds=user_thresholds,
            init_nan=True,
        ),
    ]
)
adata.uns['qc'] = df

df.to_csv(output_tsv, sep='\t', index=False)

# add QC status column to .obs with 'passed', 'failed', and 'ambiguous'
apply_thresholds(
    adata,
    thresholds=get_thresholds(
        threshold_keys,
        user_thresholds=user_thresholds,
        autoqc_thresholds=autoqc_thresholds,
    ),
    threshold_keys=threshold_keys,
    column_name='user_qc_status',
)

# set defaults for alternative thresholds
alternative_thresholds = parse_autoqc(autoqc_thresholds) | user_thresholds | alternative_thresholds
apply_thresholds(
    adata,
    thresholds=get_thresholds(
        threshold_keys,
        user_thresholds=alternative_thresholds,
        autoqc_thresholds=autoqc_thresholds,
    ),
    threshold_keys=threshold_keys,
    column_name='alternative_qc_status',
)

# get 'passed' if user_qc_status and alternative_qc_status are both True
user_status = adata.obs['user_qc_status']
alt_status = adata.obs['alternative_qc_status']
adata.obs['qc_status'] = 'ambiguous'
adata.obs.loc[user_status & alt_status, 'qc_status'] = 'passed'
adata.obs.loc[~(user_status | alt_status), 'qc_status'] = 'failed'

# convert QC status to ordered categorical
adata.obs['qc_status'] = pd.Categorical(
    adata.obs['qc_status'],
    categories=['ambiguous', 'failed', 'passed'],
    ordered=True
)

# calculate QC stats
qc_status_counts = adata.obs['qc_status'].value_counts()
qc_status_counts = pd.DataFrame(qc_status_counts).T
qc_status_counts['file_id'] = snakemake.wildcards.file_id
qc_status_counts = qc_status_counts[['file_id', 'passed', 'failed', 'ambiguous']]

qc_status_counts.to_csv(output_qc_stats, sep='\t', index=False)

write_zarr_linked(
    adata,
    input_file,
    output_file,
    files_to_keep=['obs', 'uns'],
    verbose=False,
)