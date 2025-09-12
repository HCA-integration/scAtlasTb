import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
import warnings

warnings.filterwarnings('ignore')

input_file = snakemake.input.tsv
output_file_time = snakemake.output.time
output_file_score = snakemake.output.score

metrics_df = pd.read_table(input_file)

metrics_df['metric'] = metrics_df['metric'].astype('category')

theislab_metrics = [
    'nmi',
    'ari',
    'ilisi',
    'clisi',
    'asw_label',
    'asw_batch',
    'isolated_label_asw',
    'graph_connectivity',
    'pcr'
]
yoseflab_metrics = [
    'nmi_leiden_y',
    'ari_leiden_y',
    'ilisi_y',
    'clisi_y',
    'asw_label_y',
    'asw_batch_y',
    'bras_batch',
    'isolated_label_asw_y',
    'graph_connectivity_y',
    'pcr_y'
]

metric_name_map = dict(zip(yoseflab_metrics, theislab_metrics))

metrics_df.loc[metrics_df['metric'].isin(theislab_metrics), 'implementation'] = 'scib'
metrics_df.loc[metrics_df['metric'].isin(yoseflab_metrics), 'implementation'] = 'scib-metrics'
metrics_df = metrics_df[metrics_df['implementation'].notna()]
metrics_df['metric'] = [metric_name_map[x] if x in metric_name_map.keys() else x for x in metrics_df['metric']]


# Compare run time

df = metrics_df.pivot(
    index=['dataset', 'batch', 'label', 'hyperparams', 'metric', 'method', 'output_type'],
    columns='implementation',
    values='s'
).reset_index()
df['speedup'] = df['scib-metrics'] / df['scib']

print('Speedup of scib-metrics over scib')
print(df)

g = sns.catplot(
    data=df,
    x='speedup',
    y='metric',
    margin_titles=True,
    kind='bar',
    order=df.groupby('metric')['speedup'].max().sort_values(ascending=False).index,
    errwidth=0.5,
    capsize=0.2,
    dodge=False,
).set(
    title='Speedup of scib-metrics over scib',
    xscale='log'
)
g.set_axis_labels('t scib-metrics / scib', 'Metric')
g.refline(x=1, color='black')

g.savefig(output_file_time, bbox_inches='tight', dpi=200)


# Compare scores

df_score = metrics_df.pivot(
    index=['dataset', 'batch', 'label', 'hyperparams', 'metric', 'method', 'output_type'],
    columns='implementation',
    values='score'
).reset_index()

print('Score comparison')
print(df_score)

g = sns.relplot(
    data=df_score,
    x='scib',
    y='scib-metrics',
    hue='metric',
    s=30
)
g.ax.axline(xy1=(0, 0), slope=1, color='grey', lw=0.8, alpha=0.5)
g.set(title='Comparison of scib-metrics and scib values')

g.savefig(output_file_score, bbox_inches='tight', dpi=200)