import logging
logging.basicConfig(level=logging.INFO)
import pandas as pd

from matplotlib import pyplot as plt
import seaborn as sns
from matplotlib import rcParams
# figure size in inches
rcParams['figure.figsize'] = 12, 9


input_file = snakemake.input.tsv
output_bar = snakemake.output.barplot
output_violin = snakemake.output.violinplot
dataset = snakemake.wildcards.dataset
file_id = snakemake.wildcards.file_id
n_permute = snakemake.params.n_permute

logging.info('Read TSV...')
df = pd.read_table(input_file)

if df.shape[0] == 0:
    logging.info('Empty TSV, skip plotting')
    plt.savefig(output_bar)
    plt.savefig(output_violin)
    exit(0)

title = 'Principal Component Regression scores of covariates for'
title += f'\ndataset={dataset}\nfile_id={file_id}\n{n_permute} permutations'

df = df.sort_values(
    ['pcr', 'n_covariates', 'covariate'],
    ascending=[False, True, False]
)

# Adjust figure size based on the number of covariates
num_covariates = df['covariate'].nunique()
fig, ax = plt.subplots(figsize=(10, 5 + num_covariates * 0.75))

logging.info('Barplot...')
g = sns.barplot(
    data=df,
    x='pcr',
    y='covariate',
    hue='permuted',
    errorbar='sd',
    dodge=True,
    err_kws={'linewidth': 1},
    capsize=.1,
    ax=ax,
)
g.set(title=title)

def round_values(x, prefix='', n_digits=3):
    if x < 10 ** (-n_digits):
        return f'{prefix}{x:.2e}'
    elif pd.notna(x):
        return f'{prefix}{x:.{n_digits}f}'
    return ''

# Prepend "pcr=" to the string
df['pcr_string'] = df['pcr'].apply(round_values, prefix='pcr=')
df['n_covariates'] = df['n_covariates'].apply(lambda x: f'n={x}')
# df['signif'] = df['z_score'].apply(lambda x: '**' if x > 3 else '*' if x > 1.5 else '')
df['signif'] = df['p-val'].apply(lambda x: '**' if x <= 0.01 else '*' if x <= 0.05 else '')
df['z_score'] = df['z_score'].apply(round_values, prefix='z=', n_digits=2)
df['p-val'] = df['p-val'].apply(round_values, prefix='p-val=', n_digits=3)

# create bar labels for covariate
covariate_bar_labels = df.groupby('covariate', sort=False).first()[
    ['pcr_string', 'n_covariates', 'z_score', 'p-val', 'signif']
].astype(str).agg(lambda x: ', '.join([s for s in x if s]), axis=1)
g.bar_label(g.containers[0], labels=covariate_bar_labels, padding=10)
logging.info(covariate_bar_labels)

# create bar labels for permuted covariates
if len(g.containers) > 1:
    perm_bar_labels = df.groupby('covariate', sort=False).first()['perm_std'].apply(round_values, prefix='std=')

    for rect, label in zip(g.containers[1], perm_bar_labels):
        x = rect.get_width()
        y_center = rect.get_y() + rect.get_height() / 2

        g.annotate(
            label,
            xy=(x, y_center),
            xytext=(10, 10),          # horizontal and vertical offset in POINTS (x, y)
            textcoords="offset points",
            va="center",
            ha="left",
            fontsize=10,
        )

plt.xticks(rotation=90)
sns.despine()

logging.info('Save barplot...')
plt.savefig(output_bar, bbox_inches='tight', dpi=300)

logging.info('Violin plot...')
plt.clf()
plt.grid()

# Create the violin for permutations
fig, ax = plt.subplots(figsize=(10, 4 + num_covariates * 0.25))

sns.violinplot(
    data=df[df['permuted']],
    x='pcr',
    y='covariate',
    fill=False,
    inner='quart',
    dodge=False,
    ax=ax,
)

# Overlay points for non-permuted values
sns.stripplot(
    data=df[~df['permuted']],
    x='pcr',
    y='covariate',
    color='red',      # pick a color that stands out
    size=8,
    marker='o',
    dodge=False,
    ax=ax,
)

ax2 = ax.twinx()
ax2.set_ylim(ax.get_ylim())

# set side labels directly
ax2.set_yticks(ax.get_yticks())
ax2.set_yticklabels([covariate_bar_labels.get(c, "") for c in df['covariate'].unique()])

# hide ticks and spine
ax2.tick_params(left=False, right=False, length=0)
ax2.spines['right'].set_visible(False)

sns.despine()
ax.set(title=title)

logging.info('Save violin plot...')
plt.savefig(output_violin, bbox_inches='tight', dpi=300)