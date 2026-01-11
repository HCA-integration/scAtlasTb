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

logger.info('Read TSV...')
df = pd.read_table(input_file)

if df.shape[0] == 0:
    logger.info('Empty TSV, skip plotting')
    plt.savefig(output_bar)
    plt.savefig(output_violin)
    exit(0)

title = 'Principal Component Regression scores of covariates for'
title += f'\ndataset={dataset}\nfile_id={file_id}\n{n_permute} permutations'

df = df.sort_values(
    ['pcr', 'n_covariates', 'covariate'],
    ascending=[False, True, False]
)

logger.info('Barplot...')
g = sns.barplot(
    data=df,
    x='pcr',
    y='covariate',
    hue='permuted',
    errorbar='sd',
    dodge=True,
    err_kws={'linewidth': 1},
    capsize=.1,
)
g.set(title=title)

def round_values(x, prefix='', n_digits=3):
    if x < 10 ** (-n_digits):
        f'{prefix}{x:.2e}'
    elif pd.notna(x):
        return f'{prefix}{x:.{n_digits}f}'
    return ''

df['pcr_string'] = df['pcr'].apply(round_values)
df['n_covariates'] = df['n_covariates'].apply(lambda x: f'n={x}')
df['signif'] = df['z_score'].apply(lambda x: '**' if x > 3 else '*' if x > 1.5 else '')
df['z_score'] = df['z_score'].apply(round_values, prefix='z=', n_digits=2)

# create bar labels for covariate
covariate_bar_labels = df.groupby('covariate', sort=False).first()[
    ['pcr_string', 'z_score', 'n_covariates', 'signif']
].astype(str).agg(lambda x: ', '.join([s for s in x if s]), axis=1)
print(covariate_bar_labels)
g.bar_label(g.containers[0], labels=covariate_bar_labels, padding=10)

# create bar labels for permuted covariates
if len(g.containers) > 1:
    perm_bar_labels = df.groupby('covariate', sort=False).first()['perm_std'].apply(round_values, prefix='std=')
    g.bar_label(g.containers[1], labels=perm_bar_labels, padding=25)

plt.xticks(rotation=90)
sns.despine()

logger.info('Save barplot...')
plt.savefig(output_bar, bbox_inches='tight', dpi=300)

logging.info('Violin plot...')
plt.clf()
plt.grid()
g = sns.violinplot(
    data=df,
    x='pcr',
    y='covariate',
    hue='permuted',
    dodge=False,
)
sns.despine()
ax.set(title=title)

logger.info('Save violin plot...')
plt.savefig(output_violin, bbox_inches='tight',dpi=300)