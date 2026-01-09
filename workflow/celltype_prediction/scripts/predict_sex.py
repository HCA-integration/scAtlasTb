import logging
import numpy as np
import pandas as pd

from utils.io import read_anndata, write_zarr_linked
from utils.accessors import match_genes

logging.basicConfig(level=logging.INFO)


Y_GENES = [
    "FAM197Y6",
    "FAM197Y7",
    "FAM41AY2",
    "LINC00279",
    "SRY",
    "TTTY1B",
]
X_GENES = ["XIST"]


input_file = snakemake.input[0]
output_file = snakemake.output[0]

layer = snakemake.params.get("layer")
is_normalized = snakemake.params.get("is_normalized", True)
params = snakemake.params.get("params", {})
donor_key = params.get("donor_key")
donors = params.get("donors")
x_genes = params.get("x_genes", X_GENES)
y_genes = params.get("y_genes", Y_GENES)
x_threshold = params.get("x_threshold", 0)
y_threshold = params.get("y_threshold", 4)
imbalance_frac = float(params.get("imbalance_frac", 0.1))
predict_key = params.get("predict_column", "sex")
reference_key = params.get("reference_key", predict_key)


def predict_sex_by_donor(
    adata,
    donor_key,
    predict_key,
    x_genes,
    y_genes,
    feature_column=None,
    default=None,
    x_threshold=0,
    y_threshold=4,
    frac=0.1,
):
    from fast_array_utils import stats

    x_exp = f"{predict_key}_x_exp"
    y_exp = f"{predict_key}_y_exp"

    if feature_column:
        x_genes = adata.var[feature_column].isin(x_genes)
        y_genes = adata.var[feature_column].isin(y_genes)
    else:
        x_genes = adata.var_names.isin(x_genes)
        y_genes = adata.var_names.isin(y_genes)

    # Per-cell expression (sparse-safe)
    x_exp_cell = stats.mean(adata[:, x_genes].X, axis=1)
    y_exp_cell = stats.mean(adata[:, y_genes].X, axis=1)

    # Build per-cell table
    df = pd.DataFrame(
        {
            donor_key: adata.obs[donor_key].to_numpy(),
            x_exp: x_exp_cell,
            y_exp: y_exp_cell,
        }
    )

    # Aggregate per donor
    donor_exp = df.groupby(donor_key, sort=False, observed=True).sum()

    # Vectorized classification with safe division
    x_vals = donor_exp[x_exp].to_numpy(dtype=float)
    y_vals = donor_exp[y_exp].to_numpy(dtype=float)
    y_x_frac = np.divide(
        y_vals,
        x_vals,
        out=np.full_like(y_vals, np.inf),
        where=x_vals != 0,
    )
    x_y_frac = np.divide(
        x_vals,
        y_vals,
        out=np.full_like(x_vals, np.inf),
        where=y_vals != 0,
    )
    donor_exp[predict_key] = np.select(
        [
            ((donor_exp[x_exp] > x_threshold) & (donor_exp[y_exp] <= y_threshold)) | (y_x_frac <= frac),
            ((donor_exp[x_exp] <= x_threshold) & (donor_exp[y_exp] > y_threshold)) | (x_y_frac < frac),
            ((donor_exp[x_exp] > x_threshold) & (donor_exp[y_exp] > y_threshold)) | ((y_x_frac > frac) & (x_y_frac > frac)),
        ],
        ["female", "male", "mix"],
        default=default,
    )

    return donor_exp


logging.info(f"Read file: {input_file}...")
adata = read_anndata(
    input_file,
    dask=True,
    backed=True,
    X=layer,
    obs="obs",
    var="var",
    uns="uns",
)

# parse genes
if 'feature_name' not in adata.var.columns:
    adata.var['feature_name'] = adata.var_names.astype(str)

if isinstance(x_genes, str):
    x_genes = [x_genes]

if isinstance(y_genes, str):
    y_genes = [y_genes]

x_genes_use = match_genes(adata.var, x_genes, column='feature_name', return_index=False, as_list=True)
y_genes_use = match_genes(adata.var, y_genes, column='feature_name', return_index=False, as_list=True)

assert len(x_genes_use) > 0, f'No X genes found in adata.var_names: {x_genes}'
assert len(y_genes_use) > 0, f'No Y genes found in adata.var_names: {y_genes}'

logging.info(f'Using {len(x_genes_use)} X genes: {x_genes_use}')
logging.info(f'Using {len(y_genes_use)} Y genes: {y_genes_use}')

# check donors
assert donor_key in adata.obs.columns, f'"{donor_key}" not in adata.obs.columns'
columns = [donor_key]
if reference_key in adata.obs.columns:
    columns.append(reference_key)
adata.obs = adata.obs[columns]

if not donors:
    donors = adata.obs[donor_key].unique().tolist()
else:
    adata = adata[adata.obs[donor_key].isin(donors)].copy()

# check all donors are in adata
missing_donors = set(donors) - set(adata.obs[donor_key].unique().tolist())
assert len(missing_donors) == 0, f"Donors not found in adata: {missing_donors}"

donor_exp = predict_sex_by_donor(
    adata,
    donor_key=donor_key,
    predict_key=predict_key,
    x_genes=x_genes_use,
    y_genes=y_genes_use,
    default=None,
    feature_column='feature_name',
    x_threshold=x_threshold,
    y_threshold=y_threshold,
    frac=imbalance_frac,
)

adata.obs = adata.obs.join(donor_exp[predict_key], on=donor_key)

adata.uns['predict_sex'] = {
    "x_genes": x_genes_use,
    "y_genes": y_genes_use,
    "x_threshold": x_threshold,
    "y_threshold": y_threshold,
    "imbalance_frac": imbalance_frac,
    "predictions": donor_exp,
}

logging.info(f'Predicted sex for {len(donors)} donors.')
df = adata.obs
if reference_key in df.columns:
    df = df[df[reference_key] != df[predict_key]]
    donors = df[donor_key].unique().tolist()
    donor_exp = donor_exp.query(f'{donor_key}.isin(@donors)')

logging.info(donor_exp)
logging.info(f'\n{df.value_counts(dropna=False)}')

logging.info(f'Write file: {output_file}...')
write_zarr_linked(
    adata,
    in_dir=input_file,
    out_dir=output_file,
    files_to_keep=['obs', 'uns'],
)
