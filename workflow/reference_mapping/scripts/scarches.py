import logging
logging.basicConfig(level=logging.INFO)
from pathlib import Path
import pandas as pd
import numpy as np
import scarches

from utils.io import read_anndata, write_zarr_linked


input_file = snakemake.input.zarr
model_file = snakemake.input.model
output_file = snakemake.output.zarr

logging.info('Read adata...')
adata = read_anndata(input_file, obs='obs', obsm='obsm')

def check_order(x, y):
    return((x.var.index == y.gene_symbol).all() or (x.var.index == y["Unnamed: 0"]).all())

surgery_epochs = snakemake.params['surgery_epochs']
early_stopping_kwargs_surgery = snakemake.params['early_stopping_kwargs_surgery']

reference_model_dir = str(Path(snakemake.input.model).parent)
# model_dir = str(Path(snakemake.output.model).parent)
batch_variable = snakemake.params.batch

adata_query = sc.read(snakemake.input.h5ad)
adata_query.obs["dataset"] = snakemake.wildcards.dataset
# adata_query.X = adata_query.layers['soupX_counts']

# remove bc incompatible with adapter training
del adata_query.varm
del adata_query.obsm
del adata_query.obsp
del adata_query.layers
del adata_query.uns

# map gene names to gene IDs
gene_map = pd.read_csv(snakemake.input.gene_map)
adata_query = utils.subset_and_pad_adata_object(adata_query, gene_map)
print(adata_query)


gene_set = gene_map.copy()
gene_set = gene_set.rename(columns={gene_set.columns[0]: "Unnamed: 0"})
# adata_query = pp.subset_and_pad_adata(gene_set, adata_query)
if check_order(adata_query, gene_set):
    print("Gene order is correct.")
else:
    print("WARNING: your gene order does not match the order of the HLCA reference.")
    adata_query = adata_query[:, gene_set.gene_symbol.tolist()].copy()
    print("Matching now?",  check_order(adata_query, gene_map))

adata_query.obs["scanvi_label"] = "unlabeled"


# split by batch
query_batches = adata_query.obs[batch_variable].unique()
ad_trained_list = []

for batch in query_batches:
    print(f"Query ({batch_variable}):", batch)
    ad = adata_query[adata_query.obs[batch_variable] == batch].copy()
    
    print("Shape:", ad.shape)
    
    # load model and set relevant variables:
    print('Load model...')
    model = scarches.models.SCANVI.load_query_data(
        ad,
        reference_model_dir,
        freeze_dropout=True,
    )
    model._unlabeled_indices = np.arange(adata_query.n_obs)
    model._labeled_indices = []
    print("Labelled Indices: ", len(model._labeled_indices))
    print("Unlabelled Indices: ", len(model._unlabeled_indices))

    # now train surgery model using reference model and query adata
    print('Train model...')
    model.train(
        max_epochs=surgery_epochs,
        early_stopping=True
    )

    # get latent representation
    ad.obsm['X_scarches'] = model.get_latent_representation(adata=ad)

    # add anndata to list
    ad_trained_list.append(ad)

adata_query = sc.concat(ad_trained_list)
logging.info(adata_query.__str__())

# save
# model.save(model_dir, overwrite=True)

logging.info(f'Write to {output_file}...')
write_zarr_linked(
    adata=adata,
    in_dir=input_file,
    out_dir=output_file,
    files_to_keep=['obsm']
)