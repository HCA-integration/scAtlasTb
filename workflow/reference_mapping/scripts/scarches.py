import logging
logging.basicConfig(level=logging.INFO)
from pathlib import Path
import pandas as pd
import numpy as np
import scarches
import torch

from utils.io import read_anndata, write_zarr_linked
from utils.misc import dask_compute
from openpipelines_functions import _align_query_with_registry, _detect_base_model


input_file = snakemake.input.zarr
model_input = snakemake.input.model
output_file = snakemake.output.zarr
model_output = snakemake.output.model

tmpdir = snakemake.resources.tmpdir
layer = snakemake.params.get('layer', 'X')
var_key = snakemake.params.get('var_key')

model_params = dict(
    freeze_dropout=True,
) | snakemake.params.get('model_params', {})

train_params = dict(
    max_epochs=10,
    early_stopping=True,
    check_val_every_n_epoch=1,
) | snakemake.params.get('train_kwargs', {})

batch_key = model_params.pop('batch_key')
labels_key = model_params.pop('labels_key', None)
size_factor_key = model_params.pop('size_factor_key', None)
categorical_covariate_keys = model_params.pop('categorical_covariate_keys', [])
continuous_covariate_keys = model_params.pop('continuous_covariate_keys', [])

logging.info('Get reference model...')
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
model_torch = torch.load(model_input, map_location=device, weights_only=False)
assert 'var_names' in model_torch.keys(), f'var_names not found in model file {model_input}, make sure to provide a valid scArches model'

logging.info('Read adata...')
adata = read_anndata(
    input_file,
    X=layer,
    obs='obs',
    var='var',
)

# subset to genes used by reference model
logging.info('Subsetting genes')
var_names = adata.var_names if var_key is None else adata.var[var_key]
adata.var['scarches_features'] = var_names.isin(model_torch['var_names'])
assert adata.var['scarches_features'].sum() > 0, \
    'No overlapping genes between query and reference model.' \
    f'{var_names}\n{model["var_names"]}'

# gene order must be the same
adata = adata[:, adata.var['scarches_features']].copy()
dask_compute(adata)

logging.info('Detect base model')
model_path = str(Path(model_input).parent)
model = _detect_base_model(model_path)
logging.info(model)

logging.info('Align query to model')
adata = _align_query_with_registry(
    adata,
    model_path,
    batch_key=batch_key,
    labels_key=labels_key, # Note: called "labels_key" by SCVI
    size_factor_key=size_factor_key,
    categorical_covariate_keys=categorical_covariate_keys,
    continuous_covariate_keys=continuous_covariate_keys,
)
logging.info(adata.__str__())

# prepare SCVI model
adata.varm.clear()
model.prepare_query_anndata(
    adata,
    model_path
)
## shouldn't fail, since genes were matched before
# except ValueError:
#     model_genes = model.prepare_query_anndata(
#         adata,
#         model_path,
#         return_reference_var_names=True
#     ).tolist()
#     raise ValueError(
#         "Could not perform model.prepare_model_anndata, likely because the model was trained with"\
#         "different var names then were found in the index. \n\n"\
#         f"model var_names: {model_genes} \n\n"\
#         f"query data var_names: {adata.var_names.tolist()}"
#     )

try:
    # Load query data into the model
    vae_query = model.load_query_data(
        adata,
        model_path,
        **model_params,
    )
except KeyError:
    from tempfile import TemporaryDirectory

    logging.info(
        "Older models do not have the 'setup_method_name' key saved with the model."
        "Assume these were generated from an anndata object."
    )
    model_torch["attr_dict"]["registry_"]["setup_method_name"] = "setup_anndata"

    with TemporaryDirectory(dir=tmpdir) as tempdir:
        temp_file_name = Path(tempdir) / "model.pt"
        torch.save(Path(model_torch).parent, temp_file_name)
        del model_torch
        vae_query = model.load_query_data(
            adata,
            tempdir,
            **model_params,
        )

# Train scArches model for query mapping
vae_query.train(**train_params)

# get latent representation
adata.obsm['X_emb'] = vae_query.get_latent_representation(adata=adata)
logging.info(adata.__str__())

# save
vae_query.save(model_output, overwrite=True)

logging.info(f'Write to {output_file}...')
write_zarr_linked(
    adata=adata,
    in_dir=input_file,
    out_dir=output_file,
    files_to_keep=['obsm']
)