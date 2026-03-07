import warnings
warnings.filterwarnings("ignore", message="Can't initialize NVML")

import sys
import logging
logging.basicConfig(level=logging.INFO)
from pathlib import Path
import pandas as pd
import numpy as np
import anndata as ad
import torch
import scvi
from scipy import sparse

from utils.io import read_anndata, write_zarr_linked
from utils.misc import dask_compute
from openpipelines_functions import _align_query_with_registry, _detect_base_model

# Import integration utils for plotting
sys.path.append(str(Path(__file__).parent.parent.parent / 'integration' / 'scripts' / 'methods'))
from integration_utils import plot_model_history, standardize_training_logs

torch.set_float32_matmul_precision('medium')
scvi.settings.dl_num_workers = snakemake.threads

input_file = snakemake.input.zarr
model_path = snakemake.input.model
output_file = snakemake.output.zarr
model_output = snakemake.output.model
output_plot_dir = snakemake.output.plots
Path(output_plot_dir).mkdir(parents=True, exist_ok=True)

tmpdir = snakemake.resources.tmpdir
layer = snakemake.params.get('layer', 'X')
var_key = snakemake.params.get('var_key')

model_params = snakemake.params.get('model_params', {})
train_params = snakemake.params.get('train_params', {})

# Handle train_params: False means inference only, dict means training
if train_params is False:
    inference_only = True
    train_params = None
else:
    inference_only = False
    train_params = dict(
        max_epochs=10,
        plan_kwargs=dict(weight_decay=0.0), # default behaviour, freeze all weights
        check_val_every_n_epoch=1, # needed for loss curves
    ) | train_params
logging.info(f'Inference only: {inference_only}')


batch_key = model_params.pop('batch_key')
labels_key = model_params.pop('labels_key', None)
size_factor_key = model_params.pop('size_factor_key', None)
categorical_covariate_keys = model_params.pop('categorical_covariate_keys', [])
continuous_covariate_keys = model_params.pop('continuous_covariate_keys', [])

use_gpu = torch.cuda.is_available()
logging.info(f'Use gpu: {use_gpu}')

logging.info('Get reference model...')
device = torch.device("cuda" if use_gpu else "cpu")
model_torch = torch.load(Path(model_path) / 'model.pt', map_location=device, weights_only=False)
assert 'var_names' in model_torch.keys(), f'var_names not found in model "{model_path}", make sure to provide a valid scArches model'

logging.info('Read adata...')
adata = read_anndata(
    input_file,
    X=layer,
    obs='obs',
    var='var',
    dask=True,
    backed=True,
)
logging.info(adata.__str__())

# subset to genes used by reference model
logging.info('Subsetting genes')
model_genes = pd.Index(model_torch['var_names'])
var_names = adata.var_names if var_key is None else adata.var[var_key]

adata = adata[:, var_names.isin(model_genes)].copy()
assert adata.n_vars > 0, 'No overlapping genes.'
logging.info(f'Found {adata.n_vars} overlapping genes.')
adata = dask_compute(adata, layers='X')

missing = model_genes.difference(var_names)
if len(missing) > 0:
    zero_padded = ad.AnnData(
        X=sparse.csr_matrix((adata.n_obs, len(missing)), dtype=adata.X.dtype),
        obs=adata.obs,
        var=pd.DataFrame(index=missing)
    )
    zero_padded = ad.concat([adata, zero_padded], axis=1)
    adata = ad.AnnData(
        X=zero_padded.X,
        obs=adata.obs,
        var=zero_padded.var
    )
    del zero_padded
    logging.info(adata.__str__())

logging.info('Detect base model')
model = _detect_base_model(model_path)
logging.info(f'Model class: {model.__name__}')

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
try:
    model.prepare_query_anndata(adata, model_path)
except ValueError:
    model_genes = model.prepare_query_anndata(
        adata,
        model_path,
        return_reference_var_names=True
    ).tolist()
    raise ValueError(
        "Could not perform model.prepare_query_anndata, likely because the model was trained with"
        " different var names then were found in the index. \n\n"
        f"model var_names: {model_genes} \n\n"
        f"query data var_names: {adata.var_names.tolist()}"
    )

# Inference only or training mode
if inference_only:
    logging.info('Loading reference model for inference only')
    vae_query = model.load(model_path, adata=adata)
else:
    logging.info('Loading query data for scArches training (standard approach with frozen reference)')
    
    # When training is forced, treat query batches as different from reference
    # This triggers scArches to extend embeddings for new batch categories
    if batch_key and batch_key in adata.obs.columns:
        query_batches = adata.obs[batch_key].unique()
        logging.info(f'Renaming query batches to distinguish from reference: {list(query_batches)}')
        # Prepend 'query_' to batch names to make them unique
        adata.obs[batch_key] = 'query_' + adata.obs[batch_key].astype(str)
        logging.info(f'Updated batch column: {adata.obs[batch_key].unique()}')
    
    try:
        # Standard scArches: freeze reference, train query extensions for new batches
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
            torch.save(model_torch, temp_file_name)
            del model_torch
            vae_query = model.load_query_data(
                adata,
                tempdir,
                **model_params,
            )
    
    # Train scArches model for query mapping
    logging.info(f'Training with parameters: {train_params}')
    vae_query.train(**train_params)
    
    # Plot training diagnostics
    plot_model_history(
        model=vae_query,
        output_dir=output_plot_dir,
        model_name=f"scArches on {model.__name__}",
    )
    logging.info(f'Saved training plots to {output_plot_dir}')

# get latent representation
X_emb = vae_query.get_latent_representation(adata=adata)

# get original adata to preserve layers etc.
if input_file.endswith('.h5ad'):
    obs_names = adata.obs_names.copy()
    adata = read_anndata(input_file, dask=True, backed=True)
    assert adata.obs_names.equals(obs_names), 'Cell names do not match after re-loading original data.'

adata.obsm['X_emb'] = X_emb
logging.info(adata.__str__())

# save
vae_query.save(model_output, overwrite=True)

logging.info(f'Write to {output_file}...')
logging.info(adata.__str__())
write_zarr_linked(
    adata=adata,
    in_dir=input_file,
    out_dir=output_file,
    files_to_keep=['obsm']
)