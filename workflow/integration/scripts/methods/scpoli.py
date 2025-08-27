import torch
from scarches.models.scpoli import scPoli
from pprint import pformat
from pathlib import Path

from integration_utils import add_metadata, get_hyperparams, remove_slots, \
    set_model_history_dtypes, plot_model_history, clean_categorical_column, SCPOLI_EARLY_STOPPING
from utils.io import read_anndata, write_zarr_linked
from utils.accessors import subset_hvg

input_file = snakemake.input[0]
output_file = snakemake.output[0]
output_model = snakemake.output.model
output_plot_dir = snakemake.output.plots
Path(output_plot_dir).mkdir(parents=True, exist_ok=True)

wildcards = snakemake.wildcards
params = snakemake.params
batch_key = wildcards.batch
label_key = wildcards.label
if label_key == 'None':
    label_key = None

torch.manual_seed(params.get('seed', 0))
torch.set_num_threads(snakemake.threads)

model_params, train_params = get_hyperparams(
    hyperparams=params.get('hyperparams', {}),
    model_params=[
        'share_metadata',
        'obs_metadata',
        'condition_keys',
        'conditions',
        'conditions_combined',
        'inject_condition',
        'cell_type_keys',
        'cell_types',
        'unknown_ct_names',
        'labeled_indices',
        'prototypes_labeled',
        'prototypes_unlabeled',
        'hidden_layer_sizes',
        'latent_dim',
        'embedding_dims',
        'embedding_max_norm',
        'dr_rate',
        'use_mmd',
        'mmd_on',
        'mmd_boundary',
        'recon_loss',
        'beta',
        'use_bn',
        'use_ln',
    ],
)

# set default model parameters
model_params = dict(
    condition_keys=[batch_key],
    cell_type_keys=[label_key],
    unknown_ct_names=['NA'],
) | model_params

# set default training parameters
n_epochs = train_params.get('n_epochs', 100)
if n_epochs is None:
    n_epochs = 100
train_params = train_params | {
    'early_stopping_kwargs': SCPOLI_EARLY_STOPPING,
    # set default pretrain epochs if not configured
    'pretrain_epochs': int(0.8 * n_epochs)
}

print(
    f'model parameters:\n{pformat(model_params)}\n'
    f'training parameters:\n{pformat(train_params)}',
    flush=True
)


# check GPU
print(f'GPU available: {torch.cuda.is_available()}', flush=True)

print(f'Read {input_file}...', flush=True)
adata = read_anndata(
    input_file,
    X='layers/counts',
    obs='obs',
    var='var',
    uns='uns',
    dask=True,
    backed=True,
)

clean_categorical_column(adata, batch_key)

cell_type_keys = model_params.get('cell_type_keys', [])
if isinstance(cell_type_keys, str):
    cell_type_keys = [cell_type_keys]
elif cell_type_keys is None:
    cell_type_keys = []
for key in cell_type_keys:
    if key is None:
        continue
    clean_categorical_column(adata, key)

# subset features
adata, _ = subset_hvg(adata, var_column='integration_features')

# prepare data for model
adata.X = adata.X.astype('float32')

print(f'Set up scPoli with parameters:\n{pformat(model_params)}', flush=True)
model = scPoli(adata=adata, **model_params)

print(f'Train scPoli with parameters:\n{pformat(train_params)}', flush=True)
model.train(**train_params)

print('Save model...', flush=True)
model.save(output_model, overwrite=True)

# prepare output adata
adata.obsm["X_emb"] = model.get_latent(adata, mean=True)
adata.uns[f"scpoli_{batch_key}_embeddings"] = model.get_conditional_embeddings().X
adata = remove_slots(adata=adata, output_type=params['output_type'], keep_X=True)
add_metadata(
    adata,
    wildcards,
    params,
    model_history=dict(model.trainer.logs)
)

plot_model_history(
    title='loss',
    train=model.trainer.logs['epoch_loss'],
    validation=model.trainer.logs['val_loss'],
    output_path=f'{output_plot_dir}/loss.png'
)

print(adata.__str__(), flush=True)
print(f'Write {output_file}...', flush=True)
write_zarr_linked(
    adata,
    input_file,
    output_file,
    files_to_keep=['obsm', 'uns']
)
