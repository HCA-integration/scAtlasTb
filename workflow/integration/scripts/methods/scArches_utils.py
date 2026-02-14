from pytorch_lightning.callbacks import Callback
from collections import defaultdict



SCPOLI_MODEL_PARAMS = [
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
]

SCPOLI_EARLY_STOPPING = {
    "early_stopping_metric": "val_prototype_loss",
    "mode": "min",
    "threshold": 0,
    "patience": 20,
    "reduce_lr": True,
    "lr_patience": 13,
    "lr_factor": 0.1,
}

SCGEN_MODEL_PARAMS = [
    'hidden_layer_sizes',
    'z_dimension',
    'dr_rate',
]


class ScArchesLossLogger(Callback):
    """
    Generalized Lightning callback for scArches models to capture
    all logged metrics per epoch (training and validation).
    """
    def __init__(self):
        super().__init__()
        self.history = defaultdict(lambda: {"train": [], "validation": []})

    def on_train_epoch_end(self, trainer, pl_module):
        print(f"Epoch {trainer.current_epoch} ended. Logging metrics...", flush=True)
        metrics = trainer.callback_metrics

        for key, value in metrics.items():
            if value is None:
                continue

            val = value.detach().cpu().item() if hasattr(value, "detach") else float(value)

            # Training loss keys
            if key.endswith("loss") and not key.startswith("val_"):
                metric_name = key.replace("_loss", "")
                self.history[metric_name]["train"].append(val)

            # Validation loss keys
            elif key.startswith("val_"):
                metric_name = key.replace("val_", "").replace("_loss", "")
                self.history[metric_name]["validation"].append(val)

