"""
Taken and adjusted from https://github.com/openpipelines-bio/openpipeline/blob/main/src/integrate/scarches/script.py
"""

from anndata import AnnData
import scvi
import pandas as pd
import numpy as np
from pathlib import Path
import logging
logging.basicConfig(level=logging.INFO)

# Monkey-patch torch.load to use CPU when CUDA unavailable
if not torch.cuda.is_available():
    _original_torch_load = torch.load
    torch.load = lambda *args, **kwargs: _original_torch_load(
        *args, **{**kwargs, "map_location": torch.device("cpu")}
    )


def _read_model_name_from_registry(model_path) -> str:
    """Read registry with information about the model, return the model name"""
    registry = scvi.model.base.BaseModelClass.load_registry(model_path)
    return registry["model_name"]


def _detect_base_model(model_path):
    """Read from the model's file which scvi_tools model it contains"""
    return getattr(scvi.model, _read_model_name_from_registry(model_path))


def _validate_obs_metadata_params(model_registry, model_name, **kwargs):
    """
    Validates .obs metadata parameters that are required by scvi-tools models.

    This function checks that necessary .obs metadata field names (batch, labels, size factors)
    specified in the model registry have the required corresponding parameters provided in the input.

    Parameters
    ----------
    model_registry : dict
        Dictionary containing model configuration and requirements.
    model_name : str
        Name of the model being validated.
    kwargs:
        Key-value mapping for the specific model parameters (e.g. batch_key, labels_key)
    """
    for registry_key in [
        'batch_key',
        'labels_key',
        'size_factor_key',
        'categorical_covariate_keys',
        'continuous_covariate_keys',
    ]:
        model_value = model_registry.get(registry_key)
        query_value = kwargs.get(registry_key)

        if model_value and not query_value:
            if registry_key == "labels_key" and model_registry.get("unlabeled_category"):
                # skip optional labels_key
                continue
            raise ValueError(
                f"The provided {model_name} model requires `--{registry_key}` to be provided."
            )

        elif query_value and not model_value:
            logging.warning(
                f"`--{registry_key}` was provided but is not used in the provided {model_name} model."
            )

        elif model_value and query_value and registry_key in [
            "categorical_covariate",
            "continuous_covariate",
        ]:
            if len(query_value) != len(model_value):
                raise ValueError(
                    f"The number of provided covariates in `--{registry_key}` ({query_value}) does not "\
                    f"match the number of covariates used in the provided model ({model_value})."
                )


def _align_query_with_registry(adata_query, model_path, **kwargs):
    """
    Creates a query AnnData object with the expected structure and metadata fields that are aligned with the pre-trained reference model.

    Parameters
    ----------
    adata_query : AnnData
        The query AnnData object to be aligned with the model structure.
    model_path : str
        Path to the directory containing the pre-trained model.
    kwargs:
        Key-value mapping for the specific model parameters (e.g. batch_key, labels_key)

    Returns
    -------
    AnnData
        A new AnnData object with structure and metadata aligned to match the
        requirements of the pre-trained model.
    """

    model = _detect_base_model(model_path)
    model_name = _read_model_name_from_registry(model_path)
    model_registry = model.load_registry(model_path)["setup_args"]
    _validate_obs_metadata_params(model_registry, model_name, **kwargs)

    # align observations
    query_obs = {}
    
    # most common model keys
    for registry_key in [
        'batch_key', # relevant for AUTOZI, LinearSCVI, PEAKVI, SCANVI, SCVI, TOTALVI, MULTIVI, JaxSCVI
        'size_factor_key' # relevant for SCANVI, SCVI, TOTALVI, MULTIVI
    ]:
        model_key = model_registry.get(registry_key)
        query_key = kwargs.get(registry_key)
        if model_key:
            assert query_key, f"Expected {registry_key} to be provided"
            assert query_key in adata_query.obs.columns, (
                f"The provided {registry_key} '{query_key}' was not found in the query .obs "
                f"{adata_query.obs.columns.tolist()}"
            )
            query_obs[model_key] = adata_query.obs[query_key].values

    ## labels-key, relevant for AUTOZI, CondSCVI, LinearSCVI, PEAKVI, SCANVI, SCVI
    model_labels_key = model_registry.get('labels_key')
    labels_key = kwargs.get('labels_key')
    if model_labels_key:
        if not labels_key:
            labels_key = model_labels_key
            # unlabeled by default
            adata_query.obs[labels_key] = model_registry["unlabeled_category"]
        query_obs[model_labels_key] = adata_query.obs[labels_key].values

    ## covariate keys
    for registry_key in [
        'categorical_covariate_keys',
        'continuous_covariate_keys'
    ]:
        for covariate in model_registry.get(registry_key, []):
            query_obs[covariate] = adata_query.obs[covariate].values

    adata_query.obs = pd.DataFrame(query_obs, index=adata_query.obs_names)
    adata_query.var = pd.DataFrame(index=adata_query.var_names)

    if model_registry.get("layer"):
        adata_query.layers[model_registry["layer"]] = kwargs.get("query_layer", "X")

    return adata_query
