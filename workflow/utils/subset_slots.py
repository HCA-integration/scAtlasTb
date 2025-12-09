from pathlib import Path
import shutil
from tqdm import tqdm
import numpy as np
import pandas as pd
from dask import array as da


## Writing subset masks
def set_mask_per_slot(slot, mask, out_dir, mask_dir=None, in_slot=None, in_dir=None):
    
    def _call_function_per_slot(func, path, *args, **kwargs):
        if path.is_dir() and (path / '.zattrs').exists():
            func(*args, **kwargs)

    if slot == 'uns':
        # do not set mask for uns slots
        return
    
    out_dir = Path(out_dir)
    slot_dir = out_dir / slot
    if not slot_dir.exists():
        print(f'Slot directory {slot_dir} does not exist, skipping mask for slot {slot}', flush=True)
        return

    link_slot = slot_dir.is_symlink() and in_slot is not None

    if mask_dir is None:
        mask_dir = out_dir / 'subset_mask'

    # set to the slot-specific mask directory
    mask_dir = init_mask_dir(mask_dir, slot, in_slot, in_dir)
    
    if mask is None:
        if not link_slot:
            # remove any old mask, since a new slot will be written in the correct shape
            remove_path(mask_dir)
        # if linked, keep copy of the previous mask provided by init_mask_dir
        return

    # Pick the appropriate dimension-specific mask for obs/var slots
    if slot.startswith(('obs', 'raw/obs')):
        mask = mask[0]
    elif slot.startswith(('var', 'raw/var')):
        mask = mask[1]
    
    match slot:
        
        case _ if slot in ['X', 'raw/X'] or slot.startswith(('layers/', 'raw/layers/')):
            save_feature_matrix_mask(
                mask_dir=mask_dir,
                mask=mask,
                link_slot=link_slot
            )

        case 'obs' | 'var' | 'raw/obs' | 'raw/var':
            save_slot_mask(
                mask_dir=mask_dir,
                mask=mask,
                link_slot=link_slot
            )
        
        case 'layers' | 'raw/layers':
            for path in slot_dir.iterdir():
                _call_function_per_slot(
                    func=save_feature_matrix_mask,
                    path=path,
                    mask_dir=mask_dir / path.name,
                    mask=mask,
                    link_slot=link_slot,
                )
        
        case _ if slot in {
            'obsp', 'obsm', 'varp', 'varm',
            'raw/obsp', 'raw/obsm', 'raw/varp', 'raw/varm'
        }:
            for path in slot_dir.iterdir():
                _call_function_per_slot(
                    save_slot_mask,
                    path=path,
                    mask_dir=mask_dir / path.name,
                    mask=mask,
                    link_slot=link_slot,
                )
        
        case 'raw':
            for path in slot_dir.iterdir():
                if path.name == 'raw':
                    # skip nested raw directory
                    continue
                # recursive call for subdirectories
                _call_function_per_slot(
                    set_mask_per_slot,
                    path=path,
                    slot=f'{slot}/{path.name}',
                    out_dir=out_dir,
                    mask=mask,
                    mask_dir=mask_dir.parent,
                    in_slot=f'{in_slot}/{path.name}' if in_slot else None,
                    in_dir=in_dir,
                )

        case _:
            print(f'Unknown slot, cannot subset: {slot}', flush=True)


def remove_path(path):
    if path.is_symlink():
        path.unlink()
    elif path.is_dir():
        shutil.rmtree(path)
    elif path.is_file():
        path.unlink()
    

def init_mask_dir(mask_dir, slot, in_slot, in_dir):
    """
    Create if missing or copy files from symlinked directory.
    """
    if not mask_dir.exists():
        mask_dir.mkdir(parents=True)
    elif mask_dir.is_symlink():
        # need to make a copy of the linked directory to avoid editing the original subset masks
        link_dir = mask_dir.resolve()
        mask_dir.unlink()
        # Copy the folder from the original location
        shutil.copytree(link_dir, mask_dir)
    
    mask_dir_slot = mask_dir / slot
    # remove any previous copy of the slot mask
    remove_path(mask_dir_slot)
    
    # If this slot should link (copy) an existing subset mask from another slot
    if in_slot and in_dir:
        in_dir = (in_dir / 'subset_mask' / in_slot).resolve()
        if in_dir.exists():
            # copy the entire mask directory from source to current slot-specific mask dir
            shutil.copytree(in_dir, mask_dir_slot)

    if not mask_dir_slot.exists():
        mask_dir_slot.mkdir(parents=True)

    return mask_dir_slot


def save_feature_matrix_mask(mask_dir, mask, link_slot):
    """
    Save the subset mask for the feature matrices e.g. from X, raw/X or layers/*
    """
    mask_dir.mkdir(parents=True, exist_ok=True)
    obs_mask_file = mask_dir / 'obs.npy'
    var_mask_file = mask_dir / 'var.npy'
        
    if link_slot:
        mask = (
            update_mask(obs_mask_file, mask[0]),
            update_mask(var_mask_file, mask[1])
        )
    
    np.save(obs_mask_file, mask[0])
    np.save(var_mask_file, mask[1])


def save_slot_mask(mask_dir, mask, link_slot):
    mask_dir.mkdir(parents=True, exist_ok=True)
    mask_file = mask_dir / 'mask.npy'
    
    if link_slot:
        mask = update_mask(mask_file, mask)
    
    np.save(mask_file, mask)


def update_mask(mask_file, mask):
    """
    Update the existing mask with the new mask, assuming the new mask is a subset of the old mask.
    """
    if not mask_file.exists():
        return mask
    
    mask_old = np.load(mask_file)

    if mask_old.shape == mask.shape:
        # no need to update mask, since it has the same shape
        return mask

    # check if old mask can be subsetted with the new mask
    assert mask_old[mask_old].shape == mask.shape, \
        f'Shape mismatch\n mask_old: {mask_old.shape}, mask_old[mask_old]: '\
        f'{mask_old[mask_old].shape}, mask: {mask.shape}\n{mask_file}'
    
    mask_old[mask_old] &= mask
    return mask_old
    

## Functions for subsetting slots when reading them

def subset_slot(slot_name, slot, mask_dir, chunks=('auto', -1)):
    if slot is None or not mask_dir.exists():
        return slot

    if slot_name in ('uns', 'raw/uns'):
        return slot

    if isinstance(slot, dict):
        slot = {
            key: subset_slot(
                slot_name=f'{slot_name}/{key}',
                slot=value,
                mask_dir=mask_dir,
                chunks=chunks,
            ) for key, value in slot.items()
        }
        return slot
    
    elif slot_name in ['X', 'raw/X'] or slot_name.startswith(('layers/', 'raw/layers/')):
        slot = _subset_matrix(slot, mask_dir / slot_name)
    
    elif slot_name.startswith(('obs', 'var', 'raw/obs', 'raw/var')):
        slot = _subset_slot(slot_name, slot, mask_dir / slot_name)

    # optimise data after subsetting
    if isinstance(slot, da.Array):
        slot = slot.rechunk(chunks=chunks)
    
    if isinstance(slot, pd.DataFrame):
        for col in slot.columns:
            if slot[col].dtype.name != 'category':
                continue
            slot[col] = slot[col].cat.remove_unused_categories()

    return slot


def _ensure_mask_shape(mask, expected, mask_file=None):
    if mask.shape[0] != expected:
        # Subset the mask by itself to get only the True values.
        if not np.issubdtype(mask.dtype, np.bool_):
            raise TypeError(
                f"Mask must be a boolean array to subset by itself, got dtype {mask.dtype} for mask file {mask_file}"
            )
        mask = mask[mask]
    assert mask.shape[0] == expected, \
        f'Mask shape {mask.shape} does not match expected {expected} for mask file {mask_file}'
    return mask


def _subset_matrix(slot, mask_dir):
    obs_mask_file = mask_dir / 'obs.npy'
    var_mask_file = mask_dir / 'var.npy'
    
    if not obs_mask_file.exists() or not var_mask_file.exists():
        return slot
    
    obs_mask = np.load(obs_mask_file)
    var_mask = np.load(var_mask_file)

    obs_mask = _ensure_mask_shape(obs_mask, slot.shape[0], mask_file=obs_mask_file)
    var_mask = _ensure_mask_shape(var_mask, slot.shape[1], mask_file=var_mask_file)
    
    try:
        return slot[obs_mask, :][:, var_mask].copy()
    except Exception as e:
        raise ValueError(
            f"Error subsetting matrix slot with masks from {mask_dir}\n"
            f"slot shape: {slot.shape}\n"
            f"obs_mask: shape={obs_mask.shape}, sum={obs_mask.sum()}\n"
            f"var_mask: shape={var_mask.shape}, sum={var_mask.sum()}"
        ) from e


def _subset_slot(slot_name, slot, mask_dir):
    mask_file = mask_dir / 'mask.npy'
    
    if not mask_file.exists():
        return slot
    
    mask = np.load(mask_file)
    mask = _ensure_mask_shape(mask, slot.shape[0], mask_file=mask_file)

    if slot_name.startswith('obsp/') or slot_name.startswith('varp/'):
        # subset in both dimensions
        return slot[mask, :][:, mask].copy()
    
    return slot[mask].copy()
