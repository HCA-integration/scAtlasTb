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
    
    link_slot = in_slot is not None
    
    if slot == 'uns':
        # do not set mask for uns slots
        return
    
    out_dir = Path(out_dir)
    slot_dir = out_dir / slot
    if not slot_dir.exists():
        print(f'Slot directory {slot_dir} does not exist, skipping mask for slot {slot}', flush=True)
        return

    if mask_dir is None:
        mask_dir = out_dir / 'subset_mask'

    # set to the slot-specific mask directory
    mask_dir = init_mask_dir(mask_dir, slot, in_slot, in_dir)
    
    if mask is not None:
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
                    in_slot=in_slot,
                    in_dir=in_dir,
                )

        case _:
            print(f'Unknown slot, cannot subset: {slot}', flush=True)


def remove_path(path, remove_file=False):
    if path.is_symlink():
        path.unlink()
    elif path.is_dir():
        shutil.rmtree(path)
    elif path.is_file() and remove_file:
        path.unlink()
    

def init_mask_dir(mask_dir, slot, in_slot, in_dir):
    """
    Create if missing or copy files from symlinked directory.
    """
    if not mask_dir.exists():
        mask_dir.mkdir(parents=True)
    elif mask_dir.is_symlink():
        link_dir = mask_dir.resolve()
        mask_dir.unlink()
        # Copy the folder from the original location
        shutil.copytree(link_dir, mask_dir)
    
    mask_dir_slot = mask_dir / slot
    
    # If this slot should link (copy) an existing subset mask from another slot
    if in_slot and in_dir:
        in_dir = (in_dir / 'subset_mask' / in_slot).resolve()
        if in_dir.exists():
            # Remove existing path (dir, symlink, file) so we can create a fresh copy
            remove_path(mask_dir_slot)
            # Copy the entire mask directory from source to current slot-specific mask dir
            shutil.copytree(in_dir, mask_dir_slot)

    if not mask_dir_slot.exists():
        mask_dir_slot.mkdir(parents=True)

    return mask_dir_slot


def save_feature_matrix_mask(mask_dir, mask, link_slot):
    """
    Save the subset mask for the feature matrices e.g. from X, raw/X or layers/*
    """
    obs_mask_file = mask_dir  / 'obs.npy'
    var_mask_file = mask_dir  / 'var.npy'
    
    if mask is None:
        # remove any pre-existing mask files, since slot is being overwritten with correct shape
        # -> no subset mask required
        remove_mask_file(mask_dir, link_slot)
        return
    
    if not mask_dir.exists():
        mask_dir.mkdir(parents=True)
    
    if link_slot:
        # update subset masks
        mask = (
            update_mask(obs_mask_file, mask[0]),
            update_mask(var_mask_file, mask[1])
        )
    
    np.save(obs_mask_file, mask[0])
    np.save(var_mask_file, mask[1])


def save_slot_mask(mask_dir, mask, link_slot):
     
    if mask is None:
        # # remove any pre-existing mask file, since slot is being overwritten with correct shape
        # # -> no subset mask required
        remove_mask_file(mask_dir, link_slot)
        return

    mask_file = mask_dir / 'mask.npy'
    mask_dir.mkdir(parents=True, exist_ok=True)
    if mask_file.exists() and link_slot:
        # update mask for linked mask_file slot
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


def remove_mask_file(mask_file, link_slot):
    if link_slot:
        return
    remove_path(mask_file, remove_file=True)
    

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


def _subset_matrix(slot, mask_dir):
    obs_mask_file = mask_dir / 'obs.npy'
    var_mask_file = mask_dir / 'var.npy'
    if not obs_mask_file.exists() or not var_mask_file.exists():
        return slot
    
    obs_mask = np.load(obs_mask_file)
    var_mask = np.load(var_mask_file)
    
    return slot[obs_mask, :][:, var_mask].copy()


def _subset_slot(slot_name, slot, mask_dir):
    mask_file = mask_dir / 'mask.npy'
    
    if not mask_file.exists():
        return slot
    
    mask = np.load(mask_file)
    assert mask.shape[0] == slot.shape[0], \
        f'Mask shape {mask.shape} does not match slot shape {slot.shape} for slot {slot_name}'

    if slot_name.startswith('obsp/') or slot_name.startswith('varp/'):
        # subset in both dimensions
        return slot[mask, :][:, mask].copy()
    
    return slot[mask].copy()
