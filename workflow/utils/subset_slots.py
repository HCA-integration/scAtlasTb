from pathlib import Path
import shutil
from tqdm import tqdm
import numpy as np
import pandas as pd
from dask import array as da


## Writing subset masks
def set_mask_per_slot(slot, mask, out_dir, in_slot=None, in_dir=None):
    
    def _call_function_per_slot(func, path, *args, **kwargs):
        if path.is_dir() and (path / '.zattrs').exists():
            func(*args, **kwargs)
    
    link_slot = in_slot is not None
    
    if slot == 'uns':
        # do not set mask for uns slots
        return
    
    out_dir = Path(out_dir)
    mask_dir = out_dir / 'subset_mask'
    mask_dir = init_mask_dir(mask_dir, slot, in_slot, in_dir)
    
    if mask is not None:
        if slot.startswith('obs'):
            mask = mask[0]
        elif slot.startswith('var'):
            mask = mask[1]
    
    match slot:
        
        case _ if slot == 'X' or slot.startswith('layers/'):
            save_feature_matrix_mask(
                mask_dir=mask_dir,
                mask=mask,
                link_slot=link_slot
            )

        case 'obs' | 'var':
            save_slot_mask(
                mask_dir=mask_dir,
                mask=mask,
                link_slot=link_slot
            )
        
        case 'layers':
            for path in (out_dir / slot).iterdir():
                _call_function_per_slot(
                    func=save_feature_matrix_mask,
                    path=path,
                    mask_dir=mask_dir / path.name,
                    mask=mask,
                    link_slot=link_slot,
                )
        
        case 'obsp' | 'obsm' | 'varp' | 'varm':
            for path in (out_dir / slot).iterdir():
                _call_function_per_slot(
                    save_slot_mask,
                    path=path,
                    mask_dir=mask_dir / path.name,
                    mask=mask,
                    link_slot=link_slot,
                )
        
        case 'raw':
            if not mask_dir.exists():
                return
            for path in mask_dir.iterdir():
                # recursive call for subdirectories
                _call_function_per_slot(
                    set_mask_per_slot,
                    path=path,
                    slot=slot,
                    mask=mask,
                    out_dir=out_dir / slot,
                    in_slot=in_slot,
                )

        case _:
            print(f'Unknown slot: {slot}', flush=True)


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
    
    if in_slot is not None:
        remove_path(mask_dir_slot)
        in_dir = (in_dir / 'subset_mask' / in_slot).resolve()
        if in_dir.exists():
            # Copy the folder from the original location
            shutil.copytree(in_dir, mask_dir_slot)
    
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
    assert mask_old[mask_old].shape == mask.shape, \
        f'{mask_old[mask_old].shape} != {mask.shape}\n {mask_file}'
    
    mask_old[mask_old] &= mask
    return mask_old


def remove_mask_file(mask_file, link_slot):
    if link_slot:
        return
    remove_path(mask_file, remove_file=True)
    

# def set_mask_old(mask, out_dir):
#     obs_mask, var_mask = mask
#     mask_dir = out_dir / 'mask'
    
#     if mask_dir.exists():
    
#         # update mask
#         obs_mask_old = np.load(out_dir / 'mask' / 'obs_mask.npy')
#         var_mask_old = np.load(out_dir / 'mask' / 'var_mask.npy')
        
#         assert obs_mask_old[obs_mask_old].shape == obs_mask.shape, \
#             f'{obs_mask_old[obs_mask_old].shape} != {obs_mask.shape}'
#         assert var_mask_old[var_mask_old].shape == var_mask.shape, \
#             f'{var_mask_old[var_mask_old].shape} != {var_mask.shape}'
        
#         # Update old masks
#         obs_mask_old[obs_mask_old] &= obs_mask
#         obs_mask = obs_mask_old
        
#         var_mask_old[var_mask_old] &= var_mask
#         var_mask = var_mask_old

#         # remove pre-existing symlink
#         mask_dir.unlink()
    
#     # create new directory
#     mask_dir.mkdir()
#     np.save(mask_dir / 'obs_mask.npy', obs_mask)
#     np.save(mask_dir / 'var_mask.npy', var_mask)


## Functions for subsetting slots when reading them

# def subset_slots(slots, mask_dir, chunks, verbose=True):
#     for slot_name, slot in tqdm(slots.items(), desc='subset slots', disable=not verbose):
#         if isinstance(slot, dict) and slot != 'uns':
#             for key, value in slot.items():
#                 slot = subset_slot(
#                     slot_name=slot_name,
#                     slot=value,
#                     mask_dir=mask_dir,
#                     chunks=chunks
#                 )
#                 slots[slot_name][key] = slot
#         else:
#             slot = subset_slot(
#                 slot_name=slot_name,
#                 slot=slot,
#                 mask_dir=mask_dir,
#                 chunks=chunks
#             )
#             slots[slot_name] = slot
#     return slots


def subset_slot(slot_name, slot, mask_dir, chunks=('auto', -1)):
    if slot is None or not mask_dir.exists():
        return slot
    
    if isinstance(slot, dict):
        slot = {
            key: subset_slot(
                slot_name=f'{slot_name}/{key}',
                slot=value,
                mask_dir=mask_dir / key,
                chunks=chunks
            ) for key, value in slot.items()
        }
        return slot
    
    elif slot_name == 'raw':
        # dead code
        mask_dir = mask_dir / 'raw'
        for slot_name in ALL_SLOTS.items():
            slot.X = _subset_matrix(slot.X, mask_dir)
            slot.var = _subset_slot(slot_name, slot.var, mask_dir)
            slot.varm = _subset_slot(slot_name, slot.varm, mask_dir)
    
    elif slot_name.startswith('raw/'):
        slot_name = slot_name.split('/', 1)[-1]
        mask_dir = mask_dir.parent / 'raw' / mask_dir.name
        slot = subset_slot(
            slot_name=slot_name,
            slot=slot,
            mask_dir=mask_dir
        )
    
    elif slot_name == 'X':
        slot = _subset_matrix(slot, mask_dir / slot_name)

    elif slot_name.startswith('layers/'):
        slot = _subset_matrix(slot, mask_dir / slot_name)
    
    elif slot_name != 'uns':
        slot = _subset_slot(slot_name, slot, mask_dir / slot_name)

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
        return slot # _subset_slot_old(slot_name, slot, mask_dir, chunks)
    
    obs_mask = np.load(obs_mask_file)
    var_mask = np.load(var_mask_file)
    
    return slot[obs_mask, :][:, var_mask].copy()


def _subset_slot(slot_name, slot, mask_dir):
    mask_file = mask_dir / 'mask.npy'
    
    if not mask_file.exists():
        return slot # _subset_slot_old(slot_name, slot, mask_dir, chunks)
    
    mask = np.load(mask_file)
    
    if slot_name.startswith('obsp/') or slot_name.startswith('varp/'):
        # subset in both dimensions
        return slot[mask, :][:, mask].copy()
    
    return slot[mask].copy()


# def _subset_slot_old(slot_name, slot, mask_dir, chunks):
#     obs_mask = np.load(mask_dir / 'obs_mask.npy')
#     var_mask = np.load(mask_dir / 'var_mask.npy')
    
#     n_obs = len(obs_mask)
#     n_vars = len(var_mask)
    
#     if slot_name in ['X', 'layers', 'raw']:
#         if slot.shape[0] == n_obs:
#             slot = slot[obs_mask, :]
#         if slot.shape[1] == n_vars:
#             slot = slot[:, var_mask]
#         if isinstance(slot, da.Array):
#             return slot.copy().rechunk(chunks)
#         return slot.copy()
    
#     if slot_name.startswith('obs') and slot.shape[0] == n_obs:
#         if slot_name == 'obsp':
#             slot = slot[obs_mask, :][:, obs_mask].copy()
#         slot = slot[obs_mask].copy()
#         if slot_name == 'obs':
#             for col in slot.columns:
#                 if slot[col].dtype.name == 'category':
#                     slot[col] = slot[col].cat.remove_unused_categories()
#         return slot
    
#     if slot_name.startswith('var') and slot.shape[0] == n_vars:
#         if slot_name == 'varp':
#             return slot[var_mask, :][:, var_mask].copy()
#         return slot[var_mask].copy()
    
#     return slot

