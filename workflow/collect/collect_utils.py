import pandas as pd
import numpy as np
from tqdm import tqdm
from pprint import pformat

LARGE_SLOTS = ['layers', 'obsm', 'obsp']


def check_obs_same_index(adatas):
    iterator = iter(adatas.items())
    key1, _ad1 = next(iterator)
    
    for key, _ad in iterator:
        if key == key1:
            continue
        diff = _ad.obs_names.difference(_ad1.obs_names)
        assert len(diff) == 0, \
            f'Index must be the same\n {len(diff)} differing indices \
            when comparing {key} against {key}\n' 


def check_columns_equal(adatas, col):
    def _check_columns_equal(s1, s2, col):
        assert s1.index.equals(s2.index), \
            f'Indices do not match for column comparison: {col}'
        # Fast path: check if both Series are the same object or have the same underlying data
        if s1 is s2:
            return True
        if s1._values is s2._values:
            return True
        # Use numpy for fast comparison if dtypes are compatible
        if s1.dtype == s2.dtype and pd.api.types.is_numeric_dtype(s1.dtype):
            return np.array_equal(s1.values, s2.values, equal_nan=True)
        # Fallback to pandas equals
        return s1.equals(s2)

    iterator = iter(adatas.values())
    _ad1 = next(iterator)
    if all(_check_columns_equal(_ad.obs[col], _ad1.obs[col], col) for _ad in iterator):
        return []
    return [col]


def get_same_columns(adatas, n_threads=1):
    check_obs_same_index(adatas)

    # collect obs columns that exist in all files
    iterator = iter(adatas.values())
    _ad1 = next(iterator)
    same_obs_columns = set(
        _ad1.obs.columns.tolist()
    ).intersection(
        *[_ad.obs.columns.tolist() for _ad in iterator]
    )
    
    # determine which columns do not have the same values across all files
    obs_to_remove = []
    for col in tqdm(same_obs_columns, desc='Check for equal columns'):
        obs_to_remove.extend(check_columns_equal(adatas, col))
    
    # from concurrent.futures import ThreadPoolExecutor, as_completed
    # with ThreadPoolExecutor(max_workers=n_threads) as executor, \
    #     tqdm(total=len(same_obs_columns), desc='Check for equal columns') as pbar:
    #     futures = [
    #         executor.submit(check_columns_equal, adatas=adatas, col=col)
    #         for col in same_obs_columns
    #     ]
    #     for future in as_completed(futures):
    #         obs_to_remove.extend(future.result())
    #         pbar.update(1)

    print(f'Shared columns that are not the same across datasets:\n{pformat(obs_to_remove)}', flush=True)
    return [col for col in same_obs_columns if col not in obs_to_remove]


def merge_df(
    df_current,
    file_id,
    df_previous,
    same_columns,
    sep='_',
):
    if df_previous is None:
        return df_current.rename(
            columns={
                col: f'{col}{sep}{file_id}'
                for col in df_current.columns.tolist()
                if col not in same_columns
            }
        )
    
    idx_diff = df_current.index.difference(df_previous.index).sort_values()
    n_diff = len(idx_diff)
    assert n_diff == 0, \
        f'Index must be the same\n {n_diff} differing indices: {idx_diff}\n' \
        f'current index: {df_current.index.sort_values()} \nprevious index: {df_previous.index.sort_values()}'
    unique_columns = [col for col in df_current.columns if col not in same_columns]
    df_current = df_current[unique_columns].rename(
        columns={col: f'{col}{sep}{file_id}' for col in unique_columns}
    )
    df = pd.concat([df_previous, df_current], axis=1)
    return df


def align_file_indices(adatas, files, merge_slots):
    """
    Align per-file AnnData object indices (obs_names and var_names) to match a reference AnnData,
    reloading large data containers only when a reordering is actually required.
    This function mutates the provided `adatas` mapping in-place:
    - Any AnnData whose obs_names or var_names differ only in order from the reference (`first_adata`)
        is reordered to match the reference.
    - Large data slots (e.g. layers / X-like backed attributes listed in the global LARGE_SLOTS)
        are reloaded from disk before reordering to avoid creating expensive view-based chained operations.
    - Views produced by subsetting are converted to real copies to ensure subsequent safe mutation.
    If an AnnData has a true mismatch in obs_names or var_names content (not just order),
    a ValueError is raised.
    Parameters
    ----------
    adatas : dict[str, anndata.AnnData]
            Mapping from file_id to AnnData. The first item is treated as the reference.
    files : dict[str, str]
            Mapping from file_id to file path on disk. Used to reload large slots via read_anndata.
    merge_slots : Iterable[str]
            Collection of slot names requested for merging (e.g. {'X', 'layers', 'obsm'}); only those
            intersecting LARGE_SLOTS are considered "large" for conditional reload.
    Returns
    -------
    dict[str, list[str]]
            A dictionary mapping each large slot name (subset of LARGE_SLOTS present in merge_slots)
            to a list of file_ids for which that slot's data was reloaded and must later be rewritten
            (e.g. when consolidating results). Empty lists indicate no reordering needed for that slot.
    Raises
    ------
    ValueError
            If any non-reference AnnData has obs_names or var_names whose set of values differs
            from the reference (i.e. not just a permutation).
    Side Effects
    ------------
    - Modifies entries of `adatas` in place (reordered indices, materialized copies).
    - Reloads large slots from disk using `read_anndata` with dask=True, backed=True to minimize memory.
    - Logs informational messages when reordering occurs.
    Notes
    -----
    - The global LARGE_SLOTS must be defined elsewhere; typical values might include 'layers', 'X', etc.
    - Only large slots actually requested in `merge_slots` are considered for reload.
    - Reordering is performed via slicing: `_ad[first_obs_names, :]` and `[:, first_var_names]`.
    - Ensuring slots are reloaded before reindexing avoids creating views that could lead to subtle bugs
        when later writing or concatenating data.
    """
    import logging
    from utils.io import read_anndata
    
    adata_iter = iter(adatas.items())
    first_key, first_adata = next(adata_iter)
    
    first_obs_names = first_adata.obs_names
    first_var_names = first_adata.var_names
    large_merge_slots = [slot for slot in LARGE_SLOTS if slot in merge_slots]
    file_copy_map = {slot: [] for slot in large_merge_slots}
    
    # Use iterator to skip first adata
    
    for file_id, _ad in adata_iter:
        file_name = files[file_id]
        
        # Handle obs_names alignment
        if not _ad.obs_names.equals(first_obs_names):
            if set(_ad.obs_names) != set(first_obs_names):
                missing_in_ad = set(first_obs_names) - set(_ad.obs_names)
                extra_in_ad = set(_ad.obs_names) - set(first_obs_names)
                raise ValueError(
                    f"obs_names mismatch between '{first_key}' and '{file_id}'.\n"
                    f"Missing in '{file_id}': {len(missing_in_ad)} (e.g., {list(missing_in_ad)[:5]})\n"
                    f"Extra in '{file_id}': {len(extra_in_ad)} (e.g., {list(extra_in_ad)[:5]})"
                )
            
            # Reload large slots before reordering
            _ad_tmp = read_anndata(
                file_name,
                **{slot: slot for slot in large_merge_slots}, 
                dask=True,
                backed=True,
                verbose=False
            )
            for slot in large_merge_slots:
                setattr(_ad, f'_{slot}', getattr(_ad_tmp, f'_{slot}'))
                file_copy_map[slot].append(file_id)
            
            logging.info(f"Reordering obs_names for {file_id}")
            adatas[file_id] = _ad[first_obs_names, :]
        
        # Handle var_names alignment  
        if not _ad.var_names.equals(first_var_names):
            if set(_ad.var_names) != set(first_var_names):
                raise ValueError(f"var_names mismatch between '{first_key}' and '{file_id}'")
            
            # Reload layers if needed
            if 'layers' in merge_slots:
                _ad.layers = read_anndata(
                    file_name,
                    layers='layers',
                    dask=True,
                    backed=True
                ).layers
                file_copy_map['layers'].append(file_id)
            
            logging.info(f"Reordering var_names for {file_id}")
            adatas[file_id] = adatas[file_id][:, first_var_names]
        
        # Convert views to copies
        if adatas[file_id].is_view:
            adatas[file_id] = adatas[file_id].copy()
    
    return file_copy_map

