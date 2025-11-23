import pandas as pd
import numpy as np
from tqdm import tqdm
from pprint import pformat


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
    def _check_columns_equal(s1, s2):
        assert s1.index.equals(s2.index), \
            'Indices do not match for column comparison'
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
    if all(_check_columns_equal(_ad.obs[col], _ad1.obs[col]) for _ad in iterator):
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
