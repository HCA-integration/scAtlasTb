import logging
logging.basicConfig(level=logging.INFO)
import pandas as pd
import numpy as np
from tqdm import tqdm


def mode_counts_per_row(arr, dtype=np.uint8):            
    n_rows, n_cols = arr.shape
    
    if n_rows > 255:
        dtype = np.int16

    modes = np.empty(n_rows, dtype=arr.dtype)
    mode_counts = np.empty(n_rows, dtype=dtype)

    miniters = max(1, int(n_rows) / 10)
    for i in tqdm(range(n_rows), desc='Compute majority vote per row', miniters=miniters):
        unique_vals, counts = np.unique(arr[i], return_counts=True)
        max_idx = np.argmax(counts)
        modes[i] = unique_vals[max_idx]
        mode_counts[i] = counts[max_idx]
    
    # shift negative values to non-negative to work with bincount
    arr_min = arr.min()
    shift_value = 0 if arr_min >= 0 else -arr_min
    arr = arr + shift_value
    
    # initialize mode array and mode counts
    max_val = arr.max()
    modes = np.zeros(n_rows, dtype=arr.dtype)
    mode_counts = np.zeros(n_rows, dtype=np.uint8)

    for i in range(n_rows):
        counts = np.bincount(arr[i], minlength=max_val + 1)
        modes[i] = np.argmax(counts)
        mode_counts[i] = np.max(counts)

    # restore original values of the modes
    modes = modes - shift_value

    return modes, mode_counts


def check_same_categories(df, columns):
    first_categories = df[columns[0]].cat.categories
    if not all(df[col].cat.categories.equals(first_categories) for col in columns):
        logging.warn(f"Columns {columns} have different categories, majority voting might not work correctly")


def get_majority_reference(df, reference_key, query_key, **kwargs):
    map_majority = pd.crosstab(df[reference_key], df[query_key], **kwargs).idxmax(axis=0)
    return pd.Categorical(
        df[query_key].map(map_majority),
        categories=map_majority.dropna().unique(),
    )

def get_majority_consensus(df, column_patterns, new_key='majority_consensus'):
    import re
    from functools import reduce

    # parse regex columns
    columns = [
        col for col in df.columns for pattern in column_patterns
        if re.fullmatch(pattern, col)
    ]
    columns = list(set(columns))
    assert len(columns) >= 2, f'No columns found matching patterns: {column_patterns}'
    
    agreement_col = f'{new_key}_agreement'
    low_agreement_col = f'{new_key}_low_agreement'
    n_vote_columns = len(columns)
    
    check_same_categories(df, columns)
    
    categories = reduce(pd.Index.union, [df[col].cat.categories.dropna() for col in columns])
    cat_dtype = pd.CategoricalDtype(categories=categories)
    df_cat = df[columns].astype(cat_dtype).apply(lambda x: x.cat.codes)
    
    majority_votes, mode_count = mode_counts_per_row(df_cat.to_numpy())
    
    logging.info('Assign majority votes to new column...')
    df[new_key] = pd.Categorical.from_codes(majority_votes, categories=categories)
    df[agreement_col] = mode_count / n_vote_columns
    
    min_agreement = 1 / n_vote_columns
    df.loc[df[agreement_col] <= min_agreement, new_key] = float('nan')
    # set to 0, because no agreement when all assignments different
    df.loc[df[agreement_col] <= min_agreement, agreement_col] = 0
    df[low_agreement_col] = df[agreement_col] <= 0.5
    
    return df[[new_key, agreement_col, low_agreement_col]]
