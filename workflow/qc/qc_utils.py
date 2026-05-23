import pandas as pd


def read_threshold_file(file: str):
    """
    Read user threshold file in TSV format
    """
    df = pd.read_table(file, comment='#')
    assert 'file_id' in df.columns
    prefixes = ['percent_mito', 'n_genes', 'n_counts']
    if all(not col.startswith(prefix) for col in df.columns for prefix in prefixes):
        logging.warning(f'WARNING: None of expected QC stat columns found in {file}.')
    columns = [col for col in df.columns if any(col.startswith(x) for x in prefixes) or col in ['file_id', 'threshold_type']]
    if 'threshold_type' in columns:
        df = df[df['threshold_type'].isin(['user', 'alternative'])]
    else:
        df['threshold_type'] = 'user'
    return df[columns].drop_duplicates()


def unpack_thresholds(row: pd.Series):
    thresholds = row.thresholds
    if isinstance(thresholds, str):
        thresholds = ast.literal_eval(thresholds)
    if isinstance(thresholds, dict):
         thresholds = thresholds.get(row.file_id, thresholds)
    return thresholds
