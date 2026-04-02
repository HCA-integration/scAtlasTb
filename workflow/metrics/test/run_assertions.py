import glob
import yaml
import logging
from pathlib import Path

import pandas as pd

logging.basicConfig(level=logging.INFO)


METRICS_COLUMNS = ['dataset', 'metric', 'output_type', 'metric_type', 'score']


def load_config():
    """Load test configuration from YAML file."""
    return yaml.safe_load(open('test/config.yaml', 'r'))


def find_metrics_files():
    """Find metrics output files."""
    metrics_file = Path('test/out/metrics/results/metrics.tsv')
    assert metrics_file.exists(), f'No metrics file found at {metrics_file}'
    return metrics_file


def assert_required_columns(metrics_df: pd.DataFrame, required_columns: list[str]) -> None:
    """Validate that required columns exist in metrics table."""
    missing = [column for column in required_columns if column not in metrics_df.columns]
    assert not missing, (
        f"Missing required columns: {missing}\n"
        f"Available columns: {list(metrics_df.columns)}"
    )
    logging.info(f"✅ All required columns present: {required_columns}")


def assert_no_missing_values(metrics_df: pd.DataFrame, column: str) -> None:
    """Check that a column has no missing values."""
    n_missing = metrics_df[column].isna().sum()
    assert n_missing == 0, f"Column '{column}' has {n_missing} missing values"


def validate_metrics_content(metrics_df: pd.DataFrame, config: dict) -> None:
    """Validate metrics table content and data types."""
    
    # Check that table is not empty
    assert len(metrics_df) > 0, 'Metrics table is empty'
    logging.info(f"Metrics table contains {len(metrics_df)} rows")
    
    # Check no missing values in critical columns (except score which can be NA for failed metrics)
    for col in ['dataset', 'metric']:
        assert_no_missing_values(metrics_df, col)
    
    # Check that score column is numeric where present
    assert pd.api.types.is_numeric_dtype(metrics_df['score']), \
        f"'score' column should be numeric, got {metrics_df['score'].dtype}"
    
    # Log metrics with missing scores
    missing_scores = metrics_df[metrics_df['score'].isna()]
    if len(missing_scores) > 0:
        logging.warning(f"{len(missing_scores)} metrics have missing scores")
        logging.warning(f"Metrics with NA scores: {missing_scores['metric'].unique().tolist()}")
    
    # Check that datasets match config
    config_datasets = set(config['DATASETS'].keys())
    metrics_datasets = set(metrics_df['dataset'].unique())
    unexpected = metrics_datasets - config_datasets
    if unexpected:
        logging.warning(f"Unexpected datasets in metrics: {unexpected}")
    
    # Log summary statistics for non-NA scores
    valid_scores = metrics_df[metrics_df['score'].notna()]
    if len(valid_scores) > 0:
        logging.info(f"Valid scores: {len(valid_scores)}/{len(metrics_df)}")
        logging.info(f"Score range: [{valid_scores['score'].min():.3f}, {valid_scores['score'].max():.3f}]")
    
    logging.info(f"Unique datasets: {len(metrics_datasets)}")
    logging.info(f"Unique metrics: {len(metrics_df['metric'].unique())}")


def validate_per_dataset_files(config: dict) -> None:
    """Validate that per-dataset metrics files were created."""
    dataset_files = glob.glob('test/out/metrics/results/dataset~*/label~*/metrics.tsv')
    
    if len(dataset_files) == 0:
        logging.warning("No per-dataset metrics files found")
        return
    
    logging.info(f"Found {len(dataset_files)} per-dataset metrics files")
    
    for file_path in dataset_files:
        file_path = Path(file_path)
        df = pd.read_table(file_path)
        
        # Basic validation
        assert len(df) > 0, f"Empty metrics file: {file_path}"
        assert_required_columns(df, METRICS_COLUMNS)
        
        logging.info(f"✅ Valid: {file_path.relative_to('test/out/metrics/results')}")


def main() -> None:
    """Main entry point to validate metrics outputs."""
    
    config = load_config()
    metrics_file = find_metrics_files()
    
    logging.info(f'Checking metrics file: {metrics_file}')
    
    try:
        # Load metrics table
        metrics_df = pd.read_table(metrics_file)
        
        # Run validation tests
        assert_required_columns(metrics_df, METRICS_COLUMNS)
        validate_metrics_content(metrics_df, config)
        
        # Validate per-dataset files
        validate_per_dataset_files(config)
        
        logging.info('🎉 All metrics assertion tests passed!')
        
    except AssertionError as e:
        logging.error(f"❌ Metrics assertion failed: {e}")
        if 'metrics_df' in locals():
            logging.error(f"\nFirst few rows:\n{metrics_df.head()}")
        raise e


if __name__ == '__main__':
    main()