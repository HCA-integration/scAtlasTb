import glob
import yaml
import logging
import pandas as pd
logging.basicConfig(level=logging.INFO)

from utils.io import read_anndata
from pandas.testing import assert_frame_equal


def load_config():
    """Load test configuration from YAML file."""
    return yaml.safe_load(open('test/config.yaml', 'r'))


def find_output_files():
    """Find zarr output files from collect step."""
    outputs = glob.glob('test/out/collect/dataset~*.zarr')
    assert len(outputs) > 0, 'No collect output files found'
    return outputs


def extract_dataset_name(file, config):
    """Extract dataset name from zarr file path using config datasets."""
    matching_datasets = [x for x in config['DATASETS'] if f'dataset~{x}.zarr' in file]
    return matching_datasets[0] if matching_datasets else None


def run_basic_structure_tests(adata, file):
    """Validate basic AnnData structure and index uniqueness."""
    assert adata.n_obs > 0, f"No observations in {file}"
    assert adata.n_vars > 0, f"No variables in {file}"
    assert adata.obs_names.is_unique, f"Non-unique obs_names in {file}"
    assert adata.var_names.is_unique, f"Non-unique var_names in {file}"


def check_obs_merging(adata, input_files, sep):
    """Validate obs merging when multiple input files are combined and check barcode consistency."""
    obs = adata.obs
    obs_columns = list(obs.columns)

    # Identify file-specific columns (those containing separator)
    file_specific_cols = [col for col in obs_columns if sep in col]
    assert file_specific_cols, f"No file-specific columns found with separator '{sep}'"

    # Map file id -> list of its suffixed columns
    problems = []
    checks = 0
    for file_id in input_files.keys():
        orig_obs = read_anndata(input_files[file_id], obs='obs', verbose=False).obs
        file_cols = [col for col in obs_columns if col.endswith(f'{sep}{file_id}')]
        logging.info(f"File {file_id} contributed {len(file_cols)} columns")

        # Build suffixed -> base mapping and reconstructed dataframe for this file
        stripped_map = {col: col.rsplit(sep, 1)[0] for col in file_cols if sep in col}
        if not stripped_map:
            continue

        df_suffixed = obs[list(stripped_map.keys())].rename(columns=stripped_map)
        missing_in_orig = set(df_suffixed.columns) - set(orig_obs.columns)
        if missing_in_orig:
            problems.append(f"Missing base columns in original obs for file_id={file_id}: {sorted(missing_in_orig)}")
            continue

        df_orig = orig_obs[df_suffixed.columns]

        # Align indices (allow different order; require same set)
        assert set(df_orig.index) == set(df_suffixed.index), f"Index sets differ for file_id={file_id}"

        df_orig_sorted = df_orig.sort_index()
        df_suffixed_sorted = df_suffixed.sort_index()

        # Compare with tolerance for numeric columns, ignore column order
        try:
            assert_frame_equal(
                df_orig_sorted,
                df_suffixed_sorted,
                check_like=True,
                atol=1e-9,
                rtol=1e-9,
                check_dtype=False,
            )
        except AssertionError as e:
            problems.append(f"DataFrame mismatch for file_id={file_id}: {e}")
        else:
            checks += len(df_suffixed.columns)

    assert not problems, "Inconsistencies between base and suffixed columns:\n" + "\n".join(problems)
    logging.info(f"Checked {checks} base↔suffixed column pairs for value consistency")
    logging.info(f"✅ Obs values correctly mapped to same barcodes")


def process_dataset(file, config):
    """Run all validation tests on a single dataset output file."""
    dataset = extract_dataset_name(file, config)
    if not dataset:
        logging.info(f'No matching dataset found for {file}, skipping...')
        return
    
    logging.info(f'Check dataset "{dataset}", file {file}')
    
    collect_config = config['DATASETS'][dataset].get('collect', {})
    input_files = config['DATASETS'][dataset]['input']['collect']
    
    try:
        # Load merged output data
        adata = read_anndata(file, verbose=False)
        
        # Run validation tests
        run_basic_structure_tests(adata, file)
        
        # Get merge configuration parameters
        merge_slots = collect_config.get('merge_slots', [])
        sep = collect_config.get('sep', '--')
        
        # Test obs merging consistency if multiple input files were merged
        if 'obs' in merge_slots and isinstance(input_files, dict) and len(input_files) > 1:
            check_obs_merging(adata, input_files, sep)
        
        logging.info(f"✅ All collect assertions passed for dataset={dataset}")
        
    except Exception as e:
        logging.error(f"❌ Assertion failed for {dataset}: {e}")
        raise e


def main():
    """Main entry point to validate all collect step outputs."""
    config = load_config()
    outputs = find_output_files()
    
    for file in outputs:
        process_dataset(file, config)
    
    logging.info("🎉 All collect assertion tests passed!")


if __name__ == "__main__":
    main()
