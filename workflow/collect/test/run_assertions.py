import glob
import yaml
import logging
import pandas as pd
logging.basicConfig(level=logging.INFO)

from utils.io import read_anndata


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

        for suffixed in file_cols:
            if sep not in suffixed:
                continue

            base = suffixed.rsplit(sep, 1)[0]
            assert base in orig_obs.columns, f"Base column '{base}' not found in original obs for file_id={file_id}"

            s_base = orig_obs[base]
            s_file = obs[suffixed]
            
            assert s_base.index.equals(s_file.index), f"Index mismatch between base '{base}' and suffixed '{suffixed}' for file_id={file_id}"

            mask = (~s_base.isna()) & (~s_file.isna())
            if mask.sum() == 0:
                continue
            
            v_base = s_base[mask]
            v_file = s_file[mask]
            
            if pd.api.types.is_numeric_dtype(v_base) and pd.api.types.is_numeric_dtype(v_file):
                unequal = ~((v_base.astype(float) - v_file.astype(float)).abs() <= 1e-9)
            else:
                unequal = v_base.astype(str) != v_file.astype(str)
            
            if unequal.any():
                example = list(unequal[unequal].index[:10])
                problems.append(
                    f"Mismatch base '{base}' vs '{suffixed}' (file_id={file_id}, {unequal.sum()} differing barcodes, examples: {example})"
                )
            
            checks += 1

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
