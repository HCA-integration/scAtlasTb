import glob
import yaml
import logging
logging.basicConfig(level=logging.INFO)

from utils.io import read_anndata

config = yaml.safe_load(open('test/config.yaml', 'r'))

# Find output files from filter step
single_outputs = glob.glob('test/out/filter/dataset~*/file_id~*.zarr')
assert len(single_outputs) > 0, 'No output files found'

for file in single_outputs:
    matching_datasets = [x for x in config['DATASETS'] if f'dataset~{x}/' in file]
    if not matching_datasets:
        logging.info(f'No matching dataset found for {file}, skipping...')
        continue
    
    dataset = matching_datasets[0]
    logging.info(f'Check dataset "{dataset}", file {file}')

    filter_config = config['DATASETS'][dataset]['filter']
    input_file = config['DATASETS'][dataset]['input']['filter']
    
    # Load data
    adata_out = read_anndata(file, verbose=False)
    adata_in = read_anndata(input_file, verbose=False)
    
    try:
        # Check 'filtered' column exists and is boolean
        assert 'filtered' in adata_out.obs.columns, "'filtered' column missing"
        assert adata_out.obs['filtered'].dtype == bool, "'filtered' column should be boolean"
        
        # Check subsetting behavior
        subset = filter_config.get('subset', True)
        if subset:
            # Data should be subsetted to only non-filtered cells
            assert adata_out.obs['filtered'].all(), "All cells in subsetted output should have filtered=True"
            assert adata_out.n_obs <= adata_in.n_obs, "Output has more cells than input"
            logging.info(f'Input cells: {adata_in.n_obs}, Output cells: {adata_out.n_obs}')
        else:
            # Data should retain all cells but have filtered column
            assert adata_out.n_obs == adata_in.n_obs, "subset=False but cell count changed"
        
        # # TODO: Check that wildcards exist in uns
        # assert 'wildcards' in adata_out.uns, "wildcards missing in uns"
        
        logging.info(f"✅ All assertions passed for dataset={dataset}")
        
    except AssertionError as e:
        logging.error(f"❌ Assertion failed for {dataset}: {e}")
        logging.error(adata_out)
        raise e

logging.info("🎉 All filter assertion tests passed!")