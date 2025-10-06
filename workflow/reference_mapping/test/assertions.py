from anndata.experimental import read_elem
import zarr
import glob
import yaml
import pandas as pd


files = glob.glob('test/out/label_transfer/dataset~*/file_id~*.zarr')
if not files:
    print('No files found to test')

config = yaml.safe_load(open('test/config.yaml', 'r'))

for file in files:
    dataset = [x for x in config['DATASETS'] if x in file][0]
    module_config  = config['DATASETS'][dataset]['label_transfer']
    
    print(f'read {file}...')
    with zarr.open(file) as z:
        obs = read_elem(z["obs"])
    
    print('assert majority reference...')
    key = 'majority_reference'
    maj_ref = module_config.get('majority_reference', {})
    reference = maj_ref.get('reference_key')
    cluster = maj_ref.get('cluster_key')
    categories = obs[reference].unique()
    for ref in obs[key].unique():
        assert ref in categories, f'{ref} not in "{key}":\n {categories}'