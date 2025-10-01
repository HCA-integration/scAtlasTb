from anndata.experimental import read_elem
import zarr
import re
import glob
import yaml
import pandas as pd


files = glob.glob('test/out/majority_voting/dataset~*/file_id~*.zarr')
if not files:
    print('No files found to test')

config = yaml.safe_load(open('test/config.yaml', 'r'))

for file in files:
    dataset = [x for x in config['DATASETS'] if x in file][0]
    module_config  = config['DATASETS'][dataset]['majority_voting']
    
    print(f'read {file}...')
    with zarr.open(file) as z:
        obs = read_elem(z["obs"])    
        
    print('assert majority consensus...')
    key = 'majority_consensus'
    column_patterns = module_config.get('columns', [])
    columns = [
        col for col in obs.columns for pattern in column_patterns
        if re.fullmatch(pattern, col)
    ]
    columns = list(dict.fromkeys(columns))
    print('columns:', columns)
    
    categories = pd.concat([obs[col] for col in columns]).unique()
    for ref in obs[key].unique():
        assert ref in categories, f'"{ref}" not in:\n {categories}'
