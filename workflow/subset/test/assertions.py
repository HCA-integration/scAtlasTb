import yaml
import glob
from utils.io import read_anndata


with open('test/config.yaml', 'r') as file:
    config = yaml.safe_load(file)

files = glob.glob('test/out/subset/*/*.zarr')
print(files)

for file in files:
    dataset = [x for x in config['DATASETS'] if x in file][0]
    n_cells_threshold = config['DATASETS'][dataset]['subset']['n_cells']

    print(f'read {file}...')
    adata = read_anndata(file, obs='obs', uns='uns')

    print(dataset)
    assert 'subset' in adata.uns
    n_cells = adata.obs.shape[0]
    print('\tthreshold: ', n_cells_threshold)
    print('\tadata.n_obs: ', n_cells)
    # the n_cells is the lower threshold when subsetting to complete samples
    # assert n_cells < n_cells_threshold, f'Number of cells {n_cells} exceeds threshold {n_cells_threshold} for dataset {dataset}!'
