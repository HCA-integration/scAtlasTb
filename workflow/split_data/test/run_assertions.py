import glob
from pprint import pprint
import numpy as np
from anndata.experimental import read_elem
import zarr
import logging
logging.basicConfig(level=logging.INFO)


outputs = glob.glob('test/out/split_data/dataset~*/file_id~*/key~*/value~*.zarr')
assert len(outputs) > 0, 'No output files found'

for file in outputs:
    logging.info(f'Checking {file}...')
    z = zarr.open(file)
    uns = read_elem(z["uns"])

    try:
        # Check Metadata
        assert 'wildcards' in uns
        assert 'split_data_key' in uns['wildcards']
        assert 'split_data_value' in uns['wildcards']
        assert 'X' in z
        assert 'obs' in z
        assert 'var' in z

    except AssertionError as e:
        # pprint(uns.keys())
        print('AssertionError for:', file)
        raise e
