import glob
from pprint import pformat
import yaml
from anndata.experimental import read_elem
import zarr
import logging
logging.basicConfig(level=logging.INFO)


config = yaml.safe_load(open('test/config.yaml', 'r'))

# direct integration outputs
outputs = glob.glob('test/out/clustering/dataset~*/file_id~*.zarr')
assert len(outputs) > 0, "No outputs found"

for file in outputs:
    dataset = [x for x in config['DATASETS'] if f'dataset~{x}/' in file][0]

    logging.info(f'dataset: {dataset}')
    logging.info(f'Checking {file} for...')
    z = zarr.open(file)

    clustering_config = config['DATASETS'][dataset]['clustering']
    logging.info(pformat(clustering_config))
    
    # Extract level and resolution from file path
    file_parts = file.split('/')
    level = None
    resolution = None
    for part in file_parts:
        if 'level~' in part:
            level = int(part.split('level~')[1].split('/')[0])
        if 'resolution~' in part:
            resolution = float(part.split('resolution~')[1].split('/')[0])
    
    if level is None or resolution is None:
        logging.warning(f"Could not extract level or resolution from {file}")
        continue
        
    algorithms = clustering_config.get('algorithm', ['leiden'])
    if isinstance(algorithms, str):
        algorithms = [algorithms]
        
    for algorithm in algorithms:
        # Clustering script creates columns with format: {algorithm}_{resolution}_{level}
        cluster_key = f'{algorithm}_{resolution}_{level}'
        assert cluster_key in z["obs"], f"Clustering column {cluster_key} not found in {file}. Available columns: {list(z['obs'].keys())}"
        
        # Check that clustering results contain valid cluster IDs
        cluster_data = read_elem(z["obs"][cluster_key])
        assert len(cluster_data) > 0, f"Clustering column {cluster_key} is empty"
        
        # Check for reasonable number of clusters
        unique_clusters = len(set(cluster_data))
        total_cells = len(cluster_data)
        assert unique_clusters > 0, f"No clusters found in {cluster_key}"
        assert unique_clusters < total_cells, f"Too many clusters ({unique_clusters}) for {total_cells} cells in {cluster_key}"
        
        logging.info(f"✓ Found {unique_clusters} clusters in {cluster_key} for {total_cells} cells")
        
        # For hierarchical clustering (level > 1), check naming format
        if level > 1:
            # Hierarchical clusters should have underscore-separated format
            sample_clusters = list(set(cluster_data))[:5]  # Check first 5 unique clusters
            for cluster_id in sample_clusters:
                if cluster_id != cluster_id:  # Skip NaN values
                    continue
                cluster_parts = str(cluster_id).split('_')
                assert len(cluster_parts) >= 2, f"Hierarchical cluster {cluster_id} should have format parent_child in {cluster_key}"
                9
        logging.info(f"✓ Clustering format validation passed for {cluster_key}")

logging.info("All clustering assertions passed!")
