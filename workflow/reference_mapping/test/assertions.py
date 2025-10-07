from anndata.experimental import read_elem
import zarr
import glob
import yaml
import pandas as pd
import numpy as np


files = glob.glob('test/out/reference_mapping/dataset~*/file_id~*.zarr')
if not files:
    print('No files found to test')
    exit(1)

config = yaml.safe_load(open('test/config.yaml', 'r'))

for file in files:
    dataset = [x for x in config['DATASETS'] if x in file][0]

    ## scArches
    module_config = config['DATASETS'][dataset]['reference_mapping'].get('scarches', {})
    
    
    print(f'Testing reference mapping output: {file}...')
    
    with zarr.open(file) as z:
        obs = read_elem(z["obs"])
        obsm = read_elem(z["obsm"]) if "obsm" in z else {}
        var = read_elem(z["var"])
    

    print("  ✓ Testing reference mapping outputs...")
    
    # Check for latent embedding
    assert 'X_emb' in obsm, "Expected 'X_emb' in obsm keys"
    
    # Check latent embedding dimensions
    latent_embedding = obsm['X_emb']
    assert latent_embedding.shape[0] == len(obs), f"Latent embedding rows should match obs length"
    assert latent_embedding.shape[1] > 0, f"Latent embedding should have > 0 dimensions"
    assert not np.any(np.isnan(latent_embedding)), f"Latent embedding should not contain NaN values"

    
    print("  ✓ Testing model parameter consistency...")
    model_params = module_config.get('model_params', {})
    
    # Check batch key if specified
    if 'batch_key' in model_params:
        batch_col = model_params['batch_key']
        if batch_col:  # Only check if not None/empty
            assert batch_col in obs.columns, f"Batch key '{batch_col}' not found in obs columns"
            assert obs[batch_col].notna().all(), f"Batch column '{batch_col}' should not contain NaN values"
    
    # Check labels key if specified  
    if 'labels_key' in model_params:
        labels_col = model_params['labels_key']
        if labels_col:  # Only check if not None/empty
            # Labels might be in original form or as predicted labels
            possible_label_cols = [labels_col, f"predicted_{labels_col}", "predicted_labels"]
            found_labels = [col for col in possible_label_cols if col in obs.columns]
            assert len(found_labels) > 0, f"No label columns found. Expected one of: {possible_label_cols}"
    
    # Check covariate columns if specified
    for covariate_type in ['categorical_covariate', 'continuous_covariate']:
        if covariate_type in model_params:
            covariates = model_params[covariate_type]
            if covariates:  # Only check if not None/empty
                for covariate in covariates:
                    assert covariate in obs.columns, f"Covariate '{covariate}' not found in obs columns"
    
    
    print("  ✓ Testing data integrity...")
    
    # Check for reasonable cell and gene counts
    assert len(obs) > 0, "Should have > 0 cells"
    assert len(var) > 0, "Should have > 0 genes"
    assert len(obs) < 1_000_000, "Suspiciously high number of cells (>1M)"
    
    print(f"  ✅ All tests passed for {file}")

# Check for model output directory
model_dirs = glob.glob('test/out/reference_mapping/_model/dataset~*/file_id~*')
if model_dirs:
    print(f"\n✓ Found {len(model_dirs)} model output directories")
    for model_dir in model_dirs:
        print(f"  Checking model directory: {model_dir}")
        # Check for typical scVI-tools model files
        expected_files = ['model.pt', 'attr.pkl', 'var_names.csv']
        for expected_file in expected_files:
            model_file_path = f"{model_dir}/{expected_file}"
            import os
            if os.path.exists(model_file_path):
                print(f"    ✓ Found {expected_file}")
            else:
                print(f"    ⚠ Missing {expected_file} (may be optional)")

print(f"\n🎉 Reference mapping tests completed successfully for {len(files)} output files!")

