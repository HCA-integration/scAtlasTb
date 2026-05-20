from utils.environments import get_env


rule pseudobulk:
    input:
        zarr='dataset~{dataset}.zarr',
    output:
        zarr=directory('dataset~{dataset}/per_cell.zarr'),
        bulks=directory('dataset~{dataset}/pseudobulk.zarr'),
    conda:
        get_env(config, 'scanpy')
    script:
        '../scripts/prepare.py'
