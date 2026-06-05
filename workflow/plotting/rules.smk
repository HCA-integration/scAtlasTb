from utils.environments import get_env


rule plots:
    input:
        anndata='{dataset}.h5ad',
    output:
        plots=directory('{dataset}_plot'),
    params:
        basis='X_pca'
    threads: 10
    conda:
        get_env(config, 'scanpy')
    script:
        '../preprocessing/scripts/plot.py'
