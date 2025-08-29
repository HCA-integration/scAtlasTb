rule plot_distances:
    input:
        zarr=rules.run_method.output.zarr
    output:
        histplot_path=mcfg.image_dir / paramspace.wildcard_pattern / 'distances_histplot.png'
    conda:
        get_env(config, 'sample_representation')
    resources:
        mem_mb=lambda w, attempt: mcfg.get_resource(profile='cpu', resource_key='mem_mb', attempt=attempt)
    script:
        '../scripts/distances_plot.py'


rule mds:
    input:
        zarr=rules.run_method.output.zarr
    output:
        zarr=directory(mcfg.out_dir / 'mds' / f'{paramspace.wildcard_pattern}.zarr'),
    conda:
        get_env(config, 'sample_representation')
    resources:
        mem_mb=lambda w, attempt: mcfg.get_resource(profile='cpu', resource_key='mem_mb', attempt=attempt)
    script:
        '../scripts/mds.py'


use rule plots from preprocessing as sample_representation_plot_mds with:
    input:
        zarr=rules.mds.output.zarr,
    output:
        plots=directory(mcfg.image_dir / paramspace.wildcard_pattern / 'mds')
    params:
        basis='X_mds',
        color=lambda wildcards: mcfg.get_from_parameters(wildcards, 'colors', default=[])
    resources:
        partition=mcfg.get_resource(profile='cpu',resource_key='partition'),
        qos=mcfg.get_resource(profile='cpu',resource_key='qos'),


# UMAP
use rule umap from preprocessing as sample_representation_compute_umap with:
    input:
        zarr=rules.run_method.output.zarr,
        rep=rules.run_method.output.zarr,
    output:
        zarr=directory(mcfg.out_dir / f'{paramspace.wildcard_pattern}.zarr'),
    resources:
        partition=mcfg.get_resource(profile='cpu',resource_key='partition'),
        qos=mcfg.get_resource(profile='cpu',resource_key='qos'),
        mem_mb=mcfg.get_resource(profile='cpu',resource_key='mem_mb'),
        gpu=mcfg.get_resource(profile='cpu',resource_key='gpu'),


use rule plots from preprocessing as sample_representation_plot_umap with:
    input:
        zarr=rules.sample_representation_compute_umap.output.zarr,
    output:
        plots=directory(mcfg.image_dir / paramspace.wildcard_pattern / 'umap')
    params:
        basis='X_umap',
        color=lambda wildcards: mcfg.get_from_parameters(wildcards, 'colors', default=[])
    resources:
        partition=mcfg.get_resource(profile='cpu',resource_key='partition'),
        qos=mcfg.get_resource(profile='cpu',resource_key='qos'),


# PCA

use rule plots from preprocessing as sample_representation_plot_emb with:
    input:
        zarr=rules.run_method.output.zarr,
    output:
        plots=directory(mcfg.image_dir / paramspace.wildcard_pattern / 'embeddings')
    params:
        basis='X_pca',
        color=lambda wildcards: mcfg.get_from_parameters(wildcards, 'colors', default=[])
    resources:
        partition=mcfg.get_resource(profile='cpu',resource_key='partition'),
        qos=mcfg.get_resource(profile='cpu',resource_key='qos'),
