def get_umap_file(wildcards):
    if mcfg.get_from_parameters(wildcards, 'recompute_umap', default=False):
        return rules.clustering_compute_umap.output.zarr
    return mcfg.get_input_file(**wildcards)


use rule umap from preprocessing as clustering_compute_umap with:
    input:
        zarr=get_neighbors_file,
        rep=get_neighbors_file,
    output:
        zarr=directory(mcfg.out_dir / 'umap' / f'{paramspace.wildcard_pattern}.zarr'),
    params:
        args=lambda wildcards: dict(
            neighbors_key=mcfg.get_from_parameters(wildcards, 'neighbors_key', default='neighbors'),
        ),
    threads: 10
    resources:
        partition=lambda w, attempt: mcfg.get_resource(profile='gpu',resource_key='partition',attempt=attempt, attempt_to_cpu=2),
        qos=lambda w, attempt: mcfg.get_resource(profile='gpu',resource_key='qos',attempt=attempt, attempt_to_cpu=2),
        mem_mb=lambda w, attempt: mcfg.get_resource(profile='gpu',resource_key='mem_mb',attempt=attempt, attempt_to_cpu=2),
        gpu=lambda w, attempt: mcfg.get_resource(profile='gpu',resource_key='gpu',attempt=attempt, attempt_to_cpu=2),


use rule merge from clustering as clustering_merge with:
    input:
        zarr=get_umap_file,
        cluster_anno=lambda wildcards: mcfg.get_output_files(
            rules.clustering_cluster.output.zarr,
            subset_dict=dict(wildcards),
            all_params=True
        ),
    output:
        zarr=directory(mcfg.out_dir / f'{paramspace.wildcard_pattern}.zarr'),
    params:
        neighbors_key=lambda wildcards: mcfg.get_from_parameters(wildcards, 'neighbors_key', default='neighbors'),
