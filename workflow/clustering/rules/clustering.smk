def get_neighbors_file(wildcards):
    level = int(wildcards.get('level', 1))
    if level > 1:
        return expand(
            rules.clustering_cluster.output.zarr,
            level=level - 1,
            allow_missing=True
        )[0]
    if mcfg.get_from_parameters(wildcards, 'recompute_neighbors', default=False):
        return rules.clustering_compute_neighbors.output.zarr
    return mcfg.get_input_file(**wildcards)


use rule neighbors from preprocessing as clustering_compute_neighbors with:
    input:
        zarr=lambda wildcards: mcfg.get_input_file(**wildcards),
    output:
        zarr=directory(mcfg.out_dir / 'neighbors' / f'{paramspace.wildcard_pattern}.zarr'),
    params:
        args=lambda wildcards: mcfg.get_from_parameters(
            {k: v for k, v in wildcards.items() if k not in ['algorithm', 'resolution']},
            'neighbors',
            default={},
            exclude=['algorithm', 'resolution'],
        ),
    resources:
        partition=lambda w, attempt: mcfg.get_resource(profile='gpu',resource_key='partition',attempt=attempt),
        qos=lambda w, attempt: mcfg.get_resource(profile='gpu',resource_key='qos',attempt=attempt),
        mem_mb=lambda w, attempt: mcfg.get_resource(profile='gpu',resource_key='mem_mb',attempt=attempt),
        gpu=lambda w, attempt: mcfg.get_resource(profile='gpu',resource_key='gpu',attempt=attempt),


def get_profile(wildcards):
    use_gpu = mcfg.get_from_parameters(wildcards, 'use_gpu', default=False)
    if use_gpu: # and wildcards.level == '1':
        return 'gpu'
    return 'cpu'


use rule cluster from clustering as clustering_cluster with:
    input:
        zarr=get_neighbors_file,
    output:
        zarr=directory(mcfg.out_dir / 'resolutions' / paramspace.wildcard_pattern / 'algorithm~{algorithm}' / 'resolution~{resolution}' / 'level~{level}.zarr'),
    params:
        neighbors_key=lambda wildcards: mcfg.get_from_parameters(wildcards, 'neighbors_key', default='neighbors'),
        neighbors_args=lambda wildcards: mcfg.get_from_parameters(wildcards, 'neighbors', default={}),
        clustering_args=lambda wildcards: mcfg.get_from_parameters(wildcards, 'kwargs', default={}),
        overwrite=lambda wildcards: mcfg.get_from_parameters(wildcards, 'overwrite', default=True),
    threads:
        lambda wildcards: 4 * int(wildcards.level) - 3
    conda:
        get_env(config, 'scanpy', gpu_env='rapids_singlecell')
    resources:
        partition=lambda w, attempt: mcfg.get_resource(profile=get_profile(w), resource_key='partition', attempt=attempt, attempt_to_cpu=2),
        qos=lambda w, attempt: mcfg.get_resource(profile=get_profile(w), resource_key='qos', attempt=attempt, attempt_to_cpu=2),
        mem_mb=lambda w, attempt: mcfg.get_resource(profile=get_profile(w), resource_key='mem_mb', attempt=attempt, attempt_to_cpu=2),
        gpu=lambda w, attempt: mcfg.get_resource(profile=get_profile(w), resource_key='gpu', attempt=attempt, attempt_to_cpu=2),
