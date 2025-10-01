rule majority_voting:
    input:
        zarr=lambda w: mcfg.get_input_file(**w),
    output:
        zarr=directory(mcfg.out_dir / f'{paramspace.wildcard_pattern}.zarr'),
    params:
        reference_key=lambda w: mcfg.get_from_parameters(w, 'majority_reference').get('reference_key'),
        query_key=lambda w: mcfg.get_from_parameters(w, 'majority_reference').get('query_key'),
        crosstab_kwargs=lambda w: mcfg.get_from_parameters(w, 'majority_reference').get('crosstab_kwargs', {}),
    conda:
        get_env(config, 'scanpy', env_dir='envs')
    resources:
        partition=lambda w: mcfg.get_resource(profile='cpu',resource_key='partition'),
        qos=lambda w: mcfg.get_resource(profile='cpu',resource_key='qos'),
        mem_mb=lambda w: mcfg.get_resource(profile='cpu',resource_key='mem_mb'),
    script:
        '../scripts/majority_voting.py'