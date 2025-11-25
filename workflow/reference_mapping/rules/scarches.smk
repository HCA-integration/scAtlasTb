rule scarches:
    input:
        zarr=lambda wildcards: mcfg.get_input_file(**wildcards),
        model=lambda wildcards: mcfg.get_from_parameters(wildcards, 'scarches', default={}).get('model')
    output:
        zarr=directory(mcfg.out_dir / f'{paramspace.wildcard_pattern}.zarr'),
        model=directory(mcfg.out_dir / 'model' / paramspace.wildcard_pattern),
    params:
        layer=lambda wildcards: mcfg.get_from_parameters(wildcards, 'scarches', default={}).get('layer', 'X'),
        var_key=lambda wildcards: mcfg.get_from_parameters(wildcards, 'scarches', default={}).get('var_key'),
        model_params=lambda wildcards: mcfg.get_from_parameters(wildcards, 'scarches', default={}).get('model_params', {}),
        train_params=lambda wildcards: mcfg.get_from_parameters(wildcards, 'scarches', default={}).get('train_params', {}),
    threads:
        lambda wildcards: mcfg.get_from_parameters(wildcards, 'threads', default=10, as_type=int)
    conda:
        get_env(config, 'scarches', env_dir='envs')
    resources:
        partition=lambda wildcards: mcfg.get_resource(profile='gpu',resource_key='partition'),
        gpu=lambda wildcards: mcfg.get_resource(profile='gpu',resource_key='gpu'),
        qos=lambda wildcards: mcfg.get_resource(profile='gpu',resource_key='qos'),
        mem_mb=lambda wildcards: mcfg.get_resource(profile='gpu',resource_key='mem_mb'),
    script:
        '../scripts/scarches.py'