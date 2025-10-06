# TODO: get input files

rule scarches:
    input:
        zarr=lambda w: mcfg.get_input_file(**w),
    output:
        zarr=directory(mcfg.out_dir / f'{paramspace.wildcard_pattern}.zarr'),
        # model=,
    conda:
        get_env(config, 'scarches', env_dir='envs')
    resources:
        partition=lambda w: mcfg.get_resource(profile='gpu',resource_key='partition'),
        qos=lambda w: mcfg.get_resource(profile='gpu',resource_key='qos'),
        mem_mb=lambda w: mcfg.get_resource(profile='gpu',resource_key='mem_mb'),
    script:
        '../scripts/majority_voting.py'