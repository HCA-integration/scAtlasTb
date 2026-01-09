model_dir = mcfg.out_dir / 'models' / 'celltypist'
os.environ["CELLTYPIST_FOLDER"] = str(model_dir)

envvars:
    "CELLTYPIST_FOLDER"


rule celltypist_get_model:
    output:
        model=model_dir / '{celltypist_model}.pkl'
    conda:
        lambda wildcards: get_env(config, 'celltypist')
    params:
        CELLTYPIST_FOLDER=os.environ["CELLTYPIST_FOLDER"]
    script:
        '../scripts/get_celltypist_model.py'


rule celltypist:
    input:
        zarr=lambda wildcards: mcfg.get_input_file(**wildcards),
        model=rules.celltypist_get_model.output.model,
    output:
        zarr=directory(mcfg.out_dir / 'celltypist' / paramspace.wildcard_pattern / '{celltypist_model}.zarr'),
        png=directory(mcfg.image_dir / paramspace.wildcard_pattern / 'celltypist--{celltypist_model}'),
    conda:
        lambda wildcards: get_env(config, 'celltypist')
    params:
        CELLTYPIST_FOLDER=os.environ["CELLTYPIST_FOLDER"],
        label_key=lambda wildcards: mcfg.get_from_parameters(wildcards, 'reference_label'),
        celltypist_params=lambda wildcards: mcfg.get_from_parameters(wildcards, 'celltypist').get('params', {}),
        layer=lambda wildcards: mcfg.get_from_parameters(wildcards, 'counts', default=None),
        is_normalized=lambda wildcards: mcfg.get_from_parameters(wildcards, 'is_normalized', default=False),
    resources:
        partition=lambda w, attempt: mcfg.get_resource(resource_key='partition', profile='cpu', attempt=attempt),
        qos=lambda w, attempt: mcfg.get_resource(resource_key='qos', profile='cpu', attempt=attempt),
        mem_mb=lambda w, attempt: mcfg.get_resource(resource_key='mem_mb', profile='cpu', attempt=attempt, factor=3),
        gpu=lambda w, attempt: mcfg.get_resource(resource_key='gpu', profile='cpu', attempt=attempt),
    script:
        '../scripts/celltypist.py'


rule celltypist_all:
    input: 
        celltypist=mcfg.get_output_files(rules.celltypist.output.zarr, all_params=True),
    localrule: True