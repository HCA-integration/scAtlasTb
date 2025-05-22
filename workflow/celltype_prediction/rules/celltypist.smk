model_dir = mcfg.out_dir / 'models' / 'celltypist'
os.environ["CELLTYPIST_FOLDER"] = str(model_dir)

envvars:
    "CELLTYPIST_FOLDER"


rule get_celltypist_model:
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
        model=rules.get_celltypist_model.output.model,
    output:
        zarr=directory(mcfg.out_dir / 'celltypist' / paramspace.wildcard_pattern / '{celltypist_model}.zarr'),
        png=directory(mcfg.image_dir / paramspace.wildcard_pattern / 'celltypist--{celltypist_model}'),
    conda:
        lambda wildcards: get_env(config, 'celltypist')
    params:
        CELLTYPIST_FOLDER=os.environ["CELLTYPIST_FOLDER"],
        label_key=lambda wildcards: mcfg.get_from_parameters(wildcards, 'reference_label'),
        celltypist_params=lambda wildcards: mcfg.get_from_parameters(wildcards, 'celltypist').get('params', {}),
        raw_layer=lambda wildcards: mcfg.get_from_parameters(wildcards, 'raw_counts', default=None),
        norm_layer=lambda wildcards: mcfg.get_from_parameters(wildcards, 'norm_counts', default=None),
        is_normalized=lambda wildcards: mcfg.get_from_parameters(wildcards, 'norm_counts', default=None) is not None,
    resources:
        mem_mb=mcfg.get_resource(profile='cpu',resource_key='mem_mb')
    script:
        '../scripts/celltypist.py'
