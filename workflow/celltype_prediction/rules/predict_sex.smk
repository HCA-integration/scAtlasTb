rule predict_sex:
    input:
        zarr=lambda wildcards: mcfg.get_input_file(**wildcards),
    output:
        zarr=directory(mcfg.out_dir / 'predict_sex' / f'{paramspace.wildcard_pattern}.zarr'),
    conda:
        get_env(config, 'scanpy', gpu_env='scanpy')
    params:
        layer=lambda wildcards: mcfg.get_from_parameters(wildcards, 'counts', default=None),
        is_normalized=lambda wildcards: mcfg.get_from_parameters(wildcards, 'is_normalized', default=True),
        params=lambda wildcards: mcfg.get_from_parameters(wildcards, 'predict_sex', default={}),
    resources:
        partition=lambda w, attempt: mcfg.get_resource(resource_key='partition', profile='cpu', attempt=attempt),
        qos=lambda w, attempt: mcfg.get_resource(resource_key='qos', profile='cpu', attempt=attempt),
        mem_mb=lambda w, attempt: mcfg.get_resource(resource_key='mem_mb', profile='cpu', attempt=attempt, factor=3),
        gpu=lambda w, attempt: mcfg.get_resource(resource_key='gpu', profile='cpu', attempt=attempt),
    script:
        '../scripts/predict_sex.py'


def get_predict_sex_files(wildcards):
    return {
        dataset: mcfg.get_output_files(
            rules.predict_sex.output.zarr,
            subset_dict=dict(dataset=dataset),
            all_params=True
        ) for dataset in mcfg.get_datasets()
        if mcfg.get_from_parameters(dict(dataset=dataset), 'predict_sex')
    }


rule predict_sex_all:
    input: unpack(get_predict_sex_files)
    localrule: True