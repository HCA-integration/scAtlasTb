def get_file(wildcards, rule):
    params = mcfg.get_for_dataset(
        dataset=wildcards.dataset,
        query=[mcfg.module_name],
    )
    if rule in params:
        return getattr(rules, f'batch_analysis_{rule}').output.zarr
    return mcfg.get_input_file(**wildcards)


use rule normalize from preprocessing as batch_analysis_normalize with:
    input:
        zarr=lambda wildcards: mcfg.get_input_file(**wildcards),
    output:
        zarr=directory(mcfg.out_dir / 'prepare' / paramspace.wildcard_pattern / 'normalize.zarr'),
    params:
        args=lambda wildcards: mcfg.get_from_parameters(wildcards, 'normalize', default={}),
    resources:
        partition=lambda w, attempt: mcfg.get_resource(profile='gpu',resource_key='partition',attempt=attempt),
        qos=lambda w, attempt: mcfg.get_resource(profile='gpu',resource_key='qos',attempt=attempt),
        mem_mb=lambda w, attempt: mcfg.get_resource(profile='gpu',resource_key='mem_mb',attempt=attempt),
        gpu=lambda w, attempt: mcfg.get_resource(profile='gpu',resource_key='gpu',attempt=attempt),


use rule filter_genes from preprocessing as batch_analysis_filter_genes with:
    input:
        zarr=lambda wildcards: get_file(wildcards, 'normalize'),
    output:
        zarr=directory(mcfg.out_dir / 'prepare' / paramspace.wildcard_pattern / 'filter_genes.zarr'),
    resources:
        partition=lambda w, attempt: mcfg.get_resource(profile='gpu',resource_key='partition',attempt=attempt),
        qos=lambda w, attempt: mcfg.get_resource(profile='gpu',resource_key='qos',attempt=attempt),
        mem_mb=lambda w, attempt: mcfg.get_resource(profile='gpu',resource_key='mem_mb',attempt=attempt),
        gpu=lambda w, attempt: mcfg.get_resource(profile='gpu',resource_key='gpu',attempt=attempt),


use rule highly_variable_genes from preprocessing as batch_analysis_highly_variable_genes with:
    input:
        zarr=rules.batch_analysis_filter_genes.output.zarr,
    output:
        zarr=directory(mcfg.out_dir / 'prepare' / paramspace.wildcard_pattern / 'highly_variable_genes.zarr'),
    params:
        args=lambda wildcards: mcfg.get_from_parameters(wildcards, 'highly_variable_genes', default={}),
    resources:
        partition=lambda w, attempt: mcfg.get_resource(profile='gpu',resource_key='partition',attempt=attempt),
        qos=lambda w, attempt: mcfg.get_resource(profile='gpu',resource_key='qos',attempt=attempt),
        mem_mb=lambda w, attempt: mcfg.get_resource(profile='gpu',resource_key='mem_mb',attempt=attempt),
        gpu=lambda w, attempt: mcfg.get_resource(profile='gpu',resource_key='gpu',attempt=attempt),


use rule pca from preprocessing as batch_analysis_pca with:
    input:
        zarr=lambda wildcards: get_file(wildcards, 'highly_variable_genes'),
    output:
        zarr=directory(mcfg.out_dir / 'pca' / f'{paramspace.wildcard_pattern}.zarr'),
    params:
        args=lambda wildcards: mcfg.get_from_parameters(wildcards, 'pca', default={}),
    resources:
        partition=lambda w, attempt: mcfg.get_resource(profile='gpu',resource_key='partition',attempt=attempt),
        qos=lambda w, attempt: mcfg.get_resource(profile='gpu',resource_key='qos',attempt=attempt),
        mem_mb=lambda w, attempt: mcfg.get_resource(profile='gpu',resource_key='mem_mb',attempt=attempt),
        gpu=lambda w, attempt: mcfg.get_resource(profile='gpu',resource_key='gpu',attempt=attempt),
