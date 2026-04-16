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
        raw_counts=lambda wildcards: mcfg.get_from_parameters(wildcards, 'raw_counts'),
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
        zarr=directory(mcfg.out_dir / 'prepare' / paramspace.wildcard_pattern / 'pca.zarr'),
    params:
        args=lambda wildcards: mcfg.get_from_parameters(wildcards, 'pca', default={}),
    resources:
        partition=lambda w, attempt: mcfg.get_resource(profile='gpu',resource_key='partition',attempt=attempt),
        qos=lambda w, attempt: mcfg.get_resource(profile='gpu',resource_key='qos',attempt=attempt),
        mem_mb=lambda w, attempt: mcfg.get_resource(profile='gpu',resource_key='mem_mb',attempt=attempt),
        gpu=lambda w, attempt: mcfg.get_resource(profile='gpu',resource_key='gpu',attempt=attempt),


use rule pseudobulk from sample_representation as batch_analysis_prepare with:
    input:
        zarr=lambda wildcards: mcfg.get_input_file(**wildcards),
    output:
        zarr=directory(mcfg.out_dir / 'prepare' / paramspace.wildcard_pattern / 'prepare.zarr'),
        bulks=directory(mcfg.out_dir / 'prepare' / paramspace.wildcard_pattern / 'pseudobulks.zarr'),
    params:
        sample_key=lambda wildcards: mcfg.get_from_parameters(wildcards, 'sample'),
        layer=lambda wildcards: mcfg.get_from_parameters(wildcards, 'raw_counts', default='X'),
        # aggregate=lambda wildcards: mcfg.get_from_parameters(wildcards, 'aggregate', default='sum'),
        # dask=lambda wildcards: mcfg.get_from_parameters(wildcards, 'dask', default=None),
    resources:
        partition=mcfg.get_resource(resource_key='partition'),
        qos=mcfg.get_resource(profile='cpu',resource_key='qos'),
        mem_mb=lambda w, attempt: mcfg.get_resource(profile='cpu',resource_key='mem_mb', attempt=attempt),
        gpu=mcfg.get_resource(profile='cpu',resource_key='gpu'),