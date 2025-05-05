def get_annotated_study(wildcards):
    if wildcards.study in dcp_studies:
        return dict(rules.add_dcp_metadata.output)
    return dict(rules.load_data_merge_study.output)


use rule filter from load_data_filter as load_data_filter_study with:
    message:
        """
        Filter the data for study {wildcards.study}
        input: {input}
        output: {output}
        remove_by_column: {params.remove_by_column}
        """
    input:
        unpack(get_annotated_study)
    output:
        zarr=directory(out_dir / 'filtered' / '{study}.zarr'),
    params:
        remove_by_column=lambda wildcards: config['filter_per_study'][wildcards.study].get('remove_by_column', {}),
        backed=False,
        dask=True,
        subset=False,
    resources:
        mem_mb=get_resource(config,profile='cpu',resource_key='mem_mb'),


rule filter_all:
    input: expand(rules.load_data_filter_study.output,**get_wildcards(dataset_df,['study']))
    localrule: True
