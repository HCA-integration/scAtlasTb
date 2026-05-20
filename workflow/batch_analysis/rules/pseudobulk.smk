use rule pca from preprocessing as batch_analysis_pb_pca with:
    input:
        zarr=rules.batch_analysis_prepare.output.bulks,
    output:
        zarr=directory(mcfg.out_dir / 'prepare' / paramspace.wildcard_pattern / 'pseudobulk_pca.zarr'),
    params:
        args=lambda wildcards: mcfg.get_from_parameters(wildcards, 'pca', default={}),
    resources:
        partition=lambda w, attempt: mcfg.get_resource(profile='cpu',resource_key='partition',attempt=attempt),
        qos=lambda w, attempt: mcfg.get_resource(profile='cpu',resource_key='qos',attempt=attempt),
        mem_mb=lambda w, attempt: mcfg.get_resource(profile='cpu',resource_key='mem_mb',attempt=attempt),
        gpu=lambda w, attempt: mcfg.get_resource(profile='cpu',resource_key='gpu',attempt=attempt),


use rule plots from preprocessing as batch_analysis_pb_pca_plot with:
    input:
        zarr=rules.batch_analysis_pb_pca.output.zarr,
    output:
        plots=directory(mcfg.image_dir / paramspace.wildcard_pattern / 'pca_plots'),
    params:
        basis='X_pca',
        color=lambda wildcards: mcfg.get_from_parameters(wildcards, 'covariates', default=[]),
        plot_centroids=lambda w: mcfg.get_from_parameters(w, 'pca_plot', default={}).get('plot_centroids', []),
        gene_chunk_size=lambda w: mcfg.get_from_parameters(w, 'pca_plot', default={}).get('plot_gene_chunk_size', 12),
        size=lambda w: mcfg.get_from_parameters(w, 'pca_plot', default={}).get('size', 50),
    resources:
        partition=lambda w, attempt: mcfg.get_resource(profile='cpu',resource_key='partition',attempt=attempt),
        qos=lambda w, attempt: mcfg.get_resource(profile='cpu',resource_key='qos',attempt=attempt),
        mem_mb=lambda w, attempt: mcfg.get_resource(profile='cpu',resource_key='mem_mb',attempt=attempt),
        gpu=lambda w, attempt: mcfg.get_resource(profile='cpu',resource_key='gpu',attempt=attempt),


rule pseudobulk_all:
    input:
        mcfg.get_output_files(rules.batch_analysis_prepare.output),
        mcfg.get_output_files(rules.batch_analysis_pb_pca_plot.output),
    localrule: True
