rule plot_joint:
    input:
        zarr=rules.get_thresholds.output.zarr
    output:
        joint=directory(mcfg.image_dir / params.wildcard_pattern / 'joint_plots'),
    params:
        dataset=lambda wildcards: wildcards.file_id,
        hue=lambda wildcards: mcfg.get_from_parameters(wildcards, 'hue', default=[]),
        thresholds=lambda wildcards: mcfg.get_from_parameters(wildcards, 'thresholds', default={}),
        scautoqc_metrics=lambda wildcards: mcfg.get_from_parameters(wildcards, 'scautoqc_metrics', default=['n_counts', 'n_genes', 'percent_mito']),
        max_groups=lambda wildcards: mcfg.get_from_parameters(wildcards, 'plot_params', default={}).get('max_groups', 100),
    threads:
        lambda wildcards: max(1, min(10, len(mcfg.get_from_parameters(wildcards, 'hue', default=[]))))
    conda:
        get_env(config, 'scanpy')
    script:
        '../scripts/plot_joint.py'


rule plot_removed:
    input:
        zarr=rules.get_thresholds.output.zarr
    output:
        plots=directory(mcfg.image_dir / params.wildcard_pattern / 'removed'),
    params:
        dataset=lambda wildcards: wildcards.file_id,
        hue=lambda wildcards: mcfg.get_from_parameters(wildcards, 'hue', default=[]),
        thresholds=lambda wildcards: mcfg.get_from_parameters(wildcards, 'thresholds', default={}),
        scautoqc_metrics=lambda wildcards: mcfg.get_from_parameters(wildcards, 'scautoqc_metrics', default=['n_counts', 'n_genes', 'percent_mito']),
        max_groups=lambda wildcards: mcfg.get_from_parameters(wildcards, 'plot_params', default={}).get('max_groups', 100),
    threads:
        lambda wildcards: max(1, min(10, len(mcfg.get_from_parameters(wildcards, 'hue', default=[]))))
    conda:
        get_env(config, 'scanpy')
    script:
        '../scripts/plot_removed.py'


rule joint_plots:
    input: mcfg.get_output_files(rules.plot_joint.output)
    localrule: True


rule removed_plots:
    input: mcfg.get_output_files(rules.plot_removed.output)
    localrule: True


rule plots_all:
    input: rules.joint_plots.input, rules.removed_plots.input
    localrule: True