rule prepare:
    input:
        zarr=lambda wildcards: mcfg.get_input_file(**wildcards)
    output:
        zarr=directory(mcfg.out_dir / 'prepare' / 'dataset~{dataset}' / 'file_id~{file_id}.zarr'),
    params:
        sample_key=lambda wildcards: mcfg.get_from_parameters(wildcards, 'sample_key'),
        cell_type_key=lambda wildcards: mcfg.get_from_parameters(wildcards, 'cell_type_key'),
        norm_counts=lambda wildcards: mcfg.get_from_parameters(wildcards, 'norm_counts', default='X'),
        min_cells_per_sample=lambda wildcards: mcfg.get_from_parameters(wildcards, 'min_cells_per_sample', default=1),
        min_cells_per_cell_type=lambda wildcards: mcfg.get_from_parameters(wildcards, 'min_cells_per_cell_type', default=1),
        aggregate=lambda wildcards: mcfg.get_from_parameters(wildcards, 'aggregate', default='sum'),
    conda:
        get_env(config, 'sample_representation')
    resources:
        partition=mcfg.get_resource(resource_key='partition'),
        qos=mcfg.get_resource(profile='cpu',resource_key='qos'),
        mem_mb=lambda w, attempt: mcfg.get_resource(profile='cpu',resource_key='mem_mb', attempt=attempt),
        gpu=mcfg.get_resource(profile='cpu',resource_key='gpu'),
    script:
        '../scripts/prepare.py'


rule prepare_all:
    input:
        mcfg.get_output_files(rules.prepare.output),
    localrule: True
