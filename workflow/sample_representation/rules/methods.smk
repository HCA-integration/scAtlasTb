rule run_method:
    message:
       """
       Sample representation: Run method={wildcards.method} on dataset={wildcards.dataset} file_id={wildcards.file_id}
       input: {input}
       output: {output}
       wildcards: {wildcards}
       resources: gpu={resources.gpu} mem_mb={resources.mem_mb} partition={resources.partition} qos={resources.qos}
       """
    input:
        zarr=lambda wildcards: mcfg.get_input_file(**wildcards),
        prepare=rules.prepare.output.zarr,
        script=lambda wildcards: workflow.source_path(mcfg.get_from_parameters(wildcards, 'script'))
    output:
        zarr=directory(mcfg.out_dir / 'run_method' / f'{paramspace.wildcard_pattern}.zarr'),
    params:
        sample_key=lambda wildcards: mcfg.get_from_parameters(wildcards, 'sample_key'),
        cell_type_key=lambda wildcards: mcfg.get_from_parameters(wildcards, 'cell_type_key'),
        use_rep=lambda wildcards: mcfg.get_from_parameters(wildcards, 'use_rep'),
        var_mask=lambda wildcards: mcfg.get_from_parameters(wildcards, 'var_mask'),
        hyperparams=lambda wildcards: mcfg.get_from_parameters(wildcards, 'hyperparams_dict'),
        env=lambda wildcards: mcfg.get_from_parameters(wildcards, 'env'),
        script=lambda wildcards: mcfg.get_from_parameters(wildcards, 'script'),
        cran_url=config.get('cran_url', 'https://cloud.r-project.org'),
    conda:
        lambda wildcards, params: get_env(config, params.env)
    threads:
        lambda wildcards: max(1, mcfg.get_from_parameters(wildcards, 'threads', default=1)),
    resources:
        partition=lambda w, attempt: mcfg.get_resource(resource_key='partition', profile=mcfg.get_profile(w), attempt=attempt, attempt_to_cpu=2),
        qos=lambda w, attempt: mcfg.get_resource(resource_key='qos', profile=mcfg.get_profile(w), attempt=attempt, attempt_to_cpu=2),
        mem_mb=lambda w, attempt: mcfg.get_resource(resource_key='mem_mb', profile=mcfg.get_profile(w), attempt=attempt, attempt_to_cpu=2, factor=1),
        gpu=lambda w, attempt: mcfg.get_resource(resource_key='gpu', profile=mcfg.get_profile(w), attempt=attempt, attempt_to_cpu=2),
        time="2-00:00:00",
    script:
        '{params.script}'


rule run_method_all:
    input: mcfg.get_output_files(rules.run_method.output)
    localrule: True
