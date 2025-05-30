checkpoint split_batches:
    input:
        zarr=lambda wildcards: mcfg.get_input_file(**wildcards)
    output:
        batches=directory(mcfg.out_dir / 'scatter' /  params.wildcard_pattern / '_batches'),
    params:
        layer=lambda wildcards: mcfg.get_from_parameters(wildcards, 'counts', default='X'),
        batch_key=lambda wildcards: mcfg.get_from_parameters(wildcards, 'batch_key'),
    conda:
        get_env(config, 'scanpy')
    script:
        '../scripts/split_batches.py'


def get_checkpoint_output(wildcards):
    return f'{checkpoints.split_batches.get(**wildcards).output[0]}/{{batch}}.txt'


def get_from_checkpoint(wildcards, pattern=None):
    checkpoint_output = get_checkpoint_output(wildcards)
    if pattern is None:
        pattern = checkpoint_output
    return expand(
        pattern,
        batch=glob_wildcards(checkpoint_output).batch,
        allow_missing=True
    )


def get_mem_mb(attempt, profile):
    mem_mb = mcfg.get_resource(profile=profile, resource_key='mem_mb', attempt=attempt)
    try:
        mem_mb = int(mem_mb)
    except ValueError:
        return mem_mb
    return min(100_000, int(mem_mb // 8))


rule scrublet:
    input:
        zarr=lambda wildcards: mcfg.get_input_file(**wildcards),
        batch=get_checkpoint_output,
    output:
        tsv=mcfg.out_dir / 'scatter' / params.wildcard_pattern / 'scrublet' / '{batch}.tsv',
    params:
        batch_key=lambda wildcards: mcfg.get_from_parameters(wildcards, 'batch_key', check_query_keys=False),
    conda:
        get_env(config, 'qc', gpu_env='rapids_singlecell')
    resources:
        partition=lambda w, attempt: mcfg.get_resource(profile='gpu',resource_key='partition', attempt=attempt),
        qos=lambda w, attempt: mcfg.get_resource(profile='gpu',resource_key='qos', attempt=attempt),
        gpu=lambda w, attempt: mcfg.get_resource(profile='gpu',resource_key='gpu', attempt=attempt),
        mem_mb=lambda w, attempt: get_mem_mb(attempt, profile='gpu'),
    script:
        '../scripts/scrublet.py'


rule doubletdetection:
    input:
        zarr=lambda wildcards: mcfg.get_input_file(**wildcards),
        batch=get_checkpoint_output
    output:
        tsv=mcfg.out_dir / 'scatter' / params.wildcard_pattern / 'doubletdetection' / '{batch}.tsv',
    params:
        batch_key=lambda wildcards: mcfg.get_from_parameters(wildcards, 'batch_key', check_query_keys=False),
    conda:
        get_env(config, 'qc')
    threads: 3
    resources:
        partition=mcfg.get_resource(profile='cpu',resource_key='partition'),
        qos=mcfg.get_resource(profile='cpu',resource_key='qos'),
        gpu=mcfg.get_resource(profile='cpu',resource_key='gpu'),
        mem_mb=lambda w, attempt: get_mem_mb(attempt, profile='cpu'),
    # shadow: "minimal"
    script:
        '../scripts/doubletdetection.py'


rule collect:
    input:
        zarr=lambda wildcards: mcfg.get_input_file(**wildcards),
        scrublet=lambda wildcards: get_from_checkpoint(wildcards, rules.scrublet.output.tsv),
        doubletdetection=lambda wildcards: get_from_checkpoint(wildcards, rules.doubletdetection.output.tsv),
    output:
        zarr=directory(mcfg.out_dir / f'{params.wildcard_pattern}.zarr'),
    localrule: True
    conda:
        get_env(config, 'scanpy')
    resources:
        mem_mb=mcfg.get_resource(profile='cpu',resource_key='mem_mb')
    script:
        '../scripts/collect.py'
