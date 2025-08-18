rule prepare:
    input:
        lambda wildcards: mcfg.get_input_file(**wildcards),
    output:
        zarr=directory(mcfg.out_dir / 'prepare' / paramspace.wildcard_pattern / 'prepare.zarr'),
    params:
        # labels=lambda wildcards: mcfg.get_from_parameters(wildcards, 'label', single_value=False),
        neighbor_args=lambda wildcards: mcfg.get_for_dataset(wildcards.dataset, ['preprocessing', 'neighbors'], default={}),
        unintegrated_layer=lambda wildcards: mcfg.get_from_parameters(wildcards, 'unintegrated', default='X'),
        corrected_layer=lambda wildcards: mcfg.get_from_parameters(wildcards, 'corrected', default='X'),
        var_mask=lambda wildcards: mcfg.get_from_parameters(wildcards, 'var_mask', default='highly_variable'),
    conda:
        get_env(config, 'scanpy', gpu_env='rapids_singlecell')
    resources:
        partition=lambda w: mcfg.get_resource(resource_key='partition', profile='gpu'),
        qos=lambda w: mcfg.get_resource(resource_key='qos', profile='gpu'),
        gpu=lambda w: mcfg.get_resource(resource_key='gpu', profile='gpu'),
        mem_mb=lambda w, attempt: mcfg.get_resource(resource_key='mem_mb', profile='gpu', attempt=attempt),
        time="1-00:00:00",
    script:
        '../scripts/prepare.py'


rule prepare_all:
    input:
        mcfg.get_output_files(rules.prepare.output)
    localrule: True


# Clustering for cluster-based metrics

def scale_mem_mb(wildcards, attempt, factor=3, min_mb=1_000, profile=None):
    # TODO: add to utils
    if profile is None:
        profile = mcfg.get_profile(wildcards)
    mem_mb = mcfg.get_resource(profile=profile, resource_key='mem_mb', attempt=attempt)

    try:
        mem_mb = int(mem_mb)
    except ValueError:
        return mem_mb
    
    mem_mb = int(mem_mb // factor)
    return max(min_mb, mem_mb)


use rule cluster from clustering as metrics_cluster with:
    input:
        zarr=rules.prepare.output.zarr,
    output:
        zarr=directory(mcfg.out_dir / 'prepare' / paramspace.wildcard_pattern / 'cluster_resolutions' / '{algorithm}--{resolution}--{level}.zarr'),
    params:
        clustering_args=lambda wildcards: mcfg.get_from_parameters(wildcards, 'clustering', default={}).get('kwargs', {}),
        overwrite=lambda wildcards: mcfg.get_from_parameters(wildcards, 'clustering', default={}).get('overwrite', True),
    conda:
        get_env(config, 'scanpy', gpu_env='rapids_singlecell')
    resources:
        partition=lambda w, attempt: mcfg.get_resource(profile='gpu',resource_key='partition',attempt=attempt),
        qos=lambda w, attempt: mcfg.get_resource(profile='gpu',resource_key='qos',attempt=attempt),
        # mem_mb=lambda w, attempt: mcfg.get_resource(profile='gpu',resource_key='mem_mb',attempt=attempt),
        mem_mb=lambda w, attempt: scale_mem_mb(w, attempt, factor=2, profile='gpu'),
        gpu=lambda w, attempt: mcfg.get_resource(profile='cpu',resource_key='gpu',attempt=attempt),


use rule merge from clustering as metrics_cluster_collect with:
    input:
        zarr=rules.prepare.output.zarr,
        cluster_anno=expand(
            rules.metrics_cluster.output.zarr,
            resolution=[2 * (x + 1) / 10 for x in range(10)],
            algorithm='leiden',
            level=1,
            allow_missing=True,
        ),
    output:
        zarr=directory(mcfg.out_dir / 'prepare' / paramspace.wildcard_pattern / 'clustered.zarr'),


rule cluster_all:
    input:
        mcfg.get_output_files(rules.metrics_cluster_collect.output)
    localrule: True


# Gene scoring for gene-based metrics

rule score_genes:
    input:
        zarr=rules.prepare.output.zarr,
    output:
        zarr=directory(mcfg.out_dir / 'prepare' / paramspace.wildcard_pattern / 'gene_scores.zarr'),
    params:
        unintegrated_layer=lambda wildcards: mcfg.get_from_parameters(wildcards, 'unintegrated', default='X'),
        raw_counts_layer=lambda wildcards: mcfg.get_from_parameters(wildcards, 'raw_counts', default=None),
        gene_sets=lambda wildcards: mcfg.get_gene_sets(wildcards.dataset),
        n_permutations=lambda wildcards: mcfg.get_from_parameters(wildcards, 'n_permutations', default=20),
        n_quantiles=lambda wildcards: mcfg.get_from_parameters(wildcards, 'n_quantiles', default=5),
    conda:
        get_env(config, 'scanpy', gpu_env='rapids_singlecell')
    resources:
        partition=lambda w, attempt: mcfg.get_resource(profile='gpu',resource_key='partition',attempt=attempt),
        qos=lambda w, attempt: mcfg.get_resource(profile='gpu',resource_key='qos',attempt=attempt),
        mem_mb=lambda w, attempt: mcfg.get_resource(profile='gpu',resource_key='mem_mb',attempt=attempt),
        gpu=lambda w, attempt: mcfg.get_resource(profile='gpu',resource_key='gpu',attempt=attempt),
    script:
        '../scripts/score_genes.py'


rule score_genes_all:
    input: mcfg.get_output_files(rules.score_genes.output)
    localrule: True
