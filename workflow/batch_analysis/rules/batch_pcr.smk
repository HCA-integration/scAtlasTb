checkpoint determine_covariates:
    input:
        zarr=lambda wildcards: get_file(wildcards, 'pca'),
    output:
        covariate_setup=directory(mcfg.out_dir / paramspace.wildcard_pattern / 'batch_pcr' / 'covariate_setup'),
    params:
        covariates=lambda wildcards: mcfg.get_from_parameters(wildcards, 'covariates', default=[]),
        permute_covariates=lambda wildcards: mcfg.get_from_parameters(wildcards, 'permute_covariates', default=None),
        sample_key=lambda wildcards: mcfg.get_from_parameters(wildcards, 'sample', check_query_keys=False),
        n_permute=lambda wildcards: mcfg.get_from_parameters(wildcards, 'n_permutations', default=10),
    conda:
        get_env(config, 'scanpy')
    localrule: True
    script:
        '../scripts/determine_covariates.py'


def get_checkpoint_output(wildcards):
    return f'{checkpoints.determine_covariates.get(**wildcards).output[0]}/{{covariate}}.yaml'


def get_from_checkpoint(wildcards, pattern=None):
    checkpoint_output = get_checkpoint_output(wildcards)
    if pattern is None:
        pattern = checkpoint_output
    return expand(
        pattern,
        covariate=glob_wildcards(checkpoint_output).covariate,
        allow_missing=True
    )


rule batch_pcr:
    input:
        zarr=lambda wildcards: get_file(wildcards, 'pca'),
        setup=get_checkpoint_output,
    output:
        tsv=mcfg.out_dir / paramspace.wildcard_pattern / 'batch_pcr' / '{covariate}.tsv',
    params:
        sample_key=lambda wildcards: mcfg.get_from_parameters(wildcards, 'sample', check_query_keys=False),
    conda:
        get_env(config, 'scib')
    threads:
        lambda w: max(
            1, min(
                mcfg.get_from_parameters(w, 'max_threads', check_query_keys=False, default=50),
                mcfg.get_from_parameters(w, 'n_permutations', check_query_keys=False)
            )
        )
    resources:
        partition=mcfg.get_resource(profile='cpu',resource_key='partition'),
        qos=mcfg.get_resource(profile='cpu',resource_key='qos'),
        mem_mb=lambda w, attempt: mcfg.get_resource(profile='cpu',resource_key='mem_mb', attempt=attempt),
    script:
        '../scripts/batch_pcr.py'


rule collect:
    input:
        tsv=lambda wildcards: get_from_checkpoint(wildcards, rules.batch_pcr.output.tsv),
        setup=lambda wildcards: checkpoints.determine_covariates.get(**wildcards).output[0],
    output:
        tsv=mcfg.out_dir / paramspace.wildcard_pattern / 'batch_pcr.tsv',
    run:
        dfs = [pd.read_table(file) for file in input.tsv]
        if len(dfs) == 0:
            df = pd.DataFrame(columns=['covariate', 'pcr', 'permuted', 'n_covariates', 'non_perm_z_score'])
        else:
            df = pd.concat(dfs, ignore_index=True)
        df.to_csv(output.tsv, sep='\t', index=False)


rule plot:
    input:
        tsv=rules.collect.output.tsv,
    output:
        barplot=mcfg.image_dir / paramspace.wildcard_pattern / 'batch_pcr_bar.png',
        violinplot=mcfg.image_dir / paramspace.wildcard_pattern / 'batch_pcr_violin.png',
    params:
        n_permute=lambda wildcards: mcfg.get_from_parameters(wildcards, 'n_permutations', default=10),
    conda:
        get_env(config, 'scanpy')
    script:
        '../scripts/plot.py'
