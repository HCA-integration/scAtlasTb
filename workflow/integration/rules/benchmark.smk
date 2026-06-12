from utils.wildcards import wildcards_to_str


rule benchmark_per_dataset:
    message: "[Integration] Collect benchmarks: dataset={wildcards.dataset}"
    input:
        benchmark=lambda wildcards: mcfg.get_output_files(
            rules.integration_run_method.benchmark,
            subset_dict=dict(wildcards)
        ),
    output:
        benchmark=out_dir / 'dataset~{dataset}' / 'integration.benchmark.tsv'
    params:
        wildcards=lambda wildcards: mcfg.get_wildcards(subset_dict=wildcards, as_df=True)
    localrule: True
    run:
        benchmark_df = pd.concat([pd.read_table(file) for file in input.benchmark]).reset_index(drop=True)
        benchmark_df = pd.concat([params.wildcards.reset_index(drop=True), benchmark_df],axis=1)
        print(benchmark_df)
        benchmark_df.to_csv(output.benchmark,sep='\t',index=False)


rule benchmark:
    message: "[Integration] Collect benchmarks: all datasets"
    input:
        benchmark=mcfg.get_output_files(rules.integration_run_method.benchmark),
    output:
        benchmark=out_dir / 'integration.benchmark.tsv'
    params:
        wildcards=mcfg.get_wildcards(as_df=True)
    localrule: True
    run:
        benchmark_df = pd.concat([pd.read_table(file) for file in input.benchmark]).reset_index(drop=True)
        benchmark_df = pd.concat([params.wildcards.reset_index(drop=True), benchmark_df],axis=1)
        print(benchmark_df)
        benchmark_df.to_csv(output.benchmark,sep='\t',index=False)


use rule barplot from plots as integration_barplot with:
    message: "[Integration] Benchmark barplot: metric={wildcards.metric}"
    input:
        tsv=rules.benchmark.output.benchmark
    output:
        png=image_dir / 'benchmark' / '{metric}.png'
    params:
        metric=lambda wildcards: wildcards.metric,
        category='method',
        #hue='hyperparams',
        facet_col='dataset',
        title='Integration runtime',
        description=wildcards_to_str,
        dodge=True,
    wildcard_constraints:
        metric='((?![/]).)*'


use rule barplot from plots as integration_barplot_per_dataset with:
    message: "[Integration] Benchmark barplot: dataset={wildcards.dataset} metric={wildcards.metric}"
    input:
        tsv=rules.benchmark_per_dataset.output.benchmark
    output:
        png=image_dir / 'benchmark' / 'dataset~{dataset}' / 'metric~{metric}.png'
    params:
        metric=lambda wildcards: wildcards.metric,
        category='method',
        #hue='hyperparams',
        title=lambda wildcards: f'Integration runtime for {wildcards.dataset}',
        description=wildcards_to_str,
        dodge=True,
    wildcard_constraints:
        metric='((?![/]).)*'


rule benchmark_all:
    input:
        expand(
            mcfg.get_output_files(
                rules.integration_barplot.output,
                allow_missing=True,
            ),
            metric=['s', 'max_uss', 'mean_load']
        ),
        expand(
            mcfg.get_output_files(
                rules.integration_barplot_per_dataset.output,
                allow_missing=True,
            ),
            metric=['s', 'max_uss', 'mean_load'],
        )
    localrule: True
