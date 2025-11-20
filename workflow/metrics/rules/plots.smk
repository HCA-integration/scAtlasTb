from utils.wildcards import wildcards_to_str

module plots:
    snakefile: "../../common/rules/plots.smk"
    config: config


# barplots

use rule barplot from plots as metrics_barplot with:
    input:
        tsv=rules.merge.output.tsv
    output:
        png=mcfg.image_dir / 'all' /'{metric}-barplot.png'
    params:
        metric=lambda wildcards: wildcards.metric,
        category='metric',
        hue='file_id',
        facet_row='dataset',
        #facet_col='output_type',
        title='Metrics computation',
        description=wildcards_to_str,
        dodge=True,
    group:
        'metrics_merge'


use rule barplot from plots as metrics_barplot_per_dataset with:
    input:
        tsv=rules.merge_per_dataset.output.tsv
    output:
        png=mcfg.image_dir / 'dataset~{dataset}' / '{metric}-barplot.png'
    params:
        metric=lambda wildcards: wildcards.metric,
        category='metric',
        hue='file_id',
        #facet_col='output_type',
        title='Metrics computation',
        description=wildcards_to_str,
        dodge=True,
    group:
        'metrics_merge'


use rule barplot from plots as metrics_barplot_per_file with:
    input:
        tsv=rules.merge_per_file.output.tsv
    output:
        png=mcfg.image_dir / paramspace.wildcard_pattern / '{metric}-barplot.png'
    params:
        metric=lambda wildcards: wildcards.metric,
        category='metric',
        hue='dataset',
        #facet_col='output_type',
        title='Metrics computation',
        description=wildcards_to_str,
        dodge=True,
    group:
        'metrics_merge'


# swarm plots

use rule swarmplot from plots as metrics_swarmplot with:
    input:
        tsv=rules.merge.output.tsv
    output:
        png=mcfg.image_dir / 'all' / '{metric}-swarmplot.png'
    params:
        metric=lambda wildcards: wildcards.metric,
        category='metric',
        hue='file_id',
        facet_col='dataset',
        facet_row='output_type',
        title='Metrics',
        description=wildcards_to_str,
    group:
        'metrics_merge'


use rule swarmplot from plots as metrics_swarmplot_per_dataset with:
    input:
        tsv=rules.merge_per_dataset.output.tsv
    output:
        png=mcfg.image_dir / 'dataset~{dataset}' / '{metric}-swarmplot.png'
    params:
        metric=lambda wildcards: wildcards.metric,
        category='metric',
        hue='file_id',
        facet_row='output_type',
        title='Metrics',
        description=wildcards_to_str,
        dodge=True,
    group:
        'metrics_merge'


use rule swarmplot from plots as metrics_swarmplot_per_file with:
    input:
        tsv=rules.merge_per_file.output.tsv
    output:
        png=mcfg.image_dir / 'per_file' / '{file_id}' / '{metric}-swarmplot.png'
    params:
        metric=lambda wildcards: wildcards.metric,
        category='metric',
        hue='dataset',
        facet_row='output_type',
        title='Metrics',
        description=wildcards_to_str,
        dodge=True,


rule compare_metrics:
    input:
        tsv=rules.merge.output.tsv
    output:
        time=mcfg.image_dir / 'comparison_time.png',
        score=mcfg.image_dir / 'comparison_score.png',
    conda:
        get_env(config, 'plots')
    localrule: True
    script:
        '../scripts/plots/comparison.py'


# Funkyheatmap

rule funkyheatmap:
    input:
        tsv=rules.merge.output.tsv,
        extra_columns=rules.merge.output.extra_columns,
        r_utils=workflow.source_path('../scripts/plots/r_utils.R'),
    output:
        pdf=mcfg.image_dir / 'all' / 'funky_heatmap.pdf',
        tsv=mcfg.image_dir / 'all' / 'funky_heatmap.tsv'
    params:
        id_vars=['dataset', 'output_type', 'batch', 'label'], # TODO: 'hyperparams'
        variable_var='metric',
        value_var='score',
        weight_batch=lambda wildcards: mcfg.get_from_parameters(wildcards, 'weight_batch', default=0.4),
        n_top=50,
        cran_url=config.get('cran_url', 'https://cloud.r-project.org'), #'https://ftp.fau.de/cran/'
    conda:
        get_env(config, 'funkyheatmap')  # TODO: use post-deployment script https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#providing-post-deployment-scripts
    singularity:
        'docker://ghcr.io/dynverse/funky_heatmap:latest'
    group:
        'metrics_merge'
    script:
        '../scripts/plots/funkyheatmap.R'


use rule funkyheatmap as funkyheatmap_per_dataset with:
    input:
        tsv=rules.merge_per_dataset.output.tsv,
        extra_columns=rules.merge_per_dataset.output.extra_columns,
        r_utils=workflow.source_path('../scripts/plots/r_utils.R'),
    output:
        pdf=mcfg.image_dir / 'dataset~{dataset}' / 'funky_heatmap.pdf',
        tsv=mcfg.image_dir / 'dataset~{dataset}' / 'funky_heatmap.tsv',
    params:
        id_vars=['dataset', 'output_type', 'batch', 'label'],
        variable_var='metric',
        value_var='score',
        weight_batch=lambda wildcards: mcfg.get_from_parameters(wildcards, 'weight_batch', default=0.4),
        n_top=50,
        cran_url=config.get('cran_url', 'https://cloud.r-project.org'),
    group:
        'metrics_merge'


use rule funkyheatmap as funkyheatmap_per_batch with:
    input:
        tsv=rules.merge_per_batch.output.tsv,
        extra_columns=rules.merge_per_batch.output.extra_columns,
        r_utils=workflow.source_path('../scripts/plots/r_utils.R'),
    output:
        pdf=mcfg.image_dir /  'dataset~{dataset}' / 'batch~{batch}' / 'funky_heatmap.pdf',
        tsv=mcfg.image_dir /  'dataset~{dataset}' / 'batch~{batch}' / 'funky_heatmap.tsv',
    params:
        id_vars=['dataset', 'output_type', 'batch', 'label'],
        variable_var='metric',
        value_var='score',
        weight_batch=lambda wildcards: mcfg.get_from_parameters(wildcards, 'weight_batch', default=0.4),
        n_top=50,
        cran_url=config.get('cran_url', 'https://cloud.r-project.org'),
    group:
        'metrics_merge'


use rule funkyheatmap as funkyheatmap_per_label with:
    input:
        tsv=rules.merge_per_label.output.tsv,
        extra_columns=rules.merge_per_label.output.extra_columns,
        r_utils=workflow.source_path('../scripts/plots/r_utils.R'),
    output:
        pdf=mcfg.image_dir / 'dataset~{dataset}' / 'label~{label}' / 'funky_heatmap.pdf',
        tsv=mcfg.image_dir / 'dataset~{dataset}' / 'label~{label}' / 'funky_heatmap.tsv',
    params:
        id_vars=['dataset', 'output_type', 'batch', 'label'],
        variable_var='metric',
        value_var='score',
        weight_batch=lambda wildcards: mcfg.get_from_parameters(wildcards, 'weight_batch', default=0.4),
        n_top=50,
        cran_url=config.get('cran_url', 'https://cloud.r-project.org'),
    group:
        'metrics_merge'


use rule funkyheatmap as funkyheatmap_per_file with:
    input:
        tsv=rules.merge_per_file.output.tsv,
        extra_columns=rules.merge_per_file.output.extra_columns,
        r_utils=workflow.source_path('../scripts/plots/r_utils.R'),
    output:
        pdf=mcfg.image_dir / 'per_file' / '{file_id}' / 'funky_heatmap.pdf',
        tsv=mcfg.image_dir / 'per_file' / '{file_id}' / 'funky_heatmap.tsv',
    params:
        id_vars=['dataset', 'output_type', 'batch', 'label'],
        variable_var='metric',
        value_var='score',
        weight_batch=lambda wildcards: mcfg.get_from_parameters(wildcards, 'weight_batch', default=0.4),
        n_top=50,
        cran_url=config.get('cran_url', 'https://cloud.r-project.org'),
    group:
        'metrics_merge'


rule funkyheatmap_standalone:
    input:
        tsv=rules.merge.output.tsv
    output:
        png=mcfg.image_dir / 'funky_heatmap.png'
    localrule: True
    shell:
        """
        wget https://github.com/dynverse/funkyheatmap/releases/latest/download/executable.zip
        unzip -o executable.zip -d funky_heatmap
        ./funky_heatmap/funky_heatmap --data {input.tsv} --output {output.png}
        """


rule plots_all:
    input:
        # funky heatmap
        rules.funkyheatmap.output,
        mcfg.get_output_files(rules.funkyheatmap_per_dataset.output),
        mcfg.get_output_files(rules.funkyheatmap_per_batch.output, wildcard_names=['dataset', 'batch']),
        mcfg.get_output_files(rules.funkyheatmap_per_label.output, wildcard_names=['dataset', 'label']),
        # # barplot
        expand(
            rules.metrics_barplot.output,
            metric=['s', 'max_uss', 'score']
        ),
        expand(
            rules.metrics_barplot_per_dataset.output,
            metric=['s', 'max_uss', 'score'],
            **mcfg.get_wildcards(wildcard_names=['dataset'])
        ),
        # expand(
        #     rules.metrics_barplot_per_file.output,
        #     metric=['s', 'max_uss', 'score'],
        #     **mcfg.get_wildcards(wildcard_names=['dataset', 'file_id'])
        # ),
        # # swarmplot
        # expand(rules.metrics_swarmplot.output,metric='score'),
        # expand(rules.metrics_swarmplot_per_dataset.output,metric='score',**mcfg.get_wildcards(wildcard_names=['dataset'])),
        # expand(rules.metrics_swarmplot_per_file.output,metric='score',**mcfg.get_wildcards(wildcard_names=['file_id'])),
        # # implementation comparison
        # rules.compare_metrics.output,
    localrule: True
