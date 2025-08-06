import json


class PathEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, Path):
            return str(obj)
        return super().default(obj)


rule save_config:
    output:
        json='.snakemake/{images}_{target}/config.json'
    localrule: True
    run:
        with open(output.json, 'w') as f:
            json.dump(config, f, cls=PathEncoder, indent=2)


rule rulegraph:
    input: rules.save_config.output.json
    output: '{images}/rule_graphs/{target}.png'
    localrule: True
    shell:
        """
        snakemake {wildcards.target} --configfile {input} --rulegraph | \
            sed -ne '/digraph snakemake_dag/,/}}/p' | \
            dot -Tpng -Grankdir=TB > {output}
        """


rule dag:
    input: rules.save_config.output.json
    output: '{images}/job_graphs/{target}.png'
    localrule: True
    shell:
        """
        snakemake {wildcards.target} --configfile {input} --dag | \
            sed -ne '/digraph snakemake_dag/,/}}/p' | \
            dot -Tpng -Grankdir=LR > {output}
        """


rule dependency_graph:
    input:
        rules.rulegraph.output,
        rules.dag.output,
    localrule: True
