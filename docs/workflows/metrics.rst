Metrics
#######

Given multiple representations of your data, you can evaluate the biological conservation and the batch removal to evaluate the quality of your representation with respect to batch correction.

Install dependencies
********************

For this workflow, make sure you install the following environments:

* `envs/scanpy.yaml`
* `envs/scib.yaml` for `scib` metrics https://github.com/theislab/scib
* `envs/scib_accel.yaml` for intel-accelerated `scib` metrics (TODO: deprecate)
* `envs/scib_metrics.yaml` for `scib_metrics` metrics https://github.com/YosefLab/scib-metrics
* `envs/funkyheatmap.yaml` for `funkyheatmap` plots of the metrics https://github.com/funkyheatmap/funkyheatmap

Configuring a basic metrics workflow
************************************

Given a set of files that you want to compare, configure your workflow as follows:

.. code-block:: yaml

    output_dir: test/out
    images: test/images

    use_gpu: false # only set to True if you are working with GPUs

    DATASETS:
      test_metrics:
        input:
          metrics:
            model1--model1_param=val1: test/input/pbmc68k.h5ad
            model2--model2_param=val2: test/input/pbmc68k.h5ad
            unintegrated: test/input/pbmc68k.h5ad
        metrics:
          label: bulk_labels
          batch: louvain
          overwrite_file_id: true
          metrics:
            - nmi
            - ari
            - asw_label
            - asw_batch
            - cell_cycle
            - clisi
            - ilisi
            - graph_connectivity
            - isolated_label_asw
            - isolated_label_f1
            - pcr_comparison
            - pcr_comparison
            - kbet_pg

Call the pipeline with either your runner script (e. g. called `configs/qc/run.sh`)

.. code-block:: bash

    bash configs/metrics/run.sh qc_all -nq

You should get the following dry-run output:

.. code-block::

    Config file configs/outputs.yaml is extended by additional config specified via the command line.
    Config file configs/load_data/config.yaml is extended by additional config specified via the command line.
    Config file configs/exploration/config.yaml is extended by additional config specified via the command line.
    WARNING: Duplicated columns: {'metric': ['methods', 'metrics']}
    Building DAG of jobs...
    Job stats:
    job                                 count
    --------------------------------  -------
    metrics_all                             1
    metrics_barplot                         3
    metrics_barplot_per_dataset             3
    metrics_barplot_per_file                9
    metrics_cluster                        30
    metrics_cluster_collect                 3
    metrics_collect                         3
    metrics_funkyheatmap                    1
    metrics_funkyheatmap_per_dataset        1
    metrics_merge                           1
    metrics_merge_per_batch                 1
    metrics_merge_per_dataset               1
    metrics_merge_per_file                  3
    metrics_merge_per_label                 1
    metrics_prepare                         3
    metrics_run                            36
    total                                 100

Execute the workflow as follows:

.. code-block:: bash

    bash configs/qc/run.sh metrics_all -c5

Output
******

Check the outputs under `images/metrics`.
Per dataset plots are under `images/metrics/per_dataset`

.. figure:: ../_static/images/metrics/per_dataset/test_metrics/funky_heatmap.png
    :alt: funkyheatmap
    :width: 100%
    :align: center
    
    Metrics FunkyHeatmap from
    images/metrics/per_dataset/test_metrics/funky_heatmap.pdf

.. figure:: ../_static/images/metrics/per_dataset/test_metrics/score-barplot.png
    :alt: barplot
    :width: 80%
    :align: center
    
    Metrics values in a barplot from
    images/metrics/per_dataset/test_metrics/score-barplot.png


.. figure:: ../_static/images/metrics/per_dataset/test_metrics/s-barplot.png
    :alt: barplot
    :width: 80%
    :align: center
    
    Compute duration in seconds from
    images/metrics/per_dataset/test_metrics/s-barplot.png


TL;DR Full Configuration
************************

You can find the complete configuration file and runner script under `configs/metrics/`.
Here's the final workflow configuration:

.. literalinclude:: ../../configs/metrics/example_workflow.yaml
   :language: yaml
   :caption: configs/metrics/example_workflow.yaml

.. literalinclude:: ../../configs/metrics/run.sh
   :language: yaml
   :caption: configs/metrics/run.sh