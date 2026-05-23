Data exploration and semi-automated quality control
###################################################

Exploring your data properly is an important step for deciding which cells to include to your atlas.
Based on insights gathered from analysing the data pre integration, you can determine the quality of each input dataset, as well as the data completeness of metadata and features.
If the quality of a dataset is not sufficient for the atlas you are planning to build, you might want to consider excluding that dataset from the study.
This also applies to cases where metadata is missing and/or not obtainable, or whenever important genes are missing from the datasets.


Install dependencies
********************

For this workflow, make sure you install the following environments:

* `envs/scanpy.yaml`
* `envs/qc.yaml` for `qc` and `doublets` modules


Configuring a basic QC workflow
*******************************

Cell-level quality control can tell you about the quality of the data due to experimental conditions.
Even for published studies, where you would usually obtain already filtered data, you would want to investigate the QC thresholds that the authors have chosen.
In the following is an example on test data for setting up a first iteration of QC plots.
Start creating a config file called `configs/qc/qc_config.yaml` (see example under `configs`):

.. literalinclude:: ../../configs/qc/qc_config.yaml
   :language: yaml
   :caption: configs/qc/qc_config.yaml


Call the pipeline with either your runner script (e. g. called `configs/qc/run.sh`)

.. code-block:: bash

    bash configs/qc/run.sh qc_all -nq

You should get the following dry-run output:

.. code-block::

    Using profile .profiles/local for setting default command line arguments.
    WARNING: Duplicated columns: {'metric': ['methods', 'metrics']}
    Building DAG of jobs...
    Job stats:
    job                    count
    -------------------  -------
    qc_all                     1
    qc_autoqc                  8
    qc_get_thresholds          8
    qc_merge_thresholds        1
    qc_plot_joint              8
    qc_plot_removed            8
    qc_plot_summary            1
    total                     35

Execute the QC workflow as follows:

.. code-block:: bash

    bash configs/qc/run.sh qc_all -c2


This will estimate QC thresholds using the `sctk <https://teichlab.github.io/sctk/notebooks/automatic_qc.html>`_'s `scAutoQC` function.

Output
******

Check the outputs under `images/qc`.

.. code-block:: bash

      $ ls images/qc/dataset\~example_qc_analysis/

      file_id~test--split_data_value=empty
      file_id~test1--split_data_value=1
      file_id~test1--split_data_value=2
      file_id~test1--split_data_value=3
      file_id~test1--split_data_value=4
      file_id~test2--split_data_value=1
      file_id~test2--split_data_value=2
      file_id~test2--split_data_value=3
      qc_stats.tsv
      summary
      thresholds.tsv

There are per dataset and per file outputs.
The `thresholds.tsv` contains the estimated thresholds for both input files `test1` and `test2`, while `qc_stats.tsv` contains the overall statistics of cells that passed or failed the QC thresholds.


.. literalinclude:: ../_static/images/qc/dataset~example_qc_analysis/qc_stats.tsv
   :language: tsv
   :caption: qc_stats.tsv

.. literalinclude:: ../_static/images/qc/dataset~example_qc_analysis/thresholds.tsv
   :language: tsv
   :caption: thresholds.tsv


There are also separate outputs per original input file and split, which contains per `file_id` threshold files as well as `joint_plots` and `removed` directories.

Joint plots
===========

.. code-block:: bash

    $ ls images/qc/dataset~example_qc_analysis/file_id~test1--split_data_value=1/

    joint_plots  qc_stats.tsv  removed  thresholds.tsv


.. code-block:: bash

    $ ls images/qc/dataset\~example_qc_analysis/file_id\~test1--split_data_value=1/joint_plots

    'hue=batch.svg'  'hue=bulk_labels.svg'  'hue=phase.svg'  'hue=qc_status.svg'


Each scatter plot output is now a single SVG file containing all configured QC metric pairs, with linear-scale and log-scale panels arranged side by side.
The workflow writes one scatter plot per configured hue and always adds `hue=qc_status.svg`.
If no valid hue columns are available, the workflow writes `main.svg` instead of hue-specific scatter plots.
The density plots are also combined into a single `density.svg` file rather than separate files per metric pair.


Example outputs are shown below:

.. image:: ../_static/images/qc/dataset~example_qc_analysis/file_id~test1--split_data_value=1/joint_plots/hue=bulk_labels.svg
   :alt: bulk_labels_joint_plot

Removed plots
=============

The `removed` directory contains the QC summary plots for the same `file_id`:

.. code-block:: bash

    $ ls images/qc/dataset\~example_qc_analysis/file_id\~test1--split_data_value=1/removed

    'by=batch.svg'  'by=bulk_labels.svg'  'by=phase.svg'   cells_passed_all.svg   per_metric_violin.svg

The `by={hue}.svg` plots summarise removed cells per annotation, `cells_passed_all.svg` shows the overall QC status counts, and `per_metric_violin.svg` shows the per-metric distributions split by QC status.

Example removed-plot outputs:

.. image:: ../_static/images/qc/dataset~example_qc_analysis/file_id~test1--split_data_value=1/removed/by=bulk_labels.svg
   :alt: by_bulk_labels_removed_plot

.. image:: ../_static/images/qc/dataset~example_qc_analysis/file_id~test1--split_data_value=1/removed/cells_passed_all.svg
   :alt: cells_passed_all_removed_plot

.. image:: ../_static/images/qc/dataset~example_qc_analysis/file_id~test1--split_data_value=1/removed/per_metric_violin.svg
   :alt: per_metric_violin_removed_plot


Summary plots
=============

The dataset-level `summary` directory aggregates QC behaviour across all files and splits.
It contains removed-cell heatmaps and ridge plots for each configured QC metric.

.. code-block:: bash

   $ ls images/qc/dataset\~example_qc_analysis/summary

   n_removed_heatmap.svg  removed_frac_heatmap.svg  ridge_plots

.. code-block:: bash

   $ ls images/qc/dataset\~example_qc_analysis/summary/ridge_plots

   log1p_n_counts.svg  log1p_n_genes.svg  n_counts.svg  n_genes.svg  percent_mito.svg

The ridge plots show the metric distributions per input file and split, with the configured thresholds overlaid.
The heatmaps summarise absolute and relative cell removal across the dataset.

.. image:: ../_static/images/qc/dataset~example_qc_analysis/summary/ridge_plots/percent_mito.svg
   :alt: percent_mito_ridge_plot

.. image:: ../_static/images/qc/dataset~example_qc_analysis/summary/removed_frac_heatmap.svg
   :alt: removed_fraction_heatmap


Full Configuration
******************

You can find the complete configuration file and runner script under `configs/qc/`.
Here's the final workflow configuration:

.. literalinclude:: ../../configs/qc/qc_config.yaml
   :language: yaml
   :caption: configs/qc/qc_config.yaml

.. literalinclude:: ../../configs/qc/run.sh
   :language: bash
   :caption: configs/qc/run.sh

