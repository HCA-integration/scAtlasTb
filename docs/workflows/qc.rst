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
* `envs/qc.yaml` for `qc` module
* `envs/doublets.yaml` for `doublets` module


Configuring a basic QC workflow
*******************************

Cell-level quality control can tell you about the quality of the data due to experimental conditions.
Even for published studies, where you would usually obtain already filtered data, you would want to investigate the QC thresholds that the authors have chosen.
In the following is an example on test data for setting up a first iteration of QC plots.
Start creating a config file called `configs/qc/qc_config.yaml` (see example under `configs`):

.. code-block:: yaml

    output_dir: data/out
    images: images
    user_gpu: false

    DATASETS:

      example_qc_analysis:
        input:
          qc:
            test1: data/pbmc68k.h5ad # adjust to custom file path
            test2: data/pbmc68k.h5ad # adjust to custom file path
          merge: qc
        qc:
          counts: X  # raw counts slot
          hue: # obs columns you want in your plot
            - bulk_labels
            - batch
            - phase

Call the pipeline with either your runner script (e. g. called `configs/qc/run.sh`)

.. code-block:: bash

    bash configs/qc/run.sh qc_all -nq

You should get the following dry-run output:

.. code-block::

    Config file configs/outputs.yaml is extended by additional config specified via the command line.
    Config file configs/load_data/config.yaml is extended by additional config specified via the command line.
    Config file configs/exploration/config.yaml is extended by additional config specified via the command line.
    WARNING: Duplicated columns: {'metric': ['methods', 'metrics']}
    Building DAG of jobs...
    Job stats:
    job                    count
    -------------------  -------
    qc_all                     1
    qc_autoqc                  2
    qc_get_thresholds          2
    qc_merge_thresholds        1
    qc_plot_joint              2
    qc_plot_removed            2
    total 

Execute the QC workflow as follows:

.. code-block:: bash

    bash configs/qc/run.sh qc_all -c2


This will estimate QC thresholds using the `sctk <https://teichlab.github.io/sctk/notebooks/automatic_qc.html>`_'s `scAutoQC` function.

Output
******

Check the outputs under `images/qc`.

.. code-block:: bash

    $ ls -l images/qc/dataset\~example_qc_analysis/

    total 16
    drwxr-xr-x@ 6 michaela.mueller  2125895594  192 Jun 30 16:07 file_id~test1
    drwxr-xr-x@ 6 michaela.mueller  2125895594  192 Jun 30 16:08 file_id~test2
    -rw-r--r--@ 1 michaela.mueller  2125895594   60 Jun 30 16:06 qc_stats.tsv
    -rw-r--r--@ 1 michaela.mueller  2125895594  829 Jun 30 16:06 thresholds.tsv

There are per dataset and per file outputs.
The `thresholds.tsv` contains the estimated thresholds for both input files `test1` and `test2`, while `qc_stats.tsv` contains the overall statistics of cells that passed or failed the QC thresholds.


.. literalinclude:: ../_static/images/qc/dataset~example_qc_analysis/qc_stats.tsv
   :language: tsv
   :caption: qc_stats.tsv

.. literalinclude:: ../_static/images/qc/dataset~example_qc_analysis/thresholds.tsv
   :language: tsv
   :caption: thresholds.tsv


There are also separate outputs per original input file, which contains per `file_id` threshold files as well as `joint_plots` and `removed` directories.

Joint plots
===========

.. code-block:: bash

    $ ls images/qc/dataset~example_qc_analysis/file_id~test1/

    joint_plots  qc_stats.tsv  removed  thresholds.tsv


.. code-block:: bash

    $ ls images/qc/dataset\~example_qc_analysis/file_id\~test1/joint_plots

    'hue=batch.svg'  'hue=bulk_labels.svg'  'hue=phase.svg'  'hue=qc_status.svg'


Each scatter plot output is now a single SVG file containing all configured QC metric pairs, with linear-scale and log-scale panels arranged side by side.
The workflow writes one scatter plot per configured hue and always adds `hue=qc_status.svg`.
If no valid hue columns are available, the workflow writes `main.svg` instead of hue-specific scatter plots.
The density plots are also combined into a single `density.svg` file rather than separate files per metric pair.

The `removed` directory contains the QC summary plots for the same `file_id`:

.. code-block:: bash

    $ ls images/qc/dataset\~example_qc_analysis/file_id\~test1/removed

    'by=batch.svg'  'by=bulk_labels.svg'  'by=phase.svg'   cells_passed_all.svg   per_metric_violin.svg

The `by={hue}.svg` plots summarise removed cells per annotation, `cells_passed_all.svg` shows the overall QC status counts, and `per_metric_violin.svg` shows the per-metric distributions split by QC status.

Example outputs are shown below:

.. image:: ../_static/images/qc/dataset~example_qc_analysis/file_id~test1/joint_plots/hue=bulk_labels.svg
   :alt: bulk_labels_joint_plot


TL;DR Full Configuration
************************

You can find the complete configuration file and runner script under `configs/qc/`.
Here's the final workflow configuration:

.. literalinclude:: ../../configs/qc/qc_config.yaml
   :language: yaml
   :caption: configs/qc/qc_config.yaml

.. literalinclude:: ../../configs/qc/run.sh
   :language: yaml
   :caption: configs/qc/run.sh

