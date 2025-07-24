Iterative clustering
====================

Setting up the environment
--------------------------

Make sure you have the necessary environments installed.
For the iterative clustering workflow, you will need the following:

* `envs/scanpy.yaml`
* `envs/rapids_singlecell.yaml` (optional, if you have NVIDIA GPU support)


Configuring the workflow
------------------------

Set up your configuration file


.. code-block:: yaml

    output_dir: data/out
    images: images
    user_gpu: false

    DATASETS:
        
        my_dataset: # custom task/workflow name
        
            # input specification: map of module name to map of input file name to input file path
            input:
            clustering:
                file_1: data/pbmc68k.h5ad
                file_2: data/pbmc68k.h5ad # dummy example here for more than 1 file
        
            # module configuration
            clustering:
              recompute_neighbors: true
              neighbors:
                n_neighbors: 30
                use_rep: X_pca
              recompute_clusters: true
              algorithm: leiden
              resolutions:
                - 1.0
              hierarchy:
                1: 0.1
                3: 0.2
            umap_colors:
              - bulk_labels
              - batch
              - n_genes


Calling the pipeline
--------------------

Make sure you have set up a runner script to call the Snakemake workflow.
The following command will do a dry run of all the steps in the clustering workflow:

.. code-block:: bash

    bash run_clustering.sh clustering_all -nq


There are also optional evaluation plots that you can call:

.. code-block:: bash

    bash run_clustering.sh clustering_plot_evaluation_all -nq

