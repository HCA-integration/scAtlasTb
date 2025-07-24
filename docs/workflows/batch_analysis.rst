Quantify batch effect
=====================

Understanding the existing batch effect in your data is crucial for achieving a good integration.

Using principle

.. code-block:: yaml

      example_batch_analysis:
        input:
          preprocessing:
            custom_file_name: custom/file/path.zarr # file path
          batch_analysis: preprocessing
        batch_analysis:
          sample: sample # smallest entity of a batch, e.g. <bio_sample>-<pool> if there is no 1-1 matching of sample to pool
          n_permutations: 100 # you have small data, so should be quick
          covariates:
            - # obs columns of potential batch covariates
        preprocessing:
          highly_variable_genes:
            n_top_genes: 2000
            # no batch_key bc we want to maximise batch effect
          assemble: # only pca is needed here
            - pca
          colors: # obs columns that you want in your UMAP
            - cell_type


