import numpy as np
import pandas as pd
from pprint import pformat

from utils.ModuleConfig import ModuleConfig
from utils.marker_genes import get_marker_gene_set


class MetricsConfig(ModuleConfig):

    def __init__(
        self,
        **kwargs
    ):
        """
        :param kwargs: parameters for ModuleConfig
        """
        super().__init__(**kwargs, write_output_files=False)

        wildcards_df = self.parameters.wildcards_df
        
        # check if all metrics in config are configured in parameters.tsv
        available_metrics = self.get_parameters()['metric'].unique()
        for metric in wildcards_df['metric'].unique():
            assert metric in available_metrics, f'Metric not available: "{metric}"'
        
        wildcards_df = self.set_gene_sets(wildcards_df)
        wildcards_df = self.set_covariates(wildcards_df)
        wildcards_df = self.set_celltypes_list(wildcards_df)
        wildcards_df = self.set_GT_traj_key(wildcards_df)
        
        self.update_parameters(
            wildcards_df=wildcards_df,
            wildcard_names=self.parameters.wildcard_names,
            **self.parameters.paramspace_kwargs
        )


    def set_gene_sets(self, wildcards_df, config_key='gene_set'):
        gene_set_map = dict()
        wildcards_df[config_key] = wildcards_df[config_key].astype(object)

        for dataset in self.datasets:
            gene_sets = get_marker_gene_set(
                self,
                wildcards={'dataset': dataset},
                config_key=config_key,
                flatten=True
            )
            gene_set_map[dataset] = gene_sets
            
            # set gene sets only for metrics that use it
            dataset_mask = wildcards_df['dataset'] == dataset
            use_config_mask = dataset_mask & wildcards_df[f'use_{config_key}']

            wildcards_df.loc[dataset_mask, config_key] = None
            wildcards_df.loc[use_config_mask, config_key] = pd.Series(
                [list(gene_sets.keys())] * use_config_mask.sum(),
                index=use_config_mask[use_config_mask].index,
                dtype=object,
            )

        self.gene_set_map = gene_set_map

        # explode gene sets so that each metric setting has a single gene set
        wildcards_df = wildcards_df.explode(config_key, ignore_index=True)
        wildcards_df[config_key] = wildcards_df[config_key].fillna('None')

        return wildcards_df


    def get_gene_sets(self, dataset, **kwargs):
        return self.gene_set_map[dataset]


    def set_covariates(self, wildcards_df, config_key='covariate'):
        # covariates only used for metrics that require it
        wildcards_df.loc[~wildcards_df[f'use_{config_key}'], config_key] = None

        # # explode covariates so that each metric setting has a single covariate
        # wildcards_df = wildcards_df.explode(config_key, ignore_index=True)
        
        wildcards_df[config_key] = wildcards_df[config_key].fillna('None')

        return wildcards_df

    def set_celltypes_list(self, wildcards_df, config_key='celltypes_list'):
        # celltypes_list only used for metrics that require it
        wildcards_df.loc[~wildcards_df[f'use_{config_key}'], config_key] = None

        wildcards_df[config_key] = wildcards_df[config_key].fillna('None')

        return wildcards_df

    def set_GT_traj_key(self, wildcards_df, config_key='GT_traj_key'):
        # GT_traj_key only used for metrics that require it
        wildcards_df.loc[~wildcards_df[f'use_celltypes_list'], config_key] = None

        wildcards_df[config_key] = wildcards_df[config_key].fillna('None')

        return wildcards_df

