import pandas as pd

from utils.ModuleConfig import ModuleConfig
from utils.misc import unique_dataframe, expand_dict_and_serialize


class PreprocessingConfig(ModuleConfig):
    
    HVG_PARAMS = [
        'n_top_genes',
        'min_disp',
        'max_disp',
        'min_mean',
        'max_mean',
        'span',
        'n_bins',
        'flavor',
        'batch_key'
    ]
    
    def __init__(
        self,
        **kwargs
    ):
        """
        :param kwargs: parameters for ModuleConfig
        """
        super().__init__(**kwargs)
        self.set_highly_variable_genes()
        self.set_extra_hvgs()

    
    def _process_hvg_config(self, config, allow_none_dict=True):
        """
        Process HVG config and filter by HVG_PARAMS.
        
        :param config: HVG configuration (dict, None, False, etc.)
        :param allow_none_dict: whether to allow None or non-dict to return ('None', None)
        :return: tuple of (serialized_string, config_dict) or list of tuples if expandable
        """
        # Preserve False value (no HVG)
        if config is False:
            return [('False', False)]
        elif config is None or not isinstance(config, dict):
            return [('None', None)] if allow_none_dict else []
        else:
            _dict = {
                k: v for k, v in config.items()
                if k in self.HVG_PARAMS
            }
            if _dict:  # only expand if there are parameters
                return list(expand_dict_and_serialize(_dict))
            else:
                return [('None', None)]
    
    def _update_wildcards_with_records(self, records, columns):
        """
        Helper to create DataFrame from records and update parameters.
        
        :param records: list of tuples with dataset and config data
        :param columns: column names for DataFrame
        """
        df = pd.DataFrame(records, columns=columns)
        wildcards_df = unique_dataframe(
            self.parameters.wildcards_df.merge(df, on='dataset', how='left')
        )
        self.update_parameters(wildcards_df=wildcards_df)

    
    def set_highly_variable_genes(self):
        """Create multiple HVG configs from parameter combinations."""
        wildcards_df = self.parameters.wildcards_df
        records = []
        for dataset in self.get_datasets():
            hvg_config = wildcards_df.query('dataset == @dataset')['highly_variable_genes'].iloc[0]
            processed = self._process_hvg_config(hvg_config)
            records.extend((dataset, *rec) for rec in processed)
        
        self._update_wildcards_with_records(
            records,
            columns=['dataset', 'hvg_args', 'hvg_args_dict']
        )

    
    def set_extra_hvgs(self):
        wildcards_df = self.parameters.wildcards_df
        records = []
        for dataset in self.get_datasets():
            ehvg_config = wildcards_df.query('dataset == @dataset')['extra_hvgs'].iloc[0]
            overwrite_args = ehvg_config.get('overwrite_args') if isinstance(ehvg_config, dict) else None
            processed = self._process_hvg_config(overwrite_args, allow_none_dict=True)
            records.extend((dataset, *rec) for rec in processed)
        
        self._update_wildcards_with_records(
            records,
            columns=['dataset', 'overwrite_args', 'overwrite_args_dict']
        )
                
            