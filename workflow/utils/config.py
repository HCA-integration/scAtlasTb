from pprint import pprint
from typing import Union
import warnings
import pandas as pd
import hashlib

from .misc import expand_dict_and_serialize


def set_defaults(config, modules=None, warn=False):
    if 'defaults' not in config:
        config['defaults'] = {}
    if 'datasets' not in config['defaults']:
        config['defaults']['datasets'] = list(config['DATASETS'].keys())

    if modules is None:
        modules = ['integration', 'metrics']
    elif isinstance(modules, str):
        modules = [modules]

    for module in modules:

        # initialise if module defaults not defined
        if module not in config['defaults']:
            config['defaults'][module] = {}

        # update entries for each dataset
        for dataset in config['DATASETS'].keys():
            entry = _get_or_default_from_config(
                config=config['DATASETS'],
                defaults=config['defaults'],
                key=dataset,
                value=module,
                return_missing={},
                warn=warn,
            )

            # for TSV input make sure integration methods have the proper types
            # if module == 'integration' and isinstance(entry, list):
            #     # get parameters from config
            #     entry = {k: config['defaults'][module][k] for k in entry}

            # set entry in config
            config['DATASETS'][dataset][module] = entry

    return config


def get_params_from_config(
        config,
        module_name,
        config_params,
        wildcard_names,
        defaults,
        explode_by=None,
        config_keys=None,
        warn=False,
):
    """
    Collect wildcards and parameters from an configuration instance (e.g. dataset) for a given module.
    This function assumes that the keys of the given config keys are the different instances that contain specific parameters for different modules.

    :param config: Part of the Snakemake config dictionary. The config_params must be contained per entry.
    :param module_name: Name of the module to extract the config from
    :param config_params: List of parameters for each config entry of a module
        e.g. ['integration', 'label', 'batch']
    :param wildcard_names: names of wildcards to be extracted.
        Must map to config keys, and prepended by a wildcard name for the config entries
        e.g. ['dataset', 'method', 'label', 'batch']
    :param explode_by: column to explode by, expecting list entry for that column
    :param config_keys: list of entries to subset the config by., otherwise use all keys
    :return: dataframe with wildcard mapping. Wildcard names in columns and wildcard values as entries
    """
    if not config:
        return pd.DataFrame(columns=[*wildcard_names])

    # set default config keys
    if config_keys is None:
        config_keys = config.keys()

    # collect entries for dataframe
    records = [
        (
            key,
            *[
                _get_or_default_from_config(
                    config=config[key],
                    defaults=defaults[module_name],
                    key=module_name,
                    value=param,
                    warn=warn,
                    )
                for param in config_params
            ]
        )
        for key in config_keys
    ]
    df = pd.DataFrame.from_records(records, columns=[*wildcard_names])
    if explode_by is not None:
        explode_by = [explode_by] if isinstance(explode_by, str) else explode_by
        for column in explode_by:
            df = df.explode(column)
    return df.reset_index(drop=True)

def get_wildcards_from_config(
        config,
        config_params,
        wildcard_names,
        explode_by=None,
        config_keys=None,
):
    """
    Deprecated

    :param config: Part of the Snakemake config dictionary. The config_params must be contained per entry.
    :param config_params: list oxf parameters for each config entry
        e.g. ['integration', 'label', 'batch']
    :param wildcard_names: names of wildcards to be extracted.
        Must map to config keys, and prepended by a wildcard name for the config entries
        e.g. ['dataset', 'method', 'label', 'batch']
    :param explode_by: column to explode by, expecting list entry for that column
    :param config_keys: list of entries to subset the config by., otherwise use all keys
    :return: dataframe with wildcard mapping. Wildcard names in columns and wildcard values as entries
    """
    if not config:
        return pd.DataFrame(columns=[*wildcard_names])

    if config_keys is None:
        config_keys = config.keys()
    records = [
        (key, *[_get_or_default_from_config(config, {}, key, w) for w in config_params])
        for key in config_keys
    ]
    df = pd.DataFrame.from_records(records, columns=[*wildcard_names])
    if explode_by is not None:
        explode_by = [explode_by] if isinstance(explode_by, str) else explode_by
        for column in explode_by:
            df = df.explode(column)
    return df.reset_index(drop=True)


def _get_or_default_from_config(
    config,
    defaults,
    key,
    value,
    return_missing=None,
    dont_inherit=None,
    warn=False,
):
    """
    Get entry from config or return defaults if not present

    :param config: part of the config with multiple entries of the same structure
    :param defaults: config defaults for keys with missing value
    :param key: top-level key of config dictionary
    :param value: points to entry within key
    :param return_missing: return value if no defaults for key, value
    :param dont_inherit: list of keys in value config that should not be inherited from defaults
    :param warn: warn if no defaults are defined for key, value
    :return:
    """
    if dont_inherit is None:
        dont_inherit = []
    
    if key not in config.keys():
        print('config:')
        pprint(config)
        print(f'key: {key}, value: {value}')
        raise KeyError(f'Key "{key}" not found in config')
    
    config_at_value = config.get(key, {})
    if config_at_value is None:
        config_at_value = {}
    
    default_entry = defaults.get(value, return_missing)
    if value in dont_inherit and value in config_at_value.keys():
        default_entry = return_missing
    entry = config_at_value.get(value, default_entry)
    
    # update if entry is a dict
    if isinstance(entry, dict) and isinstance(default_entry, dict):
        entry = default_entry | entry
    
    return entry


def get_hyperparams(config, module_name='integration', methods_key='methods'):
    """
    Get hyperparameters specific to each method of a module for all datasets

    :param config: config containing dataset specific information
    :param module: name of module, key must be present for each dataset entry
    :return: DataFrame with hyperparameters
    """

    records = []
    for dataset, dataset_dict in config['DATASETS'].items():
        methods_config = _get_or_default_from_config(
            config=dataset_dict,
            defaults=config['defaults'][module_name],
            key=module_name,
            value=methods_key,
            update=False,
        )
        if methods_config is None:
            continue
        for method, hyperparams_dict in methods_config.items():
            if isinstance(hyperparams_dict, dict):
                records.extend(
                    (dataset, method, *rec)
                    for rec in expand_dict_and_serialize(hyperparams_dict)
                )
            else:
                records.append((dataset, method, str(hyperparams_dict), hyperparams_dict))
    return pd.DataFrame(records, columns=['dataset', 'method', 'hyperparams', 'hyperparams_dict'])


def get_resource(config, resource_key, profile='cpu', attempt=1, factor=0.5):
    """
    Retrieve resource information from config['resources']
    
    :param config: config passed from Snakemake
    :param profile: resource profile, key under config['resources']
    :param resource_key: resource key, key under config['resources'][profile]
    """
    if 'resources' not in config or not profile:
        return ''
    resources = config['resources']
    try:
        res = resources[profile][resource_key]
    except KeyError:
        print(
            f'WARNING: Invalid profile "{profile}" or resource key "{resource_key}". '
            'Please check that your config contains the correct entries under config["resources"]'
        )
        return ''
    return int(res + (attempt - 1) * factor * res) if resource_key == 'mem_mb' else res

def get_resource_customized(config, config_customized=False, dataset_name=None, resource_key=None, profile='cpu', attempt=1):
    """
    Retrieve resource information from config['resources'_customized']
    :param config_customized=True to get the customized resources
    :param config: config passed from Snakemake
    :param profile: resource profile, key under config['resources']
    :param dataset_name: name of the dataset
    :param resource_key: resource key, key under config['resources'][profile]
    
    """
    if config_customized:
        if 'resources_customized' not in config or 'dataset_names' not in config['resources_customized']:
            return ''
        dataset_names = config['resources_customized']['dataset_names']
        if dataset_name not in dataset_names or profile not in dataset_names[dataset_name]:
            return ''
        try:
            res = dataset_names[dataset_name][profile][resource_key]
        except KeyError:
            print(
                f'WARNING: Invalid profile "{profile}" or resource key "{resource_key}". '
                'Please check that your config contains the correct entries under config["resources_customized"]'
            )
            res = ''
        return int(res * attempt * 1.2) if resource_key == 'mem_mb' else res
    else:
        if 'resources' not in config or profile not in config['resources']:
            return ''
        try:
            res = config['resources'][profile][resource_key]
        except KeyError:
            print(
                f'WARNING: Invalid profile "{profile}" or resource key "{resource_key}". '
                'Please check that your config contains the correct entries under config["resources"]'
            )
            res = ''
        return int(res * attempt * 1.2) if resource_key == 'mem_mb' else res


def get_datasets_for_module(config, module):
    """
    Collect dataset names e.g. for wildcard expansion from config["DATASETS"] for a given module
    If the "DATASETS" key is not available in the config, warn and return an empty list
    A dataset is valid if it contains an input file for the given module and it is included in config['defaults']['datasets'].
    If config['defaults']['datasets'] is not defined, all datasets are valid.

    :param config: config dictionary passed from Snakemake
    :param module: name of module to collect valid datasets for
    :return: subset of config["DATASETS"] that is valid for module
    """
    if 'DATASETS' not in config:
        warnings.warn('No datasets specified in config, cannot collect any datasets')
        return {}

    if 'defaults' in config:
        dataset_config = {
            k: v for k, v in config['DATASETS'].items()
            if k in config['defaults'].get('datasets', config['DATASETS'].keys())
        }
    else:
        dataset_config = config['DATASETS']

    return {
        dataset: entry
        for dataset, entry in dataset_config.items()
        if 'input' in dataset_config[dataset]
        and dataset_config[dataset]['input'] is not None
        and module in dataset_config[dataset]['input'].keys()
    }


def get_for_dataset(
    config: dict,
    dataset: str,
    query: list,
    default: Union[str,bool,float,int,dict,list, None] = None,
    warn: bool = False,
) -> Union[str,bool,float,int,dict,list, None]:
    """ TODO: deprecate
    Get any key from the config via query

    Args:
        config (dict): config passed to Snakemake
        dataset (str): dataset key in config['DATASETS']
        query (list): list of keys to walk down the config
        default (Union[str,bool,float,int,dict,list, None], optional): default value if key not found. Defaults to None.

    Returns:
        Union[str,bool,float,int,dict,list, None]: value of query in config
    """
    try:
        assert 'DATASETS' in config
        assert dataset in config['DATASETS']
    except AssertionError as e:
        raise ValueError(f'Assertion failed for dataset="{dataset}"') from e

    # start at top level
    value = config['DATASETS'][dataset]
    
    # walk down query
    for q in query:
        if q not in value:
            if warn:
                warnings.warn(f'key {q} not found in config for query {query}, returning default')
            return default
        value = value.get(q, default)
    return value


def get_from_config(
    config: dict,
    query: list,
    default: Union[str,bool,float,int,dict,list, None] = None,
    warn: bool = False,
) -> Union[str,bool,float,int,dict,list, None]:
    """Get any key from the config via query

    Args:
        config (str): config dictionary
        query (list): list of keys to walk down the config
        default (Union[str,bool,float,int,dict,list, None], optional): default value if key not found. Defaults to None.

    Returns:
        Union[str,bool,float,int,dict,list, None]: value of query in config
    """
    value = config # start at top level
    for q in query: # walk down query
        try:
            value = value[q]
        except (AttributeError, KeyError, TypeError):
            if warn:
                warnings.warn(f'key {q} not found in config for query {query}, returning default')
            return default
    return value


def get_input_file_per_dataset(config, dataset, module_name, digest_size=5):
    """Get input files for a given module and dataset
    This function maps an input file with its unique identifier
    
    :param config: config passed from Snakemake
    :param dataset: dataset key in config['DATASETS']
    :param module_name: name of module to collect valid datasets for
    :return: dictionary of input file ID to file
    """
    input_files = get_for_dataset(config, dataset, query=['input', module_name])
    
    if isinstance(input_files, str):
        input_files = [input_files]
    
    if isinstance(input_files, list):
        input_files = {
            hashlib.blake2b(file.encode('utf-8'), digest_size=digest_size).hexdigest(): file
            for file in input_files
        }
    
    if not isinstance(input_files, dict):
        raise ValueError(f'input_files must be a list or dict, but is {type(input_files)}')
    
    return input_files


def get_input_file(config, wildcards, module_name):
    return get_input_file_per_dataset(config, wildcards.dataset, module_name)[wildcards.file_id]


def get_input_file_wildcards(config, module_name):
    """
    Reshape the file ID to dataset mapping to a lists of wildcards
    """
    all_wildcards = dict(dataset=[], file_id=[], file_name=[])
    for dataset in get_datasets_for_module(config,module=module_name):
        for file_id, file in get_input_file_per_dataset(config, dataset, module_name).items():
            all_wildcards['dataset'].append(dataset)
            all_wildcards['file_id'].append(file_id)
            all_wildcards['file_name'].append(file)
    return all_wildcards