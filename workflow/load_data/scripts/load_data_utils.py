from pathlib import Path

# CELLxGENE columns of schema 3.0.0
# https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/3.0.0/schema.md#general-requirements

SCHEMAS = {
    "MINIMAL": [ # these are amongst CELLxGENE_OBS + EXTRA_COLUMNS
        "study", # TSV; exploration, benchmark (config)
        "dataset", # TSV; metrics, exploration, subset, integration (?), plots
        "organ", # TSV; exploration
        "donor", # TSV; benchmark (config)
        "sample", # TSV; subset
        "barcode", # TSV (depends on annotation_file); exploration
        # CxG or from author_annotation; exploration, benchmark (config), label_transfer (config)
        "cell_type",
        "author_annotation", # TSV; exploration
    ],
    "CELLxGENE_OBS": [
        'assay',
        'assay_ontology_term_id',
        'cell_type',
        'cell_type_ontology_term_id',
        'development_stage',
        'development_stage_ontology_term_id',
        'disease',
        'disease_ontology_term_id',
        'donor_id',
        'is_primary_data',
        'organism',
        'organism_ontology_term_id',
        'self_reported_ethnicity',
        'self_reported_ethnicity_ontology_term_id',
        'sex',
        'sex_ontology_term_id',
        'suspension_type',
        'tissue',
        'tissue_ontology_term_id',
    ],
    "CELLxGENE_VARS": [
        'feature_name', # load_data, exploration, label transfer
        'feature_reference',
        'feature_biotype',
    ],
    "EXTRA_COLUMNS": [
        'barcode',
        'organ',
        'tech_id',
        'batch_condition',
        'author_annotation',
        'reference',
        'study',
        'dataset',
        'modalities',
        'pipeline_version',
        'institution',
        'remapping_sample',
        'remapping_sample_barcode',
        'barcode_old',
        'age',
    ],
    "TIER1": [
        'sample_id',
        'protocol_url',
        'institute',
        'sample_collection_site',
        'sample_collection_relative_time_point',
        'library_id',
        'library_id_repository',
        'author_batch_notes',
        'manner_of_death',
        'sample_source',
        'sample_collection_method',
        'sampled_site_condition',
        'tissue_type',
        'tissue_ontology_term_id',
        'tissue_free_text',
        'sample_preservation_method',
        'cell_enrichment',
        'cell_viability_percentage',
        'cell_number_loaded',
        'sample_collection_year',
        'library_preparation_batch',
        'library_sequencing_run',
        'sequenced_fragment',
        'sequencing_platform',
        'reference_genome',
        'gene_annotation_version',
        'alignment_software',
        'intron_inclusion',
        'author_cell_type',
    ]
}


def get_union(*args: list):
    """

    :param args: lists of values
    :return: list containing union of all list values
    """
    for arg in args:
        assert isinstance(arg, list)

    sets = [set(arg) for arg in args]
    if len(sets) <= 1:
        return list(sets[0])

    union_set = sets[0]
    for _set in sets[1:]:
        union_set = union_set.union(_set)

    return list(union_set)


def get_from_dataset(dataset_df, key, value, column=None, debug=False):
    """
    Retrieve values in column from dataset_df by subset (key, value)

    :param dataset_df:
    :param key: column in dataset_df to subset to
    :param value: value to subset to
    :param column: column in dataset_df
    :param debug:
    :return:
    """
    if column is None:
        column = dataset_df.columns
    sub = dataset_df.query(f'{key} == @value')
    if debug:
        print(f'##DEBUG## {key, value, column}')
        print(sub)
        print(sub[column])
    return sub[column]


def get_url(dataset_df, wildcards):
    """
    Check if URL exists
    """
    url = get_from_dataset(dataset_df, 'dataset', wildcards.dataset, 'url').iloc[0]
    if url in ['', 'nan', None]:
        raise ValueError(
            f'File location is empty for {wildcards.dataset}. Make sure it is either a valid URL or one of the valid schemas'
        )
    #if not isinstance(url, str):
    if url in ['CxG', 'CELLxGENE', 'cxg', 'cellxgene']:
        collection_id = get_from_dataset(dataset_df, 'dataset', wildcards.dataset, 'collection_id').iloc[0]
        dataset_id = get_from_dataset(dataset_df, 'dataset', wildcards.dataset, 'dataset_id').iloc[0]
        return get_cxg_url(collection_id, dataset_id)
    elif url in ['DCP', 'dcp', 'HCA', 'hca', 'HCA DCP', 'hca dcp']:
        project_uuid = get_from_dataset(dataset_df, 'dataset', wildcards.dataset, 'project_uuid').iloc[0]
        return get_dcp_url(project_uuid)
    return url


def get_cxg_url(collection_id, dataset_id):
    """
    Adapted code from https://github.com/chanzuckerberg/single-cell-curation/blob/main/notebooks/curation_api/python_raw/get_dataset.ipynb
    """
    import requests
    print(f'Get URL for CxG collection ID "{collection_id}" and dataset ID "{dataset_id}"')

    url = f'https://api.cellxgene.cziscience.com/curation/v1/collections/{collection_id}/datasets/{dataset_id}/'
    assets = requests.get(url=url).json()['assets']
    asset = [a for a in assets if a["filetype"] == "H5AD"][0]
    url = asset["url"]

    print(f'URL: {url}')
    return url


def get_dcp_url(project_uuid, catalog='dcp25'):
    """
    Adapted code from https://github.com/DataBiosphere/azul/blob/develop/docs/download-project-matrices.ipynb
    """
    import requests
    
    
    def iterate_matrices_tree(tree, keys=()):
        if isinstance(tree, dict):
            for k, v in tree.items():
                yield from iterate_matrices_tree(v, keys=(*keys, k))
        elif isinstance(tree, list):
            for file in tree:
                yield keys, file
        else:
            assert False


    print(f'Get URL for HCA DCP project UUID "{project_uuid}"')

    endpoint_url = f'https://service.azul.data.humancellatlas.org/index/projects/{project_uuid}'
    response = requests.get(endpoint_url, params={'catalog': catalog})
    response.raise_for_status()

    response_json = response.json()
    project = response_json['projects'][0]

    file_urls = set()
    tree = project['matrices']
    for _, file_info in iterate_matrices_tree(tree):
        url = file_info['url']
        file_type = file_info['format']
        file_urls.add(url)

    if len(file_urls) > 1:
        ValueError(f'Found multiple matrices for project {project_uuid} in catalog {catalog}')
    elif len(file_urls) == 0:
        ValueError(f'Found no matrices for project {project_uuid} in catalog {catalog}')

    print(f'URL: {url}')
    print(f'Format: {file_type}')
    return url, file_type
