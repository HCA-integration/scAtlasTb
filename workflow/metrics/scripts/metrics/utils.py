import numpy as np
import scanpy as sc
import pandas as pd


def write_metrics(filename, metric_names, scores, output_types, **kwargs):
    """
    Write metrics for output type specific scores
    :param filename: file to write to
    :param metric_names: list of metric names
    :param scores: list of scores, order must match that of metric_names
    :param output_types: list of output types, order must match that of metric_names
    :param kwargs: additional information to add to output
    """
    meta_names = list(kwargs.keys())
    meta_values = [kwargs[col] for col in meta_names]

    records = [
        (*meta_values, metric, output_type, score)
        for metric, score, output_type
        in zip(metric_names, scores, output_types)
    ]
    df = pd.DataFrame.from_records(
        records,
        columns=meta_names + ['metric_name', 'output_type', 'score']
    )
    df = df.explode('score')
    df.to_csv(filename, sep='\t', index=False)


def rename_categories(adata, obs_col):
    s = adata.obs[obs_col].astype('category')
    s = s.cat.rename_categories({i for i, _ in enumerate(s.cat.categories)})
    return s.to_numpy()


def select_neighbors(adata, output_type):
    # neighbors_key = f'neighbors_{output_type}'
    neighbors_key = 'neighbors'
    adata.uns['neighbors'] = adata.uns[neighbors_key]
    
    connectivities_key = adata.uns[neighbors_key]['connectivities_key']
    assert connectivities_key in adata.obsp, f'Connectivities key "{connectivities_key}" missing from adata.obsp {adata}'
    adata.obsp['connectivities'] = adata.obsp[connectivities_key]
    
    distances_key = adata.uns[neighbors_key]['distances_key']
    assert distances_key in adata.obsp, f'Distances key "{distances_key}" missing from adata.obsp {adata}'
    adata.obsp['distances'] = adata.obsp[distances_key]
    
    return adata


def scanpy_to_neighborsresults(adata):
    from scib_metrics.nearest_neighbors import NeighborsResults
    
    n_neighbors = adata.uns["neighbors"]["params"]["n_neighbors"]
    n_obs = adata.n_obs
    
    # convert sparse to dense neighbor lists
    distances = np.full((n_obs, n_neighbors), np.inf)
    indices = np.full((n_obs, n_neighbors), -1)

    # use scanpy's stored distances (CSR sparse)
    dist_matrix = adata.obsp["distances"].tocsr()

    for i in range(n_obs):
        row = dist_matrix[i].indices
        dists = dist_matrix[i].data
        order = np.argsort(dists)[:n_neighbors]  # sort neighbors by distance
        indices[i, :len(order)] = row[order]
        distances[i, :len(order)] = dists[order]

    return NeighborsResults(indices=indices, distances=distances)
