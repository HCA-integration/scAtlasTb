import numpy as np
import scanpy as sc
from .utils import select_neighbors, rename_categories, scanpy_to_neighborsresults


def nmi(adata, output_type, batch_key, label_key, cluster_key, **kwargs):
    import scib

    adata = select_neighbors(adata, output_type)
    scib.cl.cluster_optimal_resolution(
        adata=adata,
        label_key=label_key,
        cluster_key=cluster_key,
        metric=scib.me.nmi,
        use_rep=None,
    )
    return scib.me.nmi(
        adata=adata[adata.obs[label_key].notna()],
        label_key=label_key,
        cluster_key=cluster_key,
    )


def nmi_leiden_y(adata, output_type, batch_key, label_key, **kwargs):
    import scib_metrics

    adata = select_neighbors(adata, output_type)
    labels = rename_categories(adata, label_key)

    scores = scib_metrics.nmi_ari_cluster_labels_leiden(
        X=scanpy_to_neighborsresults(adata),
        labels=labels,
        optimize_resolution=True,
    )
    return scores['nmi']


def nmi_kmeans_y(adata, output_type, batch_key, label_key, **kwargs):
    import scib_metrics

    labels = rename_categories(adata, label_key)
    adata = select_neighbors(adata, output_type)

    X = adata.obsp['connectivities']
    X = X if isinstance(X, np.ndarray) else X.toarray()

    scores = scib_metrics.nmi_ari_cluster_labels_kmeans(
        X=X,
        labels=labels,
    )
    return scores['nmi']
