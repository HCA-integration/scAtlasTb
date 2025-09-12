import numpy as np
import scanpy as sc
from .utils import select_neighbors, rename_categories, scanpy_to_neighborsresults


def ari(adata, output_type, batch_key, label_key, cluster_key, **kwargs):
    import scib

    adata = select_neighbors(adata, output_type)
    scib.cl.cluster_optimal_resolution(
        adata=adata,
        label_key=label_key,
        cluster_key=cluster_key,
        metric=scib.me.nmi,
        use_rep=None,
    )
    return scib.me.ari(
        adata[adata.obs[label_key].notna()],
        label_key,
        cluster_key
    )


def ari_leiden_y(adata, output_type, batch_key, label_key, **kwargs):
    import scib_metrics

    adata = select_neighbors(adata, output_type)
    labels = rename_categories(adata, label_key)

    scores = scib_metrics.nmi_ari_cluster_labels_leiden(
        X=scanpy_to_neighborsresults(adata),
        labels=labels,
        optimize_resolution=True,
    )
    return scores['ari']


def ari_kmeans_y(adata, output_type, batch_key, label_key, **kwargs):
    from scib_metrics import nmi_ari_cluster_labels_kmeans
    
    labels = rename_categories(adata, label_key)
    adata = select_neighbors(adata, output_type)

    scores = scib_metrics.nmi_ari_cluster_labels_kmeans(
        X=scanpy_to_neighborsresults(adata),
        labels=labels,
    )
    return scores['ari']
