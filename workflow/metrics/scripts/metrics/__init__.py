from .ari import ari, ari_kmeans_y, ari_leiden_y
from .asw import asw_batch, asw_batch_y, asw_label, asw_label_y
from .graph_connectivity import graph_connectivity, graph_connectivity_y
from .isolated_labels import isolated_label_asw, isolated_label_f1, isolated_label_asw_y
from .lisi import clisi, clisi_y, ilisi, ilisi_y
from .nmi import nmi, nmi_kmeans_y, nmi_leiden_y
from .pcr import cell_cycle, pcr_comparison, pcr_y, pcr, pcr_random, pcr_genes
from .kbet import kbet_y, kbet_pg
from .morans_i import morans_i, morans_i_random, morans_i_genes, morans_i_genescore

metric_map = {
    'ari': ari,
    'ari_kmeans_y': ari_kmeans_y,
    'ari_leiden_y': ari_leiden_y,
    'asw_batch': asw_batch,
    'asw_batch_y': asw_batch_y,
    'asw_label': asw_label,
    'asw_label_y': asw_label_y,
    'asw_label_y': asw_label_y,
    'graph_connectivity': graph_connectivity,
    'graph_connectivity_y': graph_connectivity_y,
    'isolated_label_asw': isolated_label_asw,
    'isolated_label_asw_y': isolated_label_asw_y,
    'isolated_label_f1': isolated_label_f1,
    'clisi': clisi,
    'clisi_y': clisi_y,
    'ilisi': ilisi,
    'ilisi_y': ilisi_y,
    'clisi': clisi,
    'clisi': clisi,
    'nmi': nmi,
    'nmi_kmeans_y': nmi_kmeans_y,
    'nmi_leiden_y': nmi_leiden_y,
    'cell_cycle': cell_cycle,
    'pcr_comparison': pcr_comparison,
    'pcr_y': pcr_y,
    'kbet_y': kbet_y,
    'kbet_pg': kbet_pg,
    'morans_i_random': morans_i_random,
    'morans_i': morans_i,
    'morans_i_genes': morans_i_genes,
    'morans_i_genescore': morans_i_genescore,
    # 'morans_i_random_genes': morans_i_random_genes,
    'pcr_random': pcr_random,
    'pcr': pcr,
    'pcr_genes': pcr_genes,
    # 'pcr_random_genes': pcr_random_genes,
}
