def get_threads_for_plotting(wildcards, use_centroids=False):
    colors = mcfg.get_from_parameters(wildcards, 'colors', default=[])
    centroids = mcfg.get_from_parameters(wildcards, 'plot_centroids', default=[])
    
    n_colors = len(colors)
    if use_centroids:
        n_colors += len(centroids)
    
    return min(10, max(1, n_colors))


use rule plots from preprocessing as preprocessing_plot_pca with:
    input:
        anndata=rules.preprocessing_pca.output.zarr,
    output:
        plots=directory(mcfg.image_dir / paramspace.wildcard_pattern / 'pca'),
    threads:
        get_threads_for_plotting
    params:
        color=lambda w: mcfg.get_from_parameters(w, 'colors', default=[]),
        basis='X_pca',
        ncols=1,


use rule plots from preprocessing as preprocessing_plot_umap with:
    input:
        anndata=rules.preprocessing_umap.output.zarr,
    output:
        plots=directory(mcfg.image_dir / paramspace.wildcard_pattern / 'umap'),
    threads:
        lambda wildcards: get_threads_for_plotting(wildcards, use_centroids=True),
    params:
        color=lambda w: mcfg.get_from_parameters(w, 'colors', default=[]),
        plot_centroids=lambda w: mcfg.get_from_parameters(w, 'plot_centroids', default=[]),
        gene_chunk_size=lambda w: mcfg.get_from_parameters(w, 'plot_gene_chunk_size', default=12),
        basis='X_umap',
        ncols=1,
        outlier_factor=100,
    resources:
        mem_mb=lambda w, attempt: 32000 * attempt,
