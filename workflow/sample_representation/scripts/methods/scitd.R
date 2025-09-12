options(show.error.locations = TRUE)

tryCatch({
 suppressPackageStartupMessages(library(scITD))
  }, error = function(e) {
    install.packages('scITD', repos = snakemake@params$cran_url)
  }
)

# Load required libraries
suppressPackageStartupMessages({
  library(data.table)
  library(scITD)
  library(reticulate)
  library(anndata)
})
io <- import('utils.io')
misc <- import('utils.misc')
sc <- import('scanpy')


run_scITD <- function(
    data_container,
    n_factors = 5,
    n_gene_sets = 10,
    tucker_type = 'regular',
    rotation_type = 'hybrid'
) {
  
  ranks = c(n_factors, n_gene_sets)
  
  data_container <- run_tucker_ica(
    data_container,
    ranks = ranks,
    tucker_type = tucker_type,
    rotation_type = rotation_type
  )
  
  all_scores <- list()
  all_loadings <- list()
  
  # Iterate over the factors to collect all the scores and loadings
  for (i in 1:n_factors) {
    f_data <- get_one_factor(data_container, factor_select=i)
    f_scores <- f_data[[1]]
    all_scores[[i]] <- f_scores
    
    f_loadings <- f_data[[2]]
    colnames(f_loadings) <- paste0("Factor", i, "_", colnames(f_loadings))
    all_loadings[[i]] <- f_loadings
  }
  
  # Combine all scores into a single data frame
  scores_table <- do.call(cbind, all_scores)
  colnames(scores_table) <- paste0("Factor", 1:n_factors)
  
  loadings_table <- do.call(cbind, all_loadings)
  
  list(
    data_container = data_container,
    scores = scores_table,
    loadings = loadings_table
  )
}


input_file <- snakemake@input$zarr
bulk_file <- snakemake@input$bulks
output_file <- snakemake@output$zarr

cell_type_key <- snakemake@params$cell_type_key
use_rep <- snakemake@params$use_rep
var_mask <- snakemake@params$var_mask
threads <- snakemake@params$threads
seed <- snakemake@params$seed
norm_method <- 'trim'
scale_factor <- 1e4
var_scale_power <- 2
scale_var <- TRUE

message("Read data...")
n_obs <- io$read_anndata(input_file, obs='obs')$n_obs
dask <- n_obs > 1e6
adata <- io$read_anndata(
    input_file,
    X=use_rep,
    obs='obs',
    backed=dask,
    dask=dask,
    stride=as.integer(n_obs / 5)
)

if (use_rep == 'X' | grepl('^layers/|^raw/', use_rep)) {
    adata$var <- io$read_anndata(input_file, var='var')$var
}

# subset HVGs
if (!is.null(var_mask)){
    adata <- adata[, adata$var[var_mask]$values]
}
misc$dask_compute(adata)

message("Prepare data and parameters...")
param_list <- initialize_params(
    ctypes_use = as.character(unique(adata$obs[[cell_type_key]])),
    ncores = threads,
    rand_seed = seed
)
rownames(adata$X) <- adata$obs_names
colnames(adata$X) <- adata$var_names
rownames(adata$obs) <- adata$obs_names
adata$obs[['donors']] <- adata$obs[['group']]
adata$obs[['ctypes']] <- adata$obs[[cell_type_key]]

message("Create data container...")
data_container <-  make_new_container(
    count_data = t(adata$X),
    meta_data = adata$obs,
    # gn_convert = adata$var_names,
    params = param_list,
    label_donor_sex = FALSE
  )

message("Forming a tensor")
data_container <- form_tensor(
    data_container, 
    donor_min_cells = 1,
    norm_method = norm_method,
    scale_factor = scale_factor,
    vargenes_method = 'norm_var_pvals',
    vargenes_thresh = 0.1,
    scale_var = scale_var,
    var_scale_power = var_scale_power
)
n_hvgs <- length(data_container[["all_vargenes"]])
if (n_hvgs == 0) {
    data_container[['all_vargenes']] <- rownames(adata$X)
}
message(paste("Number of HVGs:", length(data_container[["all_vargenes"]])))

# message("Running analysis for rank determination")
# data_container <- determine_ranks_tucker(
#     data_container,
#     max_ranks_test = c(10, 5),
#     shuffle_level = 'cells', 
#     num_iter = 10, 
#     norm_method = norm_method,
#     scale_factor = scale_factor,
#     scale_var = scale_var,
#     var_scale_power = var_scale_power
# )

results <- run_scITD(
    data_container,
    n_factors = 5,
    n_gene_sets = 10,
)
pairwise_distances <- dist(results$scores, method = "euclidean")

# create new anndata object
message("Creating new anndata object...")
adata <- sc$AnnData(
    obs = data.frame(row.names = unique(sample_ids)),
    obsm = list(
      distances = dist_matrix,
      X_pca = sc$tl$pca(dist_matrix)
    )
)
samples <- io$read_anndata(prepare_file, obs='obs')$obs_names
adata <- adata[samples]

# compute kNN graph
sc$pp$neighbors(adata, use_rep = 'distances', metric = 'precomputed')

message(paste("Writing to", output_file, "..."))
print(adata)
io$write_zarr_linked(
    adata,
    in_dir=prepare_file,
    out_dir=output_file,
    files_to_keep=c('obsm', 'obsp', 'uns'),
)