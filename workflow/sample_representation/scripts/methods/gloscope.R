options(show.error.locations = TRUE)

# Load required libraries
suppressPackageStartupMessages({
  library(data.table)
  library(GloScope)
  library(reticulate)
  library(anndata)
})
io <- import('utils.io')
misc <- import('utils.misc')
sc <- import('scanpy')

input_file <- snakemake@input$zarr
prepare_file <- snakemake@input$prepare
output_file <- snakemake@output$zarr

sample_key <- snakemake@params$sample_key
use_rep <- snakemake@params$use_rep
var_mask <- snakemake@params$var_mask
k <- snakemake@params$k
seed <- snakemake@params$seed
threads <- snakemake@threads


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

# subset HVGs
if (!is.null(var_mask)){
    adata$var <- io$read_anndata(input_file, var='var')$var
    adata <- adata[, adata$var[var_mask]$values]
}
misc$dask_compute(adata)

embedding <- adata$X
rownames(embedding) <- adata$obs_names
sample_ids <- adata$obs[[sample_key]] # TOOD: is this correct?
# obs_dt <- data.table(adata$obs)
# sample_dt <- obs_dt[, lapply(.SD, first), by = sample_key]

message(paste(c("Data dimensions:", paste(dim(embedding), collapse=', '))))
message(paste("Number of IDs:", uniqueN(sample_ids)))

message("Calculating distance matrix...")
dist_matrix <- gloscope(
    embedding,
    sample_ids,
    dens = "GMM",
    dist_mat = "KL",
    k = k,
    BPPARAM = BiocParallel::MulticoreParam(workers = 4, RNGseed = seed)
)
message(paste(c("Distance matrix dimensions:", paste(dim(dist_matrix), collapse=', '))))

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