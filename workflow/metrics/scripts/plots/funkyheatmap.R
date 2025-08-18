input_file <- snakemake@input$tsv
output_tsv <- snakemake@output$tsv
output_pdf <- snakemake@output$pdf

id_vars <- snakemake@params$id_vars
value_var <- snakemake@params$value_var
variable_var <- snakemake@params$variable_var
weight_batch <- snakemake@params$weight_batch
n_top <- snakemake@params$n_top

tryCatch({
 library(funkyheatmap)
  }, error = function(e) {
    install.packages('funkyheatmap', repos = snakemake@params$cran_url)
  }
)
suppressPackageStartupMessages({
  library(funkyheatmap)
  library(data.table)
  library(tidyverse)
  library(dynutils)
})
source(snakemake@input$r_utils)

# read files
dt <- fread(input_file)
dt[, (id_vars) := lapply(.SD, as.character), .SDcols = id_vars]
print(head(dt))

extra_columns <- readLines(snakemake@input$extra_columns)
print(extra_columns)
id_vars <- unique(c(id_vars, extra_columns))

# remove unintegrated output types without corresponding method
if ('unintegrated' %in% dt$method & uniqueN(dt$method) > 1) {
  ot_count <- dt[method != 'unintegrated', .(output_type, method)] %>%
    unique %>%
    .[, output_type] %>%
    table
  keep_ot <- ot_count[ot_count > 1]
  df <- dt[, output_type %in% keep_ot]
}

# define groups of interest to be plotted in funkyheatmap
integration_setup <- id_vars
bio_metrics <- unique(dt[metric_type == 'bio_conservation' & !is.na(score), metric]) 
batch_metrics <- unique(dt[metric_type == 'batch_correction' & !is.na(score), metric])
metrics <- c(bio_metrics, batch_metrics)

# subset data.table to columns of interest & transform
all_columns <- c(id_vars, variable_var, value_var, 'metric_type')
id_formula <- paste(paste(id_vars, collapse='+'), variable_var, sep = '~')
cat(paste(id_formula, '\n'))
metrics_tab <- dcast(
  subset(dt, select = all_columns),
  as.formula(id_formula),
  value.var = value_var,
  fun.aggregate = mean
)
print(head(metrics_tab))

# remove columns that are all empty
metrics_tab <- metrics_tab[, .SD, .SDcols = colSums(is.na(metrics_tab)) < nrow(metrics_tab)]

# scores should be already scaled [0,1] - however, we aim to compute the scores based on the min-max scaled metrics
scaled_metrics_tab <- as.matrix(metrics_tab[, metrics, with=FALSE])
if (nrow(scaled_metrics_tab) > 1) {
  scaled_metrics_tab <- as.data.table(apply(scaled_metrics_tab, 2, function(x) scale_minmax(x)))
} else {
  scaled_metrics_tab <- as.data.table(scaled_metrics_tab)
}
cat('Scaled metrics table: \n')
print(head(scaled_metrics_tab))

# calculate average score by group and 
score_group_batch <- rowMeans(scaled_metrics_tab[, batch_metrics, with=FALSE], na.rm = TRUE)
score_group_bio <- rowMeans(scaled_metrics_tab[, bio_metrics, with=FALSE], na.rm = TRUE)

# weighted overall score
score_all <- (weight_batch * score_group_batch + (1 - weight_batch) * score_group_bio)
if(length(score_group_batch) > 0){
  metrics_tab <- add_column(metrics_tab, "Batch Correction" = score_group_batch, .before = batch_metrics[1])
}
if(length(score_group_bio) > 0){
  metrics_tab <- add_column(metrics_tab, "Bio Conservation" = score_group_bio, .before = bio_metrics[1])
}
if(length(score_all) > 0){
  metrics_tab <- add_column(metrics_tab, "Overall Score" = score_all, .before = "Batch Correction")
  metrics_tab <- metrics_tab[order(metrics_tab$`Overall Score`,  decreasing = T), ]
}
# order methods by the overall score
metrics_tab <- metrics_tab[order(metrics_tab$`Overall Score`,  decreasing = T), ]

# write output summary.csv file before plotting ? Michaela to check the output dir?
cat('Writing output file: ', output_tsv, '\n')
fwrite(metrics_tab, output_tsv, sep='\t')

# subset data
metrics_tab <- metrics_tab[1:min(n_top, nrow(metrics_tab))]

# add funkyheatmap data
row_info <- NULL
row_groups <- NULL
# if ('split_data_value' %in% colnames(metrics_tab)) {
#   ## hardcoded!!!
#   # set split_data_value NaN values to values from file_name column
#   metrics_tab[split_data_value == '', split_data_value := file_name]
#   metrics_tab[, split_data_value := gsub("_", "+", split_data_value)]
#   metrics_tab[split_data_value == 'B+plasma', split_data_value := 'B+Plasma']
#   metrics_tab <- metrics_tab[order(get('split_data_value'))]
#   metrics_tab[, id := rownames(metrics_tab)]
#   row_info <- data.table(id = metrics_tab$id, group = metrics_tab[, get('split_data_value')]) #  metrics_tab$output_type)
#   # metrics_tab[, output_type := NULL]
#   row_groups <- data.table(group=unique(row_info$group), Group=unique(row_info$group))
#   metrics_tab[, split_data_value := NULL]
# }

# # shorten redundant information in long names
# cols <- colnames(metrics_tab)
# for (col in c('file_name', 'file_id')) {
#   col_values <- metrics_tab[[col]]
#   max_len <- max(nchar(na.omit(col_values)), na.rm = TRUE)
#   if (!is.na(max_len) && (max_len > 20) && !is.numeric(col_values)) {
#     lcs <- longest_common_substring(col_values)
#     metrics_tab[[col]] <- gsub(lcs, '', metrics_tab[[col]])
#   }
# }
# remove uninformative columns TODO: include this information in figure header
columns <- names(metrics_tab)[
  sapply(metrics_tab, uniqueN) != 1 |
  names(metrics_tab) %in% c(metrics, 'file_name', 'label', 'batch')
]
metrics_tab <- metrics_tab[, ..columns]

get_col_width <- function(col, dt, factor=0.3) {
  width <- max(nchar(as.character(dt[[col]])), na.rm = TRUE)
  if (is.infinite(width) | width < 1) {
    width <- nchar(col)
  }
  return(factor * width + 1) # +1 for padding
}

#add column info metadata for plotting using funkyheatmap
dt1 <- data.table(
  id=integration_setup,
  group="Integration Setup",
  geom='text',
  palette='setup',
  width=sapply(integration_setup, get_col_width, dt=metrics_tab)
)
dt2 <- data.table(id="Overall Score", group="Overall", geom="bar", palette="overall", width=4)
dt3 <- data.table(id="Bio Conservation", group="Bio Conservation", geom='bar', palette='bio', width=4)
dt4 <- data.table(id=bio_metrics, group="Bio Conservation", geom='funkyrect', palette='bio', width=1)
dt5 <- data.table(id="Batch Correction", group="Batch Correction", geom='bar', palette='batch', width=4)
dt6 <- data.table(id=batch_metrics, group="Batch Correction", geom='circle', palette='batch', width=1)
column_info <- rbind(dt1, dt2, dt3, dt4, dt5, dt6, fill=TRUE)
column_info <- column_info[column_info$id %in% colnames(metrics_tab)]
print("column_info")
print(column_info)

n_top <- min(n_top, nrow(metrics_tab))
g <- funky_heatmap(
  metrics_tab,
  row_info = row_info,
  row_groups = row_groups,
  column_info = column_info,
  column_groups = data.table(
    group = c("Integration Setup", "Overall", "Bio Conservation", "Batch Correction"),
    palette = c('setup', 'overall', 'bio', 'batch'),
    level1 = c("Integration Setup", "Overall", "Bio Conservation", "Batch Correction")
  ),
  palettes = list(
    setup = "Greys",
    overall = "Greens",
    bio = "YlOrBr",
    batch = "Blues"
  ),
  legends = list(
    list(title = "Bio Conservation", palette = "bio", geom = "funkyrect"),
    list(title = "Batch Correction", palette = "batch", geom = "circle")
  ),
  position_args = position_arguments(
    col_annot_offset = 7,
    col_annot_angle = 90,
  ),
  scale_column = FALSE
)
ggsave(output_pdf, g, device = cairo_pdf, width = g$width, height = g$height, dpi=300)
