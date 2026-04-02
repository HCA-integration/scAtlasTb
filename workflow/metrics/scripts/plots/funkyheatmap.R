input_file <- snakemake@input$tsv
output_tsv <- snakemake@output$tsv
output_pdf <- snakemake@output$pdf
output_pdf_overall <- snakemake@output$pdf_overall

id_vars <- snakemake@params$id_vars
value_var <- snakemake@params$value_var
variable_var <- snakemake@params$variable_var
weight_batch <- snakemake@params$weight_batch
n_top <- snakemake@params$n_top
group_col <- snakemake@params$group_col
scale <- snakemake@params$scale
scale <- ifelse(is.null(scale), FALSE, scale)
dpi <- snakemake@params$dpi
dpi <- ifelse(is.null(dpi), 300, dpi)

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
  library(RColorBrewer)
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
metrics <- c(batch_metrics, bio_metrics)

# subset data.table to columns of interest & transform
all_columns <- c(id_vars, variable_var, value_var, 'metric_type')
id_formula <- paste(paste(id_vars, collapse='+'), variable_var, sep = '~')
cat(paste(id_formula, '\n'))
metrics_tab <- dcast(
  subset(dt, select = all_columns),
  as.formula(id_formula),
  value.var = value_var,
  fun.aggregate = function(x) mean(x, na.rm = TRUE)
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
  metrics_tab <- metrics_tab[order(-`Overall Score`, na.last = TRUE)]
}

metric_columns <- c(
  "Overall Score",
  "Batch Correction",
  batch_metrics,
  "Bio Conservation",
  bio_metrics
)
metric_columns <- metric_columns[metric_columns %in% colnames(metrics_tab)]
non_metric_columns <- setdiff(colnames(metrics_tab), metric_columns)
metrics_tab <- metrics_tab[, c(non_metric_columns, metric_columns), with = FALSE]

# write output summary.csv file before plotting ? Michaela to check the output dir?
cat('Writing output file: ', output_tsv, '\n')
fwrite(metrics_tab, output_tsv, sep='\t')

# subset data
metrics_tab <- metrics_tab[1:min(n_top, nrow(metrics_tab))]

# Check if we have enough data to create a meaningful heatmap
if (nrow(metrics_tab) < 2) {
  cat("Warning: Only", nrow(metrics_tab), "row(s) available. Skipping funkyheatmap generation.\n")
  cat("Creating empty PDF placeholder...\n")
  pdf(output_pdf, width = 8, height = 4)
  plot.new()
  text(
    0.5, 0.5,
    paste(
      "Insufficient data for heatmap visualization\n",
      "Only", nrow(metrics_tab), "method(s) available\n",
      "See", basename(output_tsv), "for results"
    ),
    cex = 1.2, col = "gray40"
  )
  dev.off()
  file.copy(output_pdf, output_pdf_overall, overwrite = TRUE)
  quit(save = "no", status = 0)
}

### add funkyheatmap data

row_info <- NULL
row_groups <- NULL

# Group rows by the specified column if provided
if (!is.null(group_col) && "Overall Score" %in% colnames(metrics_tab)) {
  # Verify that the grouping column exists in the table
  assertthat::assert_that(group_col %in% colnames(metrics_tab))

  # Rank groups by best performer, then second-best performer.
  metrics_tab[, .group_value := as.character(get(group_col))]
  metrics_tab[is.na(.group_value), .group_value := "NA"]

  metrics_tab[, c("best_score", "second_best_score") := {
    s <- sort(`Overall Score`[!is.na(`Overall Score`)], decreasing = TRUE)
    list(
      if (length(s) >= 1) s[1] else -Inf,
      if (length(s) >= 2) s[2] else -Inf
    )
  }, by = .group_value]

  metrics_tab <- metrics_tab[
    order(
      -best_score,
      -second_best_score,
      .group_value,
      -`Overall Score`,
      na.last = TRUE
    )
  ]

  # Build row grouping metadata before removing the grouping column
  metrics_tab[, id := as.character(seq_len(.N))]
  row_info <- data.table(
    id = metrics_tab$id,
    group = paste(group_col, metrics_tab[[group_col]], sep = '=')
  )
  row_groups <- data.table(
    group = unique(row_info$group),
    Group = unique(row_info$group)
  )

  # remove grouping column and temporary score columns
  metrics_tab[, c(group_col, "best_score", "second_best_score", ".group_value") := NULL]
}

# remove uninformative columns
columns_to_keep <- c(metrics, 'file_name', 'label', 'batch', 'Bio Conservation', 'Batch Correction', 'Overall Score')
columns <- names(metrics_tab)[
  sapply(metrics_tab, uniqueN) != 1 |
  names(metrics_tab) %in% columns_to_keep
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
dt2 <- data.table(
  id="Overall Score",
  group="Overall",
  geom="bar",
  palette="overall",
  width=4
)
dt3 <- data.table(
  id="Batch Correction",
  group="Batch Correction",
  geom='bar',
  palette='batch',
  width=4
)
dt4 <- data.table(
  id=batch_metrics,
  name=batch_metrics,
  group="Batch Correction",
  geom='circle',
  palette='batch',
  width=1
)
dt5 <- data.table(
  id="Bio Conservation",
  group="Bio Conservation",
  geom='bar',
  palette='bio',
  width=4
)
dt6 <- data.table(
  id=bio_metrics,
  name=bio_metrics,
  group="Bio Conservation",
  geom='funkyrect',
  palette='bio',
  width=1
)
column_info <- rbind(dt1, dt2, dt3, dt4, dt5, dt6, fill=TRUE)
column_info <- column_info[column_info$id %in% colnames(metrics_tab)]
column_info$name <- column_info$id
print("column_info")
print(column_info)

n_top <- min(n_top, nrow(metrics_tab))

# Define shared plot configuration
plot_column_groups <- data.table(
  group = c("Integration Setup", "Overall", "Batch Correction", "Bio Conservation"),
  palette = c('setup', 'overall', 'batch', 'bio'),
  level1 = c("Integration Setup", "Overall", "Batch Correction", "Bio Conservation")
)

plot_palettes <- list(
  setup = rev(brewer.pal(5, "Greys")),
  overall = rev(brewer.pal(5, "Greens")),
  bio = rev(brewer.pal(5, "YlOrBr")),
  batch = rev(brewer.pal(5, "Blues"))
)

# Filter palettes to only include those referenced in column_info
used_palettes <- unique(column_info$palette[!is.na(column_info$palette)])
plot_palettes <- plot_palettes[names(plot_palettes) %in% used_palettes]

plot_legends <- list(
  list(title = "Batch Correction", palette = "batch", geom = "circle"),
  list(title = "Bio Conservation", palette = "bio", geom = "funkyrect")
)

# Adapt label spacing to annotation length to avoid overlaps.
max_col_name_len <- max(nchar(as.character(column_info$name)), na.rm = TRUE)
has_row_groups <- !is.null(row_groups) && nrow(row_groups) > 0
plot_position_args <- position_arguments(
  col_annot_offset = max(3, min(10, ceiling(max_col_name_len * 0.4))),
  col_annot_angle = 90,
)

# Extract dataset name if available
title_text <- NULL
if ('dataset' %in% colnames(metrics_tab)) {
  datasets <- unique(na.omit(metrics_tab$dataset))
  if (length(datasets) > 0) {
    title_text <- paste("Dataset:", paste(datasets, collapse=", "))
    cat('Title text:', title_text, '\n')
  }
}

# Create full heatmap
g <- funky_heatmap(
  metrics_tab,
  row_info = row_info,
  row_groups = row_groups,
  column_info = column_info,
  column_groups = plot_column_groups,
  palettes = plot_palettes,
  legends = plot_legends,
  position_args = plot_position_args,
  scale_column = scale
)

if (!is.null(title_text)) g <- g + ggtitle(title_text)
ggsave(
  output_pdf,
  g,
  device = cairo_pdf,
  width = g$width,
  height = g$height,
  dpi=dpi
)

# Create overall metrics heatmap (subset with only overall/bio/batch scores)
overall_cols <- c(integration_setup, "Overall Score", "Batch Correction", "Bio Conservation")
overall_cols <- overall_cols[overall_cols %in% colnames(metrics_tab)]
g <- funky_heatmap(
  metrics_tab[, ..overall_cols],
  row_info = row_info,
  row_groups = row_groups,
  column_info = column_info[id %in% overall_cols],
  palettes = plot_palettes,
  legends = plot_legends,
  position_args = plot_position_args,
  scale_column = scale
)

if (!is.null(title_text)) g <- g + ggtitle(title_text)
ggsave(
  output_pdf_overall,
  g,
  width = g$width,
  height = g$height,
  device = cairo_pdf,
  dpi=dpi
)
