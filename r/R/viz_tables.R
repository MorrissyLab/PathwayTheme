## Table writers: the 3 PCA summary TSVs plus the differential and category
## tables.  All are written tab-separated without row names, so they line up
## column for column with the Python implementation's output.

.pt_write_tsv <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE,
                     na = "", fileEncoding = "UTF-8")
  path
}

#' Write the three PCA summary tables
#'
#' Variance explained, top pathways per component, and the per-observation
#' pathway signatures.
#'
#' @param result a [pt_pca_result()].
#' @param out_dir directory to write into.
#' @param prefix file-name prefix.
#' @param config a PCA section, from [pt_pca_config()].
#' @return The written paths.
#' @export
pt_write_pca_tables <- function(result, out_dir, prefix, config) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  written <- character()

  # 1. variance explained
  p <- file.path(out_dir, paste0(prefix, "_pca_variance_explained.tsv"))
  written <- c(written, .pt_write_tsv(data.frame(
    PC = pt_pc_names(result),
    variance_explained = result$variance_explained,
    cumulative_variance = cumsum(result$variance_explained),
    stringsAsFactors = FALSE), p))

  # 2. top pathways per PC
  p <- file.path(out_dir, paste0(prefix, "_pca_top_pathways_per_pc.tsv"))
  written <- c(written, .pt_write_tsv(
    pt_top_pathways_per_pc(result$loadings, result$variance_explained,
                           config$top_pathways_per_pc), p))

  # 3. per-cluster signatures
  p <- file.path(out_dir, paste0(prefix, "_pca_per_cluster_signatures.tsv"))
  written <- c(written, .pt_write_tsv(result$signatures, p))

  written
}

#' Write the differential and category tables
#'
#' @param result a [pt_diff_result()].
#' @param summary a category summary from [pt_summarize_by_category()].
#' @param out_dir directory to write into.
#' @param prefix file-name prefix.
#' @return The written path, as a length-one character vector.
#' @export
pt_write_diff_table <- function(result, out_dir, prefix) {
  .pt_write_tsv(result$table,
                file.path(out_dir, paste0(prefix, "_differential_",
                                          result$method, ".tsv")))
}

#' @rdname pt_write_diff_table
#' @export
pt_write_category_table <- function(summary, out_dir, prefix) {
  .pt_write_tsv(summary,
                file.path(out_dir, paste0(prefix, "_category_summary.tsv")))
}
