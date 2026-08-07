## Input-adapter dispatcher + the generic matrix adapter.
##
## An adapter turns a source (a matrix file, or a directory of .h5ad) into a
## `pt_feature_matrix` -- features x samples plus per-sample metadata.

#' Load input data into a feature matrix
#'
#' Dispatches on `config$kind`: `"matrix"` reads a pre-quantified table,
#' `"h5ad"` aggregates single-cell files into per-cluster pseudobulk.
#'
#' @param config an input section, from [pt_input_config()].
#' @return A [pt_feature_matrix()].
#' @export
pt_load_input <- function(config) {
  kind <- tolower(config$kind)
  if (identical(kind, "matrix")) return(.pt_load_matrix_adapter(config))
  if (identical(kind, "h5ad")) return(.pt_load_h5ad_adapter(config))
  stop("Unknown input kind '", config$kind, "'; expected 'matrix' or 'h5ad'",
       call. = FALSE)
}

.pt_read_table <- function(path, sep, numeric = FALSE) {
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("parquet", "pq", "xlsx", "xls"))
    stop("'.", ext, "' input is not supported by the R implementation; ",
         "export to .tsv or .csv first", call. = FALSE)
  if (identical(ext, "csv")) sep <- ","
  df <- utils::read.table(path, sep = sep, header = TRUE, row.names = 1L,
                          check.names = FALSE, quote = "\"",
                          comment.char = "", stringsAsFactors = FALSE,
                          na.strings = c("NA", "NaN", ""))
  df
}

.pt_to_numeric_matrix <- function(df) {
  # mirrors pandas' to_numeric(errors="coerce"): unparseable entries become NA
  rn <- rownames(df); cn <- colnames(df)
  out <- vapply(df, function(col) {
    if (is.numeric(col)) return(as.double(col))
    suppressWarnings(as.double(as.character(col)))
  }, double(nrow(df)))
  out <- matrix(out, nrow = nrow(df), dimnames = list(rn, cn))
  out
}

.pt_load_matrix_adapter <- function(config) {
  if (!.pt_is_set(config$matrix_path))
    stop("input$matrix_path is required for kind='matrix'", call. = FALSE)
  raw <- .pt_read_table(config$matrix_path, config$sep)
  data <- .pt_to_numeric_matrix(raw)
  rownames(data) <- as.character(rownames(data))
  colnames(data) <- as.character(colnames(data))

  if (.pt_is_set(config$metadata_path)) {
    meta <- .pt_read_table(config$metadata_path, config$sep)
    rownames(meta) <- as.character(rownames(meta))
    missing <- setdiff(colnames(data), rownames(meta))
    if (length(missing))
      stop(length(missing), " matrix columns absent from metadata, e.g. ",
           paste(utils::head(missing, 3), collapse = ", "), call. = FALSE)
    meta <- meta[colnames(data), , drop = FALSE]
  } else {
    meta <- data.frame(sample = colnames(data), row.names = colnames(data),
                       stringsAsFactors = FALSE)
  }
  pt_feature_matrix(data, pt_sample_metadata(meta))
}
