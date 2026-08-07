## Roll pathway-level results up into broad categories.
##
## Given a long result table (a `pt_diff_result` table, or PCA
## top-pathways/loadings) and a term -> category mapping, aggregate a numeric
## column per category.

#' Read a term -> category mapping
#'
#' @param path a TSV (or CSV, by extension) with a term column and a category
#'   column.
#' @param key_col,category_col the two column names.
#' @param sep field separator; inferred from the extension when `NULL`.
#' @return A named character vector mapping term -> category, first occurrence
#'   winning on duplicates.
#' @export
pt_load_category_map <- function(path, key_col = "pathway",
                                 category_col = "category", sep = NULL) {
  if (is.null(sep))
    sep <- if (identical(tolower(tools::file_ext(path)), "csv")) "," else "\t"
  df <- utils::read.table(path, sep = sep, header = TRUE, check.names = FALSE,
                          quote = "\"", comment.char = "",
                          stringsAsFactors = FALSE)
  for (col in c(key_col, category_col)) {
    if (!col %in% colnames(df))
      stop("column '", col, "' not in ", basename(path), "; found ",
           paste(colnames(df), collapse = ", "), call. = FALSE)
  }
  keys <- as.character(df[[key_col]])
  vals <- as.character(df[[category_col]])
  keep <- !duplicated(keys)
  stats::setNames(vals[keep], keys[keep])
}

#' Aggregate a result table by broad category
#'
#' Returns one row per (group..., \[significance,\] category) with the
#' aggregated value, `n_pathways`, `n_up` and `n_down`.
#'
#' If `significance_col` is given (e.g. `"fdr"`), each category is split into
#' `significant` (value `< alpha`) vs `non_significant` rows, adding a
#' `significance` column.  `NA` thresholds count as non-significant.
#'
#' @param table a long data.frame with a pathway column.
#' @param mapping a named term -> category vector, or a list.
#' @param key_col the pathway column in `table`.
#' @param value_col the numeric column to aggregate.
#' @param group_cols extra columns to aggregate within, e.g. `"comparison"`.
#' @param stat `"mean"`, `"median"` or `"sum"`.
#' @param unmapped_label category assigned to terms absent from `mapping`.
#' @param significance_col optional column to split on.
#' @param alpha significance threshold.
#' @return A data.frame sorted by the grouping keys.
#' @export
pt_summarize_by_category <- function(table, mapping, key_col = "pathway",
                                     value_col = "effect", group_cols = NULL,
                                     stat = "mean", unmapped_label = "Other",
                                     significance_col = NULL, alpha = 0.05) {
  if (!key_col %in% colnames(table))
    stop("key_col '", key_col, "' not in table columns ",
         paste(colnames(table), collapse = ", "), call. = FALSE)
  if (!value_col %in% colnames(table))
    stop("value_col '", value_col, "' not in table columns ",
         paste(colnames(table), collapse = ", "), call. = FALSE)
  if (!stat %in% c("mean", "median", "sum"))
    stop("stat must be mean|median|sum, got '", stat, "'", call. = FALSE)

  m <- unlist(mapping)
  names(m) <- as.character(names(mapping))
  df <- as.data.frame(table, stringsAsFactors = FALSE)
  cat_vals <- unname(m[as.character(df[[key_col]])])
  cat_vals[is.na(cat_vals)] <- unmapped_label
  df$category <- as.character(cat_vals)

  keys <- intersect(.pt_or(group_cols, character()), colnames(df))
  if (!is.null(significance_col)) {
    if (!significance_col %in% colnames(df))
      stop("significance_col '", significance_col, "' not in table columns ",
           paste(colnames(table), collapse = ", "), call. = FALSE)
    is_sig <- !is.na(df[[significance_col]]) & df[[significance_col]] < alpha
    df$significance <- ifelse(is_sig, "significant", "non_significant")
    keys <- c(keys, "significance")
  }
  keys <- c(keys, "category")

  vals <- as.double(df[[value_col]])
  fn <- switch(stat, mean = mean, median = stats::median, sum = sum)
  split_key <- interaction(df[keys], drop = TRUE, sep = "\r", lex.order = TRUE)

  agg <- lapply(split(seq_len(nrow(df)), split_key), function(ix) {
    v <- vals[ix]
    data.frame(stat_value = fn(v, na.rm = TRUE),
               n_pathways = length(ix),
               n_up = sum(v > 0, na.rm = TRUE),
               n_down = sum(v < 0, na.rm = TRUE),
               stringsAsFactors = FALSE)
  })
  key_parts <- do.call(rbind, strsplit(names(agg), "\r", fixed = TRUE))
  out <- cbind(
    stats::setNames(as.data.frame(key_parts, stringsAsFactors = FALSE), keys),
    do.call(rbind, agg))
  names(out)[names(out) == "stat_value"] <- paste0(stat, "_", value_col)
  out <- out[do.call(order, c(unname(as.list(out[keys])),
                              list(method = "radix"))), , drop = FALSE]
  rownames(out) <- NULL
  out
}
