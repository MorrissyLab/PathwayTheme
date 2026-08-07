## Data contracts that flow between pipeline stages.
##
## The whole package is glued together by a small number of deliberately boring
## structures: a numeric matrix plus a metadata data.frame.  Keeping these thin
## means every stage is independently testable and swappable.
##
## Conventions
## -----------
## * A **sample** is an observation column: a bulk sample, a pseudobulk
##   `sample|cluster` column, or any other unit the user chose.
## * `pt_feature_matrix` : rows = features (genes), cols = samples.
## * `pt_score_matrix`   : rows = terms (pathways), cols = samples.  Every
##   enrichment backend returns this identical shape.
## * `pt_sample_metadata`: one row per matrix column, row names = sample ids.


# ---- sample metadata -------------------------------------------------------

#' Per-sample attributes
#'
#' A thin wrapper over a `data.frame` whose row names are the matrix column
#' labels.  Row names are coerced to character so sample ids never silently
#' become integers.
#'
#' @param table a `data.frame` with one row per sample.
#' @return An object of class `pt_sample_metadata`.
#' @export
pt_sample_metadata <- function(table) {
  if (!is.data.frame(table)) table <- as.data.frame(table, stringsAsFactors = FALSE)
  rownames(table) <- as.character(rownames(table))
  structure(list(table = table), class = "pt_sample_metadata")
}

#' @export
print.pt_sample_metadata <- function(x, ...) {
  cat("<pt_sample_metadata> ", nrow(x$table), " sample(s) x ",
      ncol(x$table), " attribute(s)\n", sep = "")
  cat("  columns: ", paste(colnames(x$table), collapse = ", "), "\n", sep = "")
  invisible(x)
}

#' Restrict and reorder metadata to a set of samples
#'
#' @param meta a [pt_sample_metadata()].
#' @param columns sample ids, in the order wanted.
#' @return A `pt_sample_metadata` holding exactly `columns`, in that order.
#' @export
pt_align_to <- function(meta, columns) {
  cols <- as.character(columns)
  missing <- setdiff(cols, rownames(meta$table))
  if (length(missing))
    stop(length(missing), " sample(s) missing from metadata, e.g. ",
         paste(utils::head(missing, 3), collapse = ", "), call. = FALSE)
  pt_sample_metadata(meta$table[cols, , drop = FALSE])
}

#' Look up one metadata column
#'
#' @param meta a [pt_sample_metadata()].
#' @param column column name.
#' @return `pt_meta_get()` returns the column as a vector named by sample id;
#'   `pt_meta_has()` returns a logical scalar.
#' @export
pt_meta_get <- function(meta, column) {
  if (!column %in% colnames(meta$table))
    stop("metadata column '", column, "' not found; available: ",
         paste(colnames(meta$table), collapse = ", "), call. = FALSE)
  out <- meta$table[[column]]
  names(out) <- rownames(meta$table)
  out
}

#' @rdname pt_meta_get
#' @export
pt_meta_has <- function(meta, column) {
  !is.null(column) && length(column) == 1L && column %in% colnames(meta$table)
}


# ---- feature matrix --------------------------------------------------------

#' Features (genes) x samples matrix
#'
#' The input to enrichment.  Row and column names are coerced to character.
#'
#' @param data a numeric matrix or data.frame, features in rows.
#' @param metadata a [pt_sample_metadata()] (or a data.frame, which is wrapped).
#' @return An object of class `pt_feature_matrix`.
#' @export
pt_feature_matrix <- function(data, metadata) {
  data <- .pt_as_matrix(data)
  if (!inherits(metadata, "pt_sample_metadata")) metadata <- pt_sample_metadata(metadata)
  structure(list(data = data, metadata = metadata), class = "pt_feature_matrix")
}

#' @export
print.pt_feature_matrix <- function(x, ...) {
  cat("<pt_feature_matrix> ", nrow(x$data), " features x ", ncol(x$data),
      " samples\n", sep = "")
  invisible(x)
}


# ---- score matrix ----------------------------------------------------------

#' Terms (pathways) x samples score matrix
#'
#' The unified output of every enrichment backend.
#'
#' @param data a numeric matrix or data.frame, terms in rows.
#' @param metadata a [pt_sample_metadata()] (or a data.frame).
#' @param backend name of the backend that produced the scores.
#' @param geneset label used for output naming.
#' @param coverage optional gene-set coverage table (see [pt_geneset_coverage()]).
#' @return An object of class `pt_score_matrix`.
#' @export
pt_score_matrix <- function(data, metadata, backend = "unknown",
                            geneset = "unknown", coverage = NULL) {
  data <- .pt_as_matrix(data)
  if (!inherits(metadata, "pt_sample_metadata")) metadata <- pt_sample_metadata(metadata)
  structure(list(data = data, metadata = metadata, backend = backend,
                 geneset = geneset, coverage = coverage),
            class = "pt_score_matrix")
}

#' @export
print.pt_score_matrix <- function(x, ...) {
  cat("<pt_score_matrix> ", nrow(x$data), " terms x ", ncol(x$data),
      " samples  (backend=", x$backend, ", geneset=", x$geneset, ")\n", sep = "")
  invisible(x)
}

#' Restrict a score matrix to a set of samples
#'
#' @param sm a [pt_score_matrix()].
#' @param columns sample ids to keep, in the order wanted.
#' @return A `pt_score_matrix`.
#' @export
pt_subset_samples <- function(sm, columns) {
  cols <- as.character(columns)
  pt_score_matrix(sm$data[, cols, drop = FALSE], pt_align_to(sm$metadata, cols),
                  backend = sm$backend, geneset = sm$geneset,
                  coverage = sm$coverage)
}


# ---- shared accessors ------------------------------------------------------

#' Row and column labels of a matrix contract
#'
#' @param x a [pt_feature_matrix()] or [pt_score_matrix()].
#' @return A character vector.
#' @export
pt_samples <- function(x) colnames(x$data)

#' @rdname pt_samples
#' @export
pt_features <- function(x) rownames(x$data)

#' @rdname pt_samples
#' @export
pt_terms <- function(x) rownames(x$data)


# ---- grouping --------------------------------------------------------------

#' Comparison labels and optional PCA scope
#'
#' `labels` drives figure colouring and comparison.  `scope` optionally
#' partitions samples into independent PCA runs (e.g. one PCA per biological
#' sample); when `scope` is `NULL` the whole matrix is a single scope.
#'
#' @param labels character vector named by sample id.
#' @param scope optional character vector named by sample id.
#' @param label_name name of the label variable (used as a legend title).
#' @param scope_name name of the scope variable.
#' @return An object of class `pt_grouping`.
#' @export
pt_grouping <- function(labels, scope = NULL, label_name = "group",
                        scope_name = NULL) {
  labels <- .pt_named_chr(labels)
  if (!is.null(scope)) scope <- .pt_named_chr(scope)
  structure(list(labels = labels, scope = scope, label_name = label_name,
                 scope_name = scope_name), class = "pt_grouping")
}

#' @export
print.pt_grouping <- function(x, ...) {
  cat("<pt_grouping> ", length(x$labels), " sample(s), ",
      length(unique(x$labels)), " label(s) of '", x$label_name, "', ",
      length(pt_scopes(x)), " scope(s)\n", sep = "")
  invisible(x)
}

#' Samples grouped by PCA scope
#'
#' @param grouping a [pt_grouping()].
#' @return A named list mapping each scope value to its sample ids.  A grouping
#'   without a scope yields a single `"__all__"` entry.
#' @export
pt_scopes <- function(grouping) {
  if (is.null(grouping$scope)) return(list(`__all__` = names(grouping$labels)))
  sc <- as.character(grouping$scope)
  samples <- names(grouping$scope)
  # preserve first-appearance order, matching the Python implementation
  out <- split(samples, factor(sc, levels = unique(sc)))
  lapply(out, as.character)
}


# ---- differential result ---------------------------------------------------

#' Differential pathway analysis result
#'
#' `table` is long -- one row per (comparison, pathway) with columns
#' `comparison`, `pathway`, `method`, `n_case`, `n_reference`, `mean_case`,
#' `mean_reference`, `effect` (= `mean_case - mean_reference`), `statistic`,
#' `p_value`, `fdr`, `direction`.  Several comparisons are stacked in one table.
#'
#' @param table the long result data.frame.
#' @param method the test that produced it.
#' @param group_col name of the variable that defined the groups.
#' @return An object of class `pt_diff_result`.
#' @export
pt_diff_result <- function(table, method = "welch", group_col = "group") {
  structure(list(table = table, method = method, group_col = group_col),
            class = "pt_diff_result")
}

#' @export
print.pt_diff_result <- function(x, ...) {
  n_sig <- sum(x$table$fdr < 0.05, na.rm = TRUE)
  cat("<pt_diff_result> method=", x$method, ", ", length(pt_comparisons(x)),
      " comparison(s), ", nrow(x$table), " row(s), ", n_sig,
      " with FDR<0.05\n", sep = "")
  invisible(x)
}

#' Comparisons and significant rows of a differential result
#'
#' @param result a [pt_diff_result()].
#' @param alpha FDR threshold.
#' @return `pt_comparisons()` returns the comparison names in table order;
#'   `pt_significant()` returns the rows below `alpha`.
#' @export
pt_comparisons <- function(result) {
  if (!"comparison" %in% names(result$table)) return(character())
  unique(as.character(result$table$comparison))
}

#' @rdname pt_comparisons
#' @export
pt_significant <- function(result, alpha = 0.05) {
  tab <- result$table
  tab[!is.na(tab$fdr) & tab$fdr < alpha, , drop = FALSE]
}


# ---- PCA result ------------------------------------------------------------

#' Output of one PCA run over a single scope
#'
#' @param scope the scope label this result belongs to.
#' @param scores observations x PCs matrix (row names are display labels).
#' @param loadings PCs x features matrix.
#' @param variance_explained numeric vector of length K.
#' @param signatures long data.frame: `cluster`, `top_pcs`, `rank`, `pathway`,
#'   `signature_score`, `direction`.
#' @param labels comparison label per observation, aligned to `rownames(scores)`.
#' @param sublabels optional secondary label per observation.
#' @param sizes optional per-observation weight used for dot sizing.
#' @param labels_name,sublabels_name names of the two label variables, used as
#'   legend titles.
#' @return An object of class `pt_pca_result`.
#' @export
pt_pca_result <- function(scope, scores, loadings, variance_explained,
                          signatures, labels, sublabels = NULL, sizes = NULL,
                          labels_name = "group", sublabels_name = "secondary") {
  structure(list(scope = scope, scores = scores, loadings = loadings,
                 variance_explained = variance_explained,
                 signatures = signatures, labels = labels,
                 sublabels = sublabels, sizes = sizes,
                 labels_name = labels_name, sublabels_name = sublabels_name),
            class = "pt_pca_result")
}

#' @export
print.pt_pca_result <- function(x, ...) {
  cat("<pt_pca_result> scope='", x$scope, "'  ", nrow(x$scores),
      " observations x ", ncol(x$scores), " PCs over ", ncol(x$loadings),
      " pathways\n", sep = "")
  invisible(x)
}

#' Principal component names and count of a PCA result
#'
#' @param result a [pt_pca_result()].
#' @return `pt_pc_names()` a character vector; `pt_k()` an integer.
#' @export
pt_pc_names <- function(result) colnames(result$scores)

#' @rdname pt_pc_names
#' @export
pt_k <- function(result) ncol(result$scores)


# ---- internal --------------------------------------------------------------

.pt_as_matrix <- function(data) {
  if (is.data.frame(data)) {
    rn <- rownames(data)
    data <- as.matrix(data)
    rownames(data) <- rn
  }
  if (!is.matrix(data)) stop("expected a matrix or data.frame", call. = FALSE)
  storage.mode(data) <- "double"
  if (is.null(rownames(data))) rownames(data) <- as.character(seq_len(nrow(data)))
  if (is.null(colnames(data))) colnames(data) <- as.character(seq_len(ncol(data)))
  rownames(data) <- as.character(rownames(data))
  colnames(data) <- as.character(colnames(data))
  data
}

## Descending order, with ties broken the way NumPy and pandas break them.
##
## Both sort descending by taking a *stable ascending* order and reversing it,
## so tied values come out in reverse of their original order -- not in original
## order, which is what `order(decreasing = TRUE)` would give.  With rank-tied
## data (a matrix with many exact zeros, say) the two conventions produce
## visibly different ssGSEA scores, so this has to match.
.pt_order_desc <- function(x) rev(order(x, method = "radix"))

.pt_named_chr <- function(x) {
  nm <- names(x)
  out <- as.character(x)
  names(out) <- as.character(nm)
  out
}
