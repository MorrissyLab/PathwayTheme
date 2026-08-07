## Core PCA computation over a term x sample score matrix.
##
##   1. z-score every term (feature) across observations (population sd);
##      drop constants.
##   2. transpose -> observations x terms; observations are the PCA points.
##   3. singular value decomposition, keeping K = min(n_obs - 1, k_max) PCs.
##
## Observations are the columns of the score matrix within one *scope*.  Labels
## / sub-labels / sizes are optional per-observation annotations used only for
## figures and the signature table (they never affect the numeric result).
##
## The component signs follow scikit-learn's `svd_flip(..., u_based_decision =
## FALSE)`: each component is flipped so that its largest-magnitude loading is
## positive.  That makes the scores and loadings directly comparable with the
## Python implementation rather than differing by an arbitrary sign per PC.

#' Run PCA on one scope's score matrix
#'
#' @param score_data a terms x observations numeric matrix.
#' @param labels comparison label per observation, named by observation id.
#' @param config a PCA section, from [pt_pca_config()].
#' @param scope the scope label carried into the result.
#' @param display_labels optional per-observation display names.
#' @param sublabels optional secondary label per observation.
#' @param sizes optional per-observation weight for dot sizing.
#' @param label_name,sublabel_name names of the two label variables, carried
#'   into the result for use as legend titles.
#' @return A [pt_pca_result()], or `NULL` when the scope has too few
#'   observations or too few variable features.
#' @export
pt_run_pca <- function(score_data, labels, config, scope = "__all__",
                       display_labels = NULL, sublabels = NULL, sizes = NULL,
                       label_name = "group", sublabel_name = "secondary") {
  n_na <- sum(is.na(score_data))
  if (n_na > 0L) {
    # a term is NA for an observation only when it was scored for some samples
    # and not others.  0 is a legitimate ssGSEA value, so filling is not
    # neutral -- say so rather than doing it silently.
    warning(sprintf(paste0("PCA scope '%s': %d missing score(s) filled with 0; ",
                           "0 is a valid enrichment score, so these observations ",
                           "are treated as average rather than unknown. Check the ",
                           "gene-set coverage table."), scope, n_na), call. = FALSE)
  }
  x <- score_data
  x[is.na(x)] <- 0
  n_features <- nrow(x)
  n_obs <- ncol(x)
  if (n_obs < config$min_observations || n_features < config$min_features)
    return(NULL)

  obs_ids <- colnames(x)
  disp <- if (!is.null(display_labels))
    as.character(.pt_lookup(display_labels, obs_ids, obs_ids)) else obs_ids

  # ---- z-score per feature across observations; drop constant features ----
  mu <- rowMeans(x)
  sds <- sqrt(rowMeans((x - mu)^2))              # population sd, ddof = 0
  keep <- sds > 1e-12
  if (sum(keep) < config$min_features) return(NULL)
  xn <- (x[keep, , drop = FALSE] - mu[keep]) / sds[keep]

  k <- min(n_obs - 1L, config$k_max)
  fit <- .pt_pca_fit(t(xn), k)

  pc_names <- paste0("PC", seq_len(k))
  scores <- fit$scores
  dimnames(scores) <- list(disp, pc_names)
  loadings <- fit$loadings
  dimnames(loadings) <- list(pc_names, rownames(xn))

  sig <- pt_compute_signatures(scores, loadings,
                               top_pcs_per_obs = config$top_pcs_per_cluster,
                               top_pathways = config$top_pathways_per_cluster)

  labels_disp <- stats::setNames(as.character(.pt_lookup(labels, obs_ids, "")), disp)
  sub_disp <- if (is.null(sublabels)) NULL else
    stats::setNames(as.character(.pt_lookup(sublabels, obs_ids, "")), disp)
  sizes_arr <- if (is.null(sizes)) NULL else
    as.double(.pt_lookup(sizes, obs_ids, 1))

  pt_pca_result(scope = scope, scores = scores, loadings = loadings,
                variance_explained = fit$variance_explained,
                signatures = sig$long, labels = labels_disp,
                sublabels = sub_disp, sizes = sizes_arr,
                labels_name = label_name, sublabels_name = sublabel_name)
}

#' Run one PCA per scope
#'
#' @param sm a [pt_score_matrix()].
#' @param grouping a [pt_grouping()].
#' @param config a PCA section, from [pt_pca_config()].
#' @param display_col,sublabel_col,size_col optional metadata columns providing
#'   per-observation display labels, secondary colour labels and dot sizes.
#' @return A list of [pt_pca_result()] objects, one per scope with enough
#'   observations.
#' @export
pt_run_pca_scopes <- function(sm, grouping, config, display_col = NULL,
                              sublabel_col = NULL, size_col = NULL) {
  meta <- sm$metadata
  col_or_null <- function(nm)
    if (pt_meta_has(meta, nm)) pt_meta_get(meta, nm) else NULL
  disp <- col_or_null(display_col)
  subl <- col_or_null(sublabel_col)
  size <- col_or_null(size_col)

  results <- list()
  scopes <- pt_scopes(grouping)
  for (scope in names(scopes)) {
    cols <- scopes[[scope]]
    cols <- cols[cols %in% colnames(sm$data)]
    if (length(cols) < config$min_observations) next
    res <- pt_run_pca(sm$data[, cols, drop = FALSE], grouping$labels[cols],
                      config, scope = scope, display_labels = disp,
                      sublabels = subl, sizes = size,
                      label_name = .pt_or(grouping$label_name, "group"),
                      sublabel_name = .pt_or(sublabel_col, "secondary"))
    if (!is.null(res)) results[[length(results) + 1L]] <- res
  }
  results
}


## ---- internals -------------------------------------------------------------

## SVD-based PCA matching scikit-learn's PCA(n_components = k).
.pt_pca_fit <- function(x, k) {
  n <- nrow(x)
  centre <- colMeans(x)
  xc <- sweep(x, 2L, centre, "-")
  sv <- svd(xc)
  # sign convention: flip each component so its largest-|loading| entry is
  # positive (scikit-learn's svd_flip with u_based_decision = FALSE)
  signs <- apply(sv$v, 2L, function(col) {
    s <- sign(col[[which.max(abs(col))]])
    if (s == 0) 1 else s
  })
  u <- sweep(sv$u, 2L, signs, "*")
  v <- sweep(sv$v, 2L, signs, "*")

  explained <- sv$d^2 / (n - 1)
  ratio <- explained / sum(explained)

  list(scores = u[, seq_len(k), drop = FALSE] *
         rep(sv$d[seq_len(k)], each = n),
       loadings = t(v[, seq_len(k), drop = FALSE]),
       variance_explained = ratio[seq_len(k)])
}

## Look up `keys` in a named vector, falling back to `default`.
.pt_lookup <- function(src, keys, default) {
  if (is.null(src)) return(rep(default, length(keys)))
  out <- src[keys]
  miss <- is.na(names(out)) | !(keys %in% names(src))
  if (any(miss)) {
    fallback <- if (length(default) == length(keys)) default[miss] else default
    out[miss] <- fallback
  }
  unname(out)
}
