## Resolve a comparison grouping from one of three sources.
##
## * `target`   : colour/compare by a named metadata column (e.g. classification)
## * `existing` : use cluster labels already present in the metadata column
## * `auto`     : cluster the score matrix when no labels exist
##
## An optional `scope_col` partitions samples into independent PCA runs (e.g.
## one PCA per biological sample).  Without it, all samples form a single scope.

#' Resolve a comparison grouping
#'
#' @param sm a [pt_score_matrix()].
#' @param config a grouping section, from [pt_grouping_config()].
#' @return A [pt_grouping()].
#' @export
pt_resolve_grouping <- function(sm, config) {
  meta <- sm$metadata$table
  samples <- colnames(sm$data)

  mode <- tolower(config$mode)
  if (mode %in% c("target", "existing")) {
    if (!.pt_is_set(config$target_col))
      stop("grouping$mode='", mode, "' requires grouping$target_col", call. = FALSE)
    if (!config$target_col %in% colnames(meta))
      stop("target_col '", config$target_col, "' not in metadata columns ",
           paste(colnames(meta), collapse = ", "), call. = FALSE)
    labels <- as.character(meta[samples, config$target_col])
    names(labels) <- samples
    label_name <- config$target_col
  } else if (mode == "auto") {
    labels <- pt_auto_cluster(sm$data, config)
    labels <- labels[samples]
    names(labels) <- samples
    label_name <- "auto_cluster"
  } else {
    stop("Unknown grouping$mode '", config$mode,
         "'; expected 'target' | 'existing' | 'auto'", call. = FALSE)
  }

  scope <- NULL
  scope_name <- NULL
  if (.pt_is_set(config$scope_col)) {
    if (!config$scope_col %in% colnames(meta))
      stop("scope_col '", config$scope_col, "' not in metadata", call. = FALSE)
    scope <- as.character(meta[samples, config$scope_col])
    names(scope) <- samples
    scope_name <- config$scope_col
  }

  pt_grouping(labels = labels, scope = scope, label_name = label_name,
              scope_name = scope_name)
}


## ---- auto clustering -------------------------------------------------------

#' Cluster samples in term space
#'
#' Observations are the score-matrix columns (samples); features are the terms.
#' Features are z-scored across samples, then k-means or Ward agglomerative
#' clustering is run.  `k` is chosen by mean silhouette width over
#' `k_min:k_max` unless a fixed `n_clusters` is given.
#'
#' Cluster *identities* are arbitrary and will not match the Python
#' implementation label for label -- the partitions are what carry meaning.
#'
#' @param score_data a terms x samples numeric matrix.
#' @param config a grouping section, from [pt_grouping_config()].
#' @return A character vector of labels (`"c0"`, `"c1"`, ...) named by sample.
#' @export
pt_auto_cluster <- function(score_data, config) {
  samples <- colnames(score_data)
  n <- length(samples)
  if (n < 2L) return(stats::setNames(rep("c0", n), samples))

  x <- t(score_data)
  x[is.na(x)] <- 0
  x <- .pt_standard_scale(x)

  fit <- function(k) {
    k <- max(1L, min(as.integer(k), n))
    if (identical(config$method, "hierarchical")) {
      stats::cutree(stats::hclust(stats::dist(x), method = "ward.D2"), k = k)
    } else {
      stats::kmeans(x, centers = k, nstart = 10L)$cluster
    }
  }

  if (.pt_is_set(config$n_clusters)) {
    labels <- fit(min(config$n_clusters, n))
  } else {
    k_lo <- max(2L, config$k_min)
    k_hi <- min(config$k_max, n - 1L)
    best_score <- -1
    labels <- NULL
    if (k_hi >= k_lo) {
      for (k in k_lo:k_hi) {
        lab <- fit(k)
        if (length(unique(lab)) < 2L) next
        s <- pt_silhouette(x, lab)
        if (is.finite(s) && s > best_score) {
          best_score <- s
          labels <- lab
        }
      }
    }
    if (is.null(labels)) labels <- fit(min(2L, n))
  }
  stats::setNames(paste0("c", as.integer(labels) - 1L), samples)
}

## z-score every column; zero-variance columns are left at zero, as
## sklearn's StandardScaler does
.pt_standard_scale <- function(x) {
  mu <- colMeans(x)
  sdv <- sqrt(colMeans((x - rep(mu, each = nrow(x)))^2))
  sdv[sdv < .Machine$double.eps] <- 1
  sweep(sweep(x, 2L, mu, "-"), 2L, sdv, "/")
}

#' Mean silhouette width of a clustering
#'
#' @param x an observations x features numeric matrix.
#' @param labels cluster assignment per observation.
#' @return The mean silhouette width, or `NA` when fewer than two clusters.
#' @export
pt_silhouette <- function(x, labels) {
  labels <- as.character(labels)
  clusters <- unique(labels)
  if (length(clusters) < 2L) return(NA_real_)
  d <- as.matrix(stats::dist(x))
  n <- nrow(d)
  sil <- numeric(n)
  for (i in seq_len(n)) {
    own <- labels[[i]]
    same <- which(labels == own)
    same <- same[same != i]
    if (!length(same)) { sil[[i]] <- 0; next }
    a <- mean(d[i, same])
    b <- min(vapply(setdiff(clusters, own),
                    function(cl) mean(d[i, labels == cl]), 1))
    sil[[i]] <- if (max(a, b) > 0) (b - a) / max(a, b) else 0
  }
  mean(sil)
}
