## ssGSEA backend (continuous per-sample NES matrix).
##
## A direct implementation of the ssGSEAProjection enrichment score used by
## GenePattern / ssGSEA2.0 and by gseapy's `ssgsea(sample_norm_method="rank",
## correl_norm_type="rank")`, which is what the Python implementation calls:
##
##   1. within each sample, replace values by their average ranks, scaled to
##      `10000 * rank / n_features`;
##   2. order genes by decreasing rank metric;
##   3. walk the ordered list accumulating `|r|^weight / sum(|r|^weight)` at set
##      members and `-1/(N - Nh)` elsewhere;
##   4. ES is the area under that running sum;
##   5. NES rescales every ES by the range of ES over the whole result table.
##
## Step 3-4 are evaluated in closed form.  Writing the running sum as
## `RES_i = sum_{j<=i} t_j` gives `sum_i RES_i = sum_j (N - j + 1) * t_j`, so the
## area needs only the positions of the set's own members -- O(set size) per
## sample instead of O(n_features).

#' Rank-normalise a matrix column-wise, the way ssGSEA does
#'
#' Average ranks within each sample, missing values ranked last (sharing the
#' average of the trailing positions), scaled to `10000 * rank / n_features`.
#'
#' @param data numeric matrix, features x samples.
#' @return A matrix of the same shape.
#' @export
pt_rank_normalise <- function(data) {
  n <- nrow(data)
  out <- apply(data, 2L, function(col) {
    na <- is.na(col)
    r <- numeric(n)
    n_valid <- sum(!na)
    if (n_valid) r[!na] <- rank(col[!na], ties.method = "average")
    # every NA is tied at the bottom: they share the mean trailing position
    if (any(na)) r[na] <- mean(seq.int(n_valid + 1L, n))
    r
  })
  out <- matrix(out, nrow = n, dimnames = dimnames(data))
  out * 10000 / n
}

## Enrichment score (area under the running sum) for one ordered sample.
## `hit_pos` are the 1-based positions of the set's members in the ordered list.
.pt_ssgsea_es <- function(hit_pos, cv_sorted, n) {
  n_h <- length(hit_pos)
  n_m <- n - n_h
  if (n_h == 0L || n_m == 0L) return(NA_real_)
  w <- n - hit_pos + 1                     # position weights of the hits
  cv <- cv_sorted[hit_pos]
  s <- sum(cv)
  if (!is.finite(s) || s <= 0) return(NA_real_)
  total_w <- n * (n + 1) / 2
  sum(w * cv) / s - (total_w - sum(w)) / n_m
}

#' Single-sample GSEA scores for a matrix
#'
#' @param data numeric matrix, features x samples (raw values, not yet ranked).
#' @param gene_sets named list of gene sets.
#' @param weight exponent applied to the rank metric.
#' @param min_size,max_size bounds on the number of set members present in
#'   `data`; sets outside the range are skipped.
#' @param normalise divide every score by the range of scores over the whole
#'   table, reproducing gseapy's `NES` column.  Set `FALSE` to get raw `ES`.
#' @return A terms x samples numeric matrix, terms in sorted order.
#' @export
pt_ssgsea <- function(data, gene_sets, weight = 0.25, min_size = 5L,
                      max_size = 5000L, normalise = TRUE) {
  n <- nrow(data)
  features <- rownames(data)

  # restrict every set to measurable genes, then apply the size bounds
  members <- lapply(gene_sets, function(g) unique(match(g, features)))
  members <- lapply(members, function(ix) ix[!is.na(ix)])
  sizes <- lengths(members)
  keep <- sizes >= min_size & sizes <= max_size
  members <- members[keep]
  if (!length(members))
    stop("no gene set has between ", min_size, " and ", max_size,
         " members present in the feature matrix", call. = FALSE)
  members <- members[sort(names(members), method = "radix")]

  ranked <- pt_rank_normalise(data)
  out <- matrix(NA_real_, nrow = length(members), ncol = ncol(data),
                dimnames = list(names(members), colnames(data)))

  for (j in seq_len(ncol(ranked))) {
    v <- ranked[, j]
    ord <- .pt_order_desc(v)
    cv <- abs(v[ord])^weight
    # position of each feature in the ordered list
    pos <- integer(n)
    pos[ord] <- seq_len(n)
    for (i in seq_along(members)) {
      out[i, j] <- .pt_ssgsea_es(pos[members[[i]]], cv, n)
    }
  }

  if (normalise) {
    rng <- suppressWarnings(max(out, na.rm = TRUE) - min(out, na.rm = TRUE))
    if (is.finite(rng) && rng > 0) out <- out / rng
  }
  out
}

.pt_backend_ssgsea <- function(fm, config, verbose = TRUE) {
  if (!.pt_is_set(config$gmt_path))
    stop("ssGSEA requires enrichment$gmt_path (a .gmt file)", call. = FALSE)
  gs <- pt_parse_gmt(config$gmt_path)

  expr <- fm$data
  harmonised <- pt_harmonize_case(rownames(expr), gs)
  rownames(expr) <- harmonised$transform(rownames(expr))
  gs <- harmonised$gene_sets
  expr <- .pt_collapse_duplicate_rows(expr)

  scores <- pt_ssgsea(expr, gs, weight = config$weight,
                      min_size = config$min_gene_set_size,
                      max_size = config$max_gene_set_size)
  keep <- colnames(scores) %in% pt_samples(fm)
  scores[, colnames(scores)[keep], drop = FALSE]
}

## Average rows that share a name, as pandas' groupby(level=0).mean() does.
.pt_collapse_duplicate_rows <- function(m) {
  if (!anyDuplicated(rownames(m))) return(m)
  rn <- rownames(m)
  out <- rowsum(m, group = rn, reorder = TRUE, na.rm = FALSE)
  counts <- as.vector(table(rn)[rownames(out)])
  out / counts
}
