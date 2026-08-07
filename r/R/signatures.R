## Per-observation pathway signatures from a fitted PCA.
##
## For each observation, pick its top-N PCs by |score|, then
##
##     signature(obs, pathway) = sum_{k in top PCs} score[obs, k] * loading[k, pathway]
##
## Rank pathways by |signature|, keep the top M (sign preserved: + = up,
## - = down).

#' Per-observation pathway signatures
#'
#' @param scores observations x PCs matrix.
#' @param loadings PCs x pathways matrix.
#' @param top_pcs_per_obs PCs combined into one observation's signature.
#' @param top_pathways pathways kept per observation.
#' @return A list with `per_obs` (a list of `top_pcs` / `values` per
#'   observation) and `long` (a data.frame with `cluster`, `top_pcs`, `rank`,
#'   `pathway`, `signature_score`, `direction`).
#' @export
pt_compute_signatures <- function(scores, loadings, top_pcs_per_obs,
                                  top_pathways) {
  pc_names <- colnames(scores)
  pathways <- colnames(loadings)
  obs_labels <- rownames(scores)

  per_obs <- vector("list", nrow(scores))
  names(per_obs) <- obs_labels
  records <- vector("list", nrow(scores))

  for (i in seq_len(nrow(scores))) {
    s_vec <- scores[i, ]
    top_pcs <- pc_names[utils::head(.pt_order_desc(abs(s_vec)), top_pcs_per_obs)]
    sig <- as.vector(s_vec[top_pcs] %*% loadings[top_pcs, , drop = FALSE])
    names(sig) <- pathways

    keep <- utils::head(.pt_order_desc(abs(sig)), top_pathways)
    top_paths <- sig[keep]
    per_obs[[i]] <- list(top_pcs = top_pcs, values = top_paths)
    records[[i]] <- data.frame(
      cluster = obs_labels[[i]],
      top_pcs = paste(top_pcs, collapse = ","),
      rank = seq_along(top_paths),
      pathway = names(top_paths),
      signature_score = as.double(top_paths),
      direction = ifelse(top_paths > 0, "up", "down"),
      stringsAsFactors = FALSE)
  }

  long <- if (length(records)) do.call(rbind, records) else
    data.frame(cluster = character(), top_pcs = character(), rank = integer(),
               pathway = character(), signature_score = double(),
               direction = character(), stringsAsFactors = FALSE)
  rownames(long) <- NULL
  list(per_obs = per_obs, long = long)
}

#' Top-loading pathways per principal component
#'
#' @param loadings PCs x pathways matrix.
#' @param variance_explained variance explained per PC.
#' @param top_n pathways tabulated per component.
#' @return A long data.frame with `PC`, `variance_explained`, `rank`,
#'   `pathway`, `loading`, `direction`.
#' @export
pt_top_pathways_per_pc <- function(loadings, variance_explained, top_n) {
  pc_names <- rownames(loadings)
  rows <- vector("list", length(pc_names))
  for (i in seq_along(pc_names)) {
    lv <- loadings[i, ]
    keep <- utils::head(.pt_order_desc(abs(lv)), top_n)
    rows[[i]] <- data.frame(
      PC = pc_names[[i]],
      variance_explained = as.double(variance_explained[[i]]),
      rank = seq_along(keep),
      pathway = colnames(loadings)[keep],
      loading = as.double(lv[keep]),
      direction = ifelse(lv[keep] > 0, "up", "down"),
      stringsAsFactors = FALSE)
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
