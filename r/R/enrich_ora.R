## EnrichR backend: per-sample over-representation -> term x sample matrix.
##
## For each sample we take a gene list (see [pt_select_genes()]) and run
## enrichment:
##
## * **offline** (default, deterministic, testable): a hypergeometric test of
##   each set in a local `.gmt` (`enrichment$gmt_path`) using the full feature
##   set as the statistical background.  This reproduces `gseapy.enrich`,
##   including its Haldane-Anscombe-corrected odds ratio and the Enrichr
##   combined score `-log(p) * OR`.
## * **online**: the Enrichr web service via the optional `enrichR` package
##   (`enrichment$enrichr_library`) -- requires internet.
##
## The per-sample results are pivoted into `terms x samples` using either
## `-log10(adjusted p)` (default) or the combined score.  Terms absent for a
## sample are 0 (not enriched).

#' Over-representation of a gene list against local gene sets
#'
#' A hypergeometric test per gene set, matching `gseapy.enrich`: the p-value is
#' the upper tail `phyper(x - 1, m, bg - m, k, lower.tail = FALSE)`, the odds
#' ratio carries the Haldane-Anscombe 0.5 correction, and the combined score is
#' `-log(p) * OR`.
#'
#' @param genes the query gene list.
#' @param gene_sets named list of gene sets.
#' @param background character vector of all measurable features.
#' @return A data.frame with `Term`, `Overlap`, `P.value`, `Adjusted.P.value`,
#'   `Odds.Ratio`, `Combined.Score`, one row per set with at least one hit.
#' @export
pt_ora <- function(genes, gene_sets, background) {
  bg_set <- unique(as.character(background))
  query <- intersect(unique(as.character(genes)), bg_set)
  bg <- length(bg_set)
  k <- length(query)
  if (!k) return(.pt_empty_ora())

  terms <- sort(names(gene_sets), method = "radix")
  rows <- vector("list", length(terms))
  for (i in seq_along(terms)) {
    category <- intersect(unique(as.character(gene_sets[[terms[[i]]]])), bg_set)
    m <- length(category)
    x <- length(intersect(query, category))
    if (x < 1L) next
    pval <- stats::phyper(x - 1L, m, bg - m, k, lower.tail = FALSE)
    bu <- 0.5   # Haldane-Anscombe correction, as in gseapy
    oddr <- ((x + bu) * (bg - m - k + x + bu)) / ((m - x + bu) * (k - x + bu))
    rows[[i]] <- data.frame(Term = terms[[i]],
                            Overlap = paste0(x, "/", m),
                            P.value = pval,
                            Odds.Ratio = oddr,
                            Combined.Score = -log(pval) * oddr,
                            stringsAsFactors = FALSE)
  }
  rows <- rows[!vapply(rows, is.null, TRUE)]
  if (!length(rows)) return(.pt_empty_ora())
  out <- do.call(rbind, rows)
  out$Adjusted.P.value <- stats::p.adjust(out$P.value, method = "BH")
  out[, c("Term", "Overlap", "P.value", "Adjusted.P.value",
          "Odds.Ratio", "Combined.Score")]
}

.pt_empty_ora <- function() {
  data.frame(Term = character(), Overlap = character(), P.value = double(),
             Adjusted.P.value = double(), Odds.Ratio = double(),
             Combined.Score = double(), stringsAsFactors = FALSE)
}

## One enrichment result table -> a named term -> score vector.
.pt_ora_score_column <- function(res, how) {
  if (!nrow(res)) return(stats::setNames(numeric(), character()))
  vals <- if (identical(how, "combined_score")) {
    res$Combined.Score
  } else {
    -log10(pmax(res$Adjusted.P.value, 1e-300))
  }
  # duplicate terms collapse to their maximum, matching the Python backend
  tapply(vals, res$Term, max)
}

.pt_backend_enrichr <- function(fm, config, verbose = TRUE) {
  offline <- .pt_is_set(config$gmt_path)
  if (offline) {
    gene_sets <- pt_parse_gmt(config$gmt_path)
  } else if (.pt_is_set(config$enrichr_library)) {
    if (!requireNamespace("enrichR", quietly = TRUE))
      stop("online Enrichr needs the 'enrichR' package; supply ",
           "enrichment$gmt_path to run offline instead", call. = FALSE)
    gene_sets <- NULL
  } else {
    stop("EnrichR needs enrichment$gmt_path (offline) or ",
         "enrichment$enrichr_library (online)", call. = FALSE)
  }

  background <- pt_features(fm)
  per_sample <- list()
  for (sample in pt_samples(fm)) {
    col <- fm$data[, sample]
    names(col) <- rownames(fm$data)
    genes <- pt_select_genes(col, config)
    if (!length(genes)) next
    res <- tryCatch({
      if (offline) pt_ora(genes, gene_sets, background)
      else .pt_enrichr_online(genes, config)
    }, error = function(e) {
      if (verbose) message("  [enrichr] ", sample, ": ", conditionMessage(e))
      NULL
    })
    if (is.null(res) || !nrow(res)) next
    per_sample[[sample]] <- .pt_ora_score_column(res, config$enrichr_score)
  }
  if (!length(per_sample))
    stop("EnrichR produced no results for any sample", call. = FALSE)
  .pt_bind_score_columns(per_sample, pt_samples(fm))
}

.pt_enrichr_online <- function(genes, config) {
  res <- enrichR::enrichr(genes, databases = config$enrichr_library)
  df <- res[[1L]]
  data.frame(Term = as.character(df$Term),
             Overlap = as.character(df$Overlap),
             P.value = as.double(df$P.value),
             Adjusted.P.value = as.double(df$Adjusted.P.value),
             Odds.Ratio = as.double(df$Odds.Ratio),
             Combined.Score = as.double(df$Combined.Score),
             stringsAsFactors = FALSE)
}

## Align a list of named term -> score vectors into a terms x samples matrix,
## zero-filling terms a sample did not hit.
.pt_bind_score_columns <- function(per_sample, sample_order) {
  cols <- sample_order[sample_order %in% names(per_sample)]
  terms <- sort(unique(unlist(lapply(per_sample, names), use.names = FALSE)),
                method = "radix")
  out <- matrix(0, nrow = length(terms), ncol = length(cols),
                dimnames = list(terms, cols))
  for (j in seq_along(cols)) {
    v <- per_sample[[cols[[j]]]]
    out[names(v), j] <- as.double(v)
  }
  out
}
