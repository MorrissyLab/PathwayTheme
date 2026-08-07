## Per-sample gene-list selection shared by list-based backends (EnrichR/GoSlim).
##
## Continuous omic data has no intrinsic "gene list", so for each sample we
## derive one from its column: the top-N features by value, or all features
## above a threshold.  This is the input to over-representation style enrichment.

#' Pick the gene list for one sample column
#'
#' @param col a numeric vector named by feature.
#' @param config an enrichment section, from [pt_enrichment_config()].
#' @return A character vector of feature names.
#' @export
pt_select_genes <- function(col, config) {
  s <- col[!is.na(col)]
  if (identical(config$gene_selection, "threshold")) {
    if (!.pt_is_set(config$threshold))
      stop("gene_selection='threshold' requires enrichment$threshold", call. = FALSE)
    return(as.character(names(s)[s >= config$threshold]))
  }
  n <- min(config$top_n, length(s))
  if (n < 1L) return(character())
  as.character(names(s)[.pt_order_desc(s)[seq_len(n)]])
}
