#' pathwaytheme: pathway-level summaries of any omic data
#'
#' The pipeline is five stages, each one function, each output feeding the next:
#'
#' \preformatted{
#' fm      <- pt_load_matrix("expr.tsv", metadata = "meta.tsv")
#' scores  <- pt_enrich(fm, backend = "ssgsea", gmt = "go_bp.gmt", geneset = "GO_BP")
#' groups  <- pt_group(scores, mode = "target", target_col = "treatment")
#' results <- pt_pca(scores, groups)
#' pt_figures(results, "out", geneset = "GO_BP")
#' }
#'
#' or all at once with [pt_run_pipeline()].
#'
#' This package is the R implementation of PathwayTheme and mirrors the Python
#' package stage for stage.  See `vignette`-free docs in the repository
#' `README.md` and `docs/METHODS.md`.
#'
#' @keywords internal
"_PACKAGE"

## quiet R CMD check for ggplot2 non-standard evaluation
utils::globalVariables(c(
  "x", "y", "xend", "yend", "fill", "value", "label", "colour",
  "pathway", "effect", "neglog", "sig", "row_i", "col_i", "lab", "grp", "size"
))
