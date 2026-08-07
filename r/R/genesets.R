## Gene-set (.gmt) parsing + gene-universe helpers.

#' Parse a GMT gene-set file
#'
#' GMT format is one gene set per line:
#' `term <tab> description <tab> gene1 <tab> gene2 ...`.
#'
#' @param path path to a `.gmt` file.
#' @return A named list mapping each term to its character vector of genes.
#' @export
pt_parse_gmt <- function(path) {
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  gs <- list()
  for (line in lines) {
    sp <- strsplit(sub("\r$", "", line), "\t", fixed = TRUE)[[1]]
    if (length(sp) >= 3L) {
      genes <- sp[-c(1L, 2L)]
      gs[[sp[[1]]]] <- genes[nzchar(genes)]
    }
  }
  if (!length(gs)) stop("No gene sets parsed from ", path, call. = FALSE)
  gs
}

#' Union of all genes across every gene set
#'
#' @param gs a named list of gene sets, e.g. from [pt_parse_gmt()].
#' @return A character vector of unique gene identifiers.
#' @export
pt_gene_universe <- function(gs) unique(unlist(gs, use.names = FALSE))

#' Decide whether upper-casing improves gene overlap
#'
#' Mirrors the Python implementation: if upper-casing both sides matches more
#' genes than the raw comparison, both the expression index and the gene sets
#' are upper-cased; otherwise both are left alone.
#'
#' @param feature_names character vector of matrix row names.
#' @param gs a named list of gene sets.
#' @return A list with `transform` (a function to apply to feature names) and
#'   `gene_sets` (the possibly upper-cased sets).
#' @export
pt_harmonize_case <- function(feature_names, gs) {
  universe <- pt_gene_universe(gs)
  raw_ov <- length(intersect(as.character(feature_names), universe))
  up_ov <- length(intersect(unique(toupper(feature_names)), unique(toupper(universe))))
  if (up_ov > raw_ov) {
    list(transform = function(s) toupper(as.character(s)),
         gene_sets = lapply(gs, toupper))
  } else {
    list(transform = function(s) as.character(s), gene_sets = gs)
  }
}
