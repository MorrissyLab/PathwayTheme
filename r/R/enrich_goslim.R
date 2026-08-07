## GoSlim backend: coarse GO-slim category scores per sample.
##
## Two modes, both producing a `slim-term x sample` matrix:
##
## * **gmt mode** (default, offline, testable): a GO-slim `.gmt` where each term
##   is a slim category and members are its genes.  Per sample, the score for a
##   slim term is the fraction (or count) of the sample's selected genes that
##   are members of that term.
## * **obo mode**: with `obo_path` + `slim_obo_path` + `gene2go_path`, genes are
##   mapped to their GO annotations and rolled up to slim ancestors via
##   [pt_mapslim()]; the score is again the fraction/count of the sample's genes
##   hitting each slim term.

.pt_score_from_membership <- function(genes, term_members, how) {
  gset <- unique(as.character(genes))
  n <- max(1L, length(gset))
  hits <- vapply(term_members, function(members) length(intersect(gset, members)),
                 1L)
  hits <- hits[hits > 0L]
  if (!length(hits)) return(stats::setNames(numeric(), character()))
  # as.double() would drop the term names the caller pivots on
  if (identical(how, "fraction")) hits / n
  else stats::setNames(as.double(hits), names(hits))
}

#' Read an id-to-GOs association file
#'
#' The `IdToGos` format used by goatools: two tab-separated columns,
#' `id <tab> GO:0000001;GO:0000002`.
#'
#' @param path the association file.
#' @return A named list mapping each id to its character vector of GO ids.
#' @export
pt_read_id2gos <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines <- lines[nzchar(trimws(lines)) & !startsWith(lines, "#")]
  sp <- strsplit(lines, "\t", fixed = TRUE)
  ids <- vapply(sp, function(x) x[[1L]], "")
  gos <- lapply(sp, function(x) {
    if (length(x) < 2L) return(character())
    g <- trimws(strsplit(x[[2L]], ";", fixed = TRUE)[[1L]])
    g[nzchar(g)]
  })
  stats::setNames(gos, ids)
}

## Map genes -> GO-slim terms; returns list(slim_term_name = character(genes)).
.pt_build_obo_membership <- function(config, universe) {
  obo <- pt_read_obo(config$obo_path)
  slim <- pt_read_obo(config$slim_obo_path)
  assoc <- pt_read_id2gos(config$gene2go_path)
  assoc <- assoc[names(assoc) %in% universe]

  cache <- new.env(parent = emptyenv())
  term_members <- list()
  for (gene in names(assoc)) {
    slims <- character()
    for (go in assoc[[gene]]) {
      slims <- union(slims, pt_mapslim(go, obo, slim, cache)$direct)
    }
    for (s in slims) {
      nm <- if (!is.null(slim$names[[s]])) unname(slim$names[[s]]) else s
      term_members[[nm]] <- c(term_members[[nm]], gene)
    }
  }
  lapply(term_members, unique)
}

.pt_backend_goslim <- function(fm, config, verbose = TRUE) {
  universe <- unique(pt_features(fm))
  if (.pt_is_set(config$obo_path) && .pt_is_set(config$slim_obo_path) &&
      .pt_is_set(config$gene2go_path)) {
    term_members <- .pt_build_obo_membership(config, universe)
  } else if (.pt_is_set(config$gmt_path)) {
    term_members <- lapply(pt_parse_gmt(config$gmt_path), unique)
  } else {
    stop("GoSlim needs either enrichment$gmt_path (a GO-slim gmt) or ",
         "obo_path + slim_obo_path + gene2go_path", call. = FALSE)
  }

  per_sample <- list()
  for (sample in pt_samples(fm)) {
    col <- fm$data[, sample]
    names(col) <- rownames(fm$data)
    genes <- pt_select_genes(col, config)
    if (!length(genes)) next
    ser <- .pt_score_from_membership(genes, term_members, config$goslim_score)
    if (length(ser)) per_sample[[sample]] <- ser
  }
  if (!length(per_sample))
    stop("GoSlim produced no results for any sample", call. = FALSE)
  .pt_bind_score_columns(per_sample, pt_samples(fm))
}
