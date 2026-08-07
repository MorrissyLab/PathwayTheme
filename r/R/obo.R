## Minimal OBO ontology reader + the GO-slim mapping used by goatools.
##
## Enough of the Gene Ontology graph to reproduce `goatools.mapslim`: term ids,
## names, and the `is_a` / `part_of` parent edges that `paths_to_top` walks.
## Obsolete terms are skipped, alt_ids are aliased to their primary id.

#' Read an OBO ontology
#'
#' @param path a `.obo` file (e.g. `go-basic.obo`, `goslim_generic.obo`).
#' @return An object of class `pt_obo` with `parents` (named list of parent ids),
#'   `names` (id -> term name) and `alt` (alt_id -> primary id).
#' @export
pt_read_obo <- function(path) {
  lines <- readLines(path, warn = FALSE)
  starts <- which(lines == "[Term]")
  if (!length(starts)) stop("no [Term] stanzas found in ", path, call. = FALSE)
  stanza_ends <- c(starts[-1L] - 1L, length(lines))

  ids <- character(length(starts))
  nms <- character(length(starts))
  parents <- vector("list", length(starts))
  alt <- list()
  obsolete <- logical(length(starts))

  for (i in seq_along(starts)) {
    block <- lines[(starts[[i]] + 1L):stanza_ends[[i]]]
    block <- block[nzchar(block)]
    id <- .pt_obo_field(block, "id: ")
    if (!length(id)) next
    ids[[i]] <- id[[1L]]
    nm <- .pt_obo_field(block, "name: ")
    nms[[i]] <- if (length(nm)) nm[[1L]] else id[[1L]]
    obsolete[[i]] <- any(startsWith(block, "is_obsolete: true"))
    is_a <- sub(" !.*$", "", .pt_obo_field(block, "is_a: "))
    rel <- .pt_obo_field(block, "relationship: part_of ")
    rel <- sub(" !.*$", "", rel)
    parents[[i]] <- unique(trimws(c(is_a, rel)))
    for (a in .pt_obo_field(block, "alt_id: ")) alt[[a]] <- id[[1L]]
  }

  keep <- nzchar(ids) & !obsolete
  ids <- ids[keep]; nms <- nms[keep]; parents <- parents[keep]
  names(parents) <- ids
  names(nms) <- ids
  structure(list(parents = parents, names = nms,
                 alt = unlist(alt, use.names = TRUE)), class = "pt_obo")
}

#' @export
print.pt_obo <- function(x, ...) {
  cat("<pt_obo> ", length(x$names), " terms\n", sep = "")
  invisible(x)
}

.pt_obo_field <- function(block, prefix) {
  hit <- block[startsWith(block, prefix)]
  if (!length(hit)) return(character())
  substring(hit, nchar(prefix) + 1L)
}

.pt_obo_resolve <- function(obo, id) {
  if (id %in% names(obo$parents)) return(id)
  if (!is.null(obo$alt) && id %in% names(obo$alt)) return(unname(obo$alt[[id]]))
  NA_character_
}

## Ancestors of `id` including itself, memoised in `cache` (an environment).
.pt_obo_ancestors <- function(obo, id, cache) {
  hit <- cache[[id]]
  if (!is.null(hit)) return(hit)
  # iterative walk so deep ontologies cannot blow the R stack
  seen <- character()
  stack <- id
  while (length(stack)) {
    cur <- stack[[1L]]; stack <- stack[-1L]
    if (cur %in% seen) next
    seen <- c(seen, cur)
    up <- obo$parents[[cur]]
    if (length(up)) stack <- c(stack, setdiff(up, seen))
  }
  cache[[id]] <- seen
  seen
}

#' Map a GO term onto its GO-slim ancestors
#'
#' Reproduces `goatools.mapslim`: `all` is every slim term on a path from the
#' term to the root (the term itself included), and `direct` is the subset that
#' no other member of `all` sits below -- the most specific slim ancestors.
#'
#' @param go_id a GO id.
#' @param obo the full ontology, from [pt_read_obo()].
#' @param slim the GO-slim ontology, from [pt_read_obo()].
#' @param cache an environment used to memoise ancestor sets across calls.
#' @return A list with `direct` and `all` character vectors of slim ids.
#' @export
pt_mapslim <- function(go_id, obo, slim, cache = new.env(parent = emptyenv())) {
  id <- .pt_obo_resolve(obo, go_id)
  if (is.na(id)) return(list(direct = character(), all = character()))
  anc <- .pt_obo_ancestors(obo, id, cache)
  all_slims <- intersect(anc, names(slim$names))
  if (length(all_slims) <= 1L) return(list(direct = all_slims, all = all_slims))
  # a slim term is "covered" when another slim ancestor sits below it
  covered <- character()
  for (b in all_slims) {
    anc_b <- setdiff(.pt_obo_ancestors(obo, b, cache), b)
    covered <- union(covered, intersect(anc_b, all_slims))
  }
  list(direct = setdiff(all_slims, covered), all = all_slims)
}

#' Map GO IDs to a representative GO-slim term
#'
#' The downstream counterpart of the GO-slim enrichment backend: turns pathway
#' names that are GO ids into broad categories for [pt_summarize_categories()].
#' Terms with no slim ancestor are omitted, so callers fall back to the
#' unmapped label.
#'
#' @param go_ids character vector of GO ids.
#' @param obo_path path to the full ontology (`go-basic.obo`).
#' @param slim_obo_path path to the slim ontology (`goslim_generic.obo`).
#' @return A named character vector mapping GO id -> slim term name.
#' @export
pt_goslim_category_map <- function(go_ids, obo_path, slim_obo_path) {
  obo <- pt_read_obo(obo_path)
  slim <- pt_read_obo(slim_obo_path)
  cache <- new.env(parent = emptyenv())
  slim_cache <- new.env(parent = emptyenv())   # keyed by id too, so keep separate
  out <- character()
  for (gid in unique(as.character(go_ids))) {
    m <- pt_mapslim(gid, obo, slim, cache)
    if (!length(m$all)) next
    # the broadest representative: shallowest slim term, ties broken by id
    depth <- vapply(m$all, function(s)
      length(.pt_obo_ancestors(slim, s, slim_cache)), 1L)
    best <- m$all[order(depth, m$all, method = "radix")][[1L]]
    out[[gid]] <- unname(slim$names[[best]])
  }
  out
}
