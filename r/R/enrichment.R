## Enrichment-backend interface, registry, caching, and dispatcher.
##
## Every backend consumes a `pt_feature_matrix` and returns a `pt_score_matrix`
## with the identical `terms x samples` shape, so nothing downstream branches on
## the backend.  Results are cached to a TSV keyed on backend + geneset.

## Registry entries carry `run` and `per_sample`.  `per_sample` is TRUE when a
## sample's scores depend only on that sample -- ssGSEA ranks genes within a
## sample, the ORA backend takes that sample's top-n genes, GO-slim scores that
## sample's own profile.  Only then may the cache be filled in sample by sample;
## a cross-sample-normalised backend must set it FALSE so a partial cache hit
## recomputes everything.
.pt_backends <- new.env(parent = emptyenv())

#' Register an enrichment backend
#'
#' @param name backend name used in `enrichment$backend`.
#' @param run a function `(fm, config, verbose)` returning a terms x samples
#'   numeric matrix whose columns align to `pt_samples(fm)`.
#' @param per_sample `TRUE` when a sample's scores do not depend on which other
#'   samples were run; only such backends may fill a partial cache.
#' @return The name, invisibly.
#' @export
pt_register_backend <- function(name, run, per_sample = TRUE) {
  assign(tolower(name), list(run = run, per_sample = per_sample),
         envir = .pt_backends)
  invisible(name)
}

#' Look up a registered enrichment backend
#'
#' @param name backend name.
#' @return A list with `run` and `per_sample`.
#' @export
pt_get_backend <- function(name) {
  key <- tolower(name)
  if (!exists(key, envir = .pt_backends, inherits = FALSE))
    stop("Unknown enrichment backend '", name, "'; available: ",
         paste(sort(ls(.pt_backends)), collapse = ", "), call. = FALSE)
  get(key, envir = .pt_backends, inherits = FALSE)
}

.pt_register_builtin_backends <- function() {
  pt_register_backend("ssgsea", .pt_backend_ssgsea, per_sample = TRUE)
  pt_register_backend("enrichr", .pt_backend_enrichr, per_sample = TRUE)
  pt_register_backend("goslim", .pt_backend_goslim, per_sample = TRUE)
}

.onLoad <- function(libname, pkgname) {
  .pt_register_builtin_backends()
}


# ---- caching ---------------------------------------------------------------

.pt_cache_path <- function(config) {
  if (!.pt_is_set(config$cache_dir)) return(NULL)
  # the missing-value policy changes the scores, so it is part of the key --
  # but only when non-default, so existing caches stay valid
  suffix <- ""
  if (!identical(.pt_or(config$missing, "as_is"), "as_is"))
    suffix <- paste0(suffix, "_miss-", config$missing)
  if (.pt_or(config$min_observed_fraction, 0) > 0)
    suffix <- paste0(suffix, "_obs", format(config$min_observed_fraction))
  file.path(config$cache_dir,
            paste0(config$geneset, "_", config$backend, suffix, "_scores.tsv"))
}

.pt_read_cache <- function(path) {
  raw <- utils::read.table(path, sep = "\t", header = TRUE, row.names = 1L,
                           check.names = FALSE, quote = "\"",
                           comment.char = "", stringsAsFactors = FALSE)
  .pt_to_numeric_matrix(raw)
}

.pt_write_cache <- function(data, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(data, path, sep = "\t", quote = FALSE, col.names = NA,
                     na = "")
  invisible(path)
}


# ---- dispatcher ------------------------------------------------------------

.pt_prepare_features <- function(fm, config, verbose = TRUE) {
  data <- pt_apply_missing_policy(fm$data, .pt_or(config$missing, "as_is"),
                                  .pt_or(config$min_observed_fraction, 0),
                                  verbose = verbose)
  if (identical(dim(data), dim(fm$data)) && identical(data, fm$data)) return(fm)
  pt_feature_matrix(data, pt_align_to(fm$metadata, colnames(data)))
}

.pt_coverage_for <- function(fm, config, verbose = TRUE) {
  if (!.pt_is_set(config$gmt_path)) return(NULL)
  gs <- tryCatch(pt_parse_gmt(config$gmt_path), error = function(e) NULL)
  if (is.null(gs)) return(NULL)
  cov <- pt_geneset_coverage(pt_features(fm), gs, config$min_gene_set_size,
                             config$max_gene_set_size)
  pt_warn_on_dropped(cov, geneset = config$geneset, verbose = verbose)
  cov
}

#' Run (or load from cache) the configured enrichment backend
#'
#' The cache is keyed on backend + geneset only, so two configs that differ in
#' which samples they select share one file.  A requested sample that is not in
#' the cache is therefore scored and added, never silently dropped: with a
#' per-sample backend a sample's scores do not depend on which other samples
#' were run, so filling the cache in is exact.
#'
#' @param fm a [pt_feature_matrix()].
#' @param config an enrichment section, from [pt_enrichment_config()].
#' @param verbose print progress and coverage warnings.
#' @return A [pt_score_matrix()].
#' @export
pt_run_enrichment <- function(fm, config, verbose = TRUE) {
  cache <- .pt_cache_path(config)
  backend <- pt_get_backend(config$backend)
  wanted <- colnames(fm$data)
  meta <- fm$metadata            # fm is narrowed below on a partial cache hit

  # resolve missing values and account for gene-set coverage up front, so a
  # cache hit reports the same coverage a fresh run would
  fm <- .pt_prepare_features(fm, config, verbose = verbose)
  coverage <- .pt_coverage_for(fm, config, verbose = verbose)

  cached <- NULL
  if (!is.null(cache) && file.exists(cache)) {
    cached <- .pt_read_cache(cache)
    have <- intersect(wanted, colnames(cached))
    if (length(have) == length(wanted)) {
      data <- cached[, wanted, drop = FALSE]
      data <- data[rowSums(!is.na(data)) > 0L, , drop = FALSE]
      return(pt_score_matrix(data, pt_align_to(meta, wanted),
                             backend = config$backend, geneset = config$geneset,
                             coverage = coverage))
    }
    missing <- setdiff(wanted, colnames(cached))
    if (!isTRUE(backend$per_sample)) {
      cached <- NULL                       # must recompute as a whole
    } else {
      if (verbose)
        message(sprintf("      cache: %d of %d samples hit, scoring %d more",
                        length(have), length(wanted), length(missing)))
      fm <- pt_feature_matrix(fm$data[, missing, drop = FALSE],
                              pt_align_to(fm$metadata, missing))
    }
  }

  data <- backend$run(fm, config, verbose)
  data <- .pt_as_matrix(data)
  # keep only real samples, drop all-NA rows/cols
  data <- data[rowSums(!is.na(data)) > 0L, , drop = FALSE]
  data <- data[, colSums(!is.na(data)) > 0L, drop = FALSE]
  data <- data[, colnames(data) %in% rownames(fm$metadata$table), drop = FALSE]

  if (!is.null(cached)) data <- .pt_outer_join_terms(cached, data)

  if (!is.null(cache)) .pt_write_cache(data, cache)

  cols <- intersect(wanted, colnames(data))
  missing <- setdiff(wanted, colnames(data))
  if (length(missing))
    stop(length(missing), " requested sample(s) produced no scores, e.g. ",
         paste(utils::head(missing, 3), collapse = ", "), call. = FALSE)
  data <- data[, cols, drop = FALSE]
  data <- data[rowSums(!is.na(data)) > 0L, , drop = FALSE]
  pt_score_matrix(data, pt_align_to(meta, cols),
                  backend = config$backend, geneset = config$geneset,
                  coverage = coverage)
}

## Outer join two score matrices on their term names: a newly scored sample may
## carry a term the cached samples did not, and vice versa.
.pt_outer_join_terms <- function(left, right) {
  terms <- union(rownames(left), rownames(right))
  cols <- c(colnames(left), colnames(right))
  out <- matrix(NA_real_, nrow = length(terms), ncol = length(cols),
                dimnames = list(terms, cols))
  out[rownames(left), colnames(left)] <- left
  out[rownames(right), colnames(right)] <- right
  out
}
