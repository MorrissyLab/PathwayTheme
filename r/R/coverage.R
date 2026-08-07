## Gene-set coverage accounting and missing-value policy.
##
## Two things the pipeline must not do silently, both of which first showed up
## when the same gene sets were run over a mass-spectrometry proteome rather
## than a transcriptome:
##
## **Sets can vanish.**  A `.gmt` set whose members are largely absent from the
## feature matrix falls below `min_gene_set_size` and is dropped by the backend.
## `pt_geneset_coverage()` records the fate of every requested set.
##
## **Values can be missing.**  A TMT proteome is distributed as signed
## log-ratios with ~20% of entries absent.  ssGSEA is rank-based, so an unstated
## fill puts those entries at an arbitrary point in the ranking.
## `pt_apply_missing_policy()` makes the choice explicit and leaves the default
## a no-op so existing results are unchanged.

#' Fates a requested gene set can meet
#' @export
PT_SCORED <- "scored"
#' @rdname PT_SCORED
#' @export
PT_NO_OVERLAP <- "dropped_no_overlap"
#' @rdname PT_SCORED
#' @export
PT_TOO_SMALL <- "dropped_below_min_size"
#' @rdname PT_SCORED
#' @export
PT_TOO_LARGE <- "dropped_above_max_size"

PT_MISSING_POLICIES <- c("as_is", "zero", "low")

#' Per requested gene set: how many members are measurable, and its fate
#'
#' `features` is the feature index of the matrix actually handed to the backend.
#' Matching is done on the upper-cased symbol, mirroring [pt_harmonize_case()],
#' so the count here is the count the backend will see.
#'
#' @param features character vector of feature (gene) names.
#' @param gene_sets named list of gene sets.
#' @param min_size,max_size size bounds applied by the backend.
#' @return A data.frame with `pathway`, `n_genes`, `n_present`,
#'   `fraction_present` and `status`, sorted by status then pathway.
#' @export
pt_geneset_coverage <- function(features, gene_sets, min_size, max_size) {
  present <- unique(toupper(as.character(features)))
  n_genes <- integer(length(gene_sets))
  n_present <- integer(length(gene_sets))
  for (i in seq_along(gene_sets)) {
    uniq <- unique(toupper(as.character(gene_sets[[i]])))
    n_genes[[i]] <- length(uniq)
    n_present[[i]] <- length(intersect(uniq, present))
  }
  status <- ifelse(n_present == 0L, PT_NO_OVERLAP,
            ifelse(n_present < min_size, PT_TOO_SMALL,
            ifelse(n_present > max_size, PT_TOO_LARGE, PT_SCORED)))
  out <- data.frame(
    pathway = names(gene_sets),
    n_genes = n_genes,
    n_present = n_present,
    fraction_present = ifelse(n_genes > 0L, n_present / n_genes, 0),
    status = status,
    stringsAsFactors = FALSE)
  out <- out[order(out$status, out$pathway, method = "radix"), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Warn about gene sets that will not be scored
#'
#' @param coverage a table from [pt_geneset_coverage()].
#' @param geneset label used in the warning text.
#' @param verbose set `FALSE` to stay silent.
#' @return The dropped rows, invisibly usable by the caller.
#' @export
pt_warn_on_dropped <- function(coverage, geneset = "", verbose = TRUE) {
  dropped <- coverage[coverage$status != PT_SCORED, , drop = FALSE]
  if (!nrow(dropped) || !verbose) return(dropped)
  label <- if (nzchar(geneset)) paste0(" (", geneset, ")") else ""
  head_rows <- utils::head(dropped, 10L)
  lines <- sprintf("        %s: %d/%d members measured [%s]",
                   head_rows$pathway, head_rows$n_present, head_rows$n_genes,
                   head_rows$status)
  more <- if (nrow(dropped) <= 10L) "" else
    sprintf("\n        ... and %d more", nrow(dropped) - 10L)
  warning(sprintf("%d of %d gene sets%s will not be scored:\n%s%s",
                  nrow(dropped), nrow(coverage), label,
                  paste(lines, collapse = "\n"), more),
          call. = FALSE)
  dropped
}

#' Drop under-observed features, then resolve the remaining missing values
#'
#' `min_observed_fraction` drops any feature measured in fewer than that
#' fraction of samples (0, the default, drops nothing beyond all-missing rows).
#' `policy` then decides the rest:
#'
#' * `as_is` leave them for the backend, warning that something else will decide.
#'   The default -- it changes no existing result.
#' * `zero` fill with 0 explicitly, making the de-facto behaviour a stated choice.
#' * `low` fill with one unit below each sample's own observed minimum, so an
#'   unmeasured feature ranks below every measured one.
#'
#' @param data a numeric matrix, features x samples.
#' @param policy one of `"as_is"`, `"zero"`, `"low"`.
#' @param min_observed_fraction minimum fraction of samples a feature must be
#'   measured in.
#' @param verbose print/warn about what was done.
#' @return The resolved matrix.
#' @export
pt_apply_missing_policy <- function(data, policy = "as_is",
                                    min_observed_fraction = 0,
                                    verbose = TRUE) {
  if (!policy %in% PT_MISSING_POLICIES)
    stop("enrichment$missing must be one of ",
         paste(PT_MISSING_POLICIES, collapse = ", "), ", got '", policy, "'",
         call. = FALSE)
  out <- data
  n_feat_before <- nrow(out)

  frac_obs <- rowMeans(!is.na(out))
  keep <- if (min_observed_fraction <= 0) frac_obs > 0 else
    frac_obs >= min_observed_fraction
  if (!all(keep)) {
    out <- out[keep, , drop = FALSE]
    if (verbose)
      message(sprintf("      missing data: dropped %d of %d features observed in <%s of samples",
                      n_feat_before - nrow(out), n_feat_before,
                      paste0(format(round(max(min_observed_fraction, 0) * 100)), "%")))
  }

  n_missing <- sum(is.na(out))
  if (n_missing == 0L) return(out)
  frac <- n_missing / length(out)
  pct <- sprintf("%.1f%%", frac * 100)

  if (policy == "as_is") {
    if (verbose)
      warning(sprintf(paste0("%d missing value(s) (%s of the matrix) reach the enrichment ",
                             "backend with enrichment$missing='as_is'; unmeasured features ",
                             "are ranked as if they were average. Set enrichment$missing='low' ",
                             "to rank unmeasured features below measured ones, or 'zero' to ",
                             "state the current behaviour explicitly."),
                      n_missing, pct), call. = FALSE)
    return(out)
  }
  if (policy == "zero") {
    out[is.na(out)] <- 0
    if (verbose)
      message(sprintf("      missing data: filled %d value(s) (%s) with 0",
                      n_missing, pct))
    return(out)
  }
  # policy == "low": one unit below each sample's own observed minimum
  floor_vals <- suppressWarnings(apply(out, 2L, min, na.rm = TRUE)) - 1
  floor_vals[!is.finite(floor_vals)] <- 0
  for (j in seq_len(ncol(out))) {
    na_j <- is.na(out[, j])
    if (any(na_j)) out[na_j, j] <- floor_vals[[j]]
  }
  if (verbose)
    message(sprintf(paste0("      missing data: filled %d value(s) (%s) one unit ",
                           "below each sample's observed minimum"), n_missing, pct))
  out
}
