## Simple high-level API -- one call per stage, or one call for everything.
##
##     library(pathwaytheme)
##
##     # ---- full pipeline ----
##     pt_run_pipeline("examples/moh_sm.yaml")
##
##     # ---- or stage by stage ----
##     fm      <- pt_load_matrix("expr.tsv", metadata = "meta.tsv")
##     scores  <- pt_enrich(fm, backend = "ssgsea", gmt = "go_bp.gmt",
##                          geneset = "GO_BP")
##     groups  <- pt_group(scores, mode = "target", target_col = "grp",
##                         scope_col = "sample_id")
##     results <- pt_pca(scores, groups, size_col = "n_cells")
##     pt_figures(results, "out", geneset = "GO_BP")
##
## Every wrapper builds the underlying config for you; drop down to
## `pt_config()` plus the stage functions only when you need full control.

# ---- stage 1: input --------------------------------------------------------

#' Load a features x samples table
#'
#' @param matrix_path a pre-quantified features x samples table (`.tsv`,
#'   `.csv`, or any delimited text).
#' @param metadata optional samples x attributes table.
#' @param sep field separator for delimited text.
#' @return A [pt_feature_matrix()].
#' @export
pt_load_matrix <- function(matrix_path, metadata = NULL, sep = "\t") {
  pt_load_input(pt_input_config(kind = "matrix", matrix_path = matrix_path,
                                metadata_path = metadata, sep = sep))
}

#' Aggregate single-cell `.h5ad` files into per-cluster pseudobulk
#'
#' @param directory a directory of `.h5ad` files.
#' @param files an explicit vector of `.h5ad` paths.
#' @param cluster_col,sample_col obs columns naming the cluster and the sample.
#' @param aggregate `"mean"` or `"sum"`.
#' @return A [pt_feature_matrix()].
#' @export
pt_load_h5ad <- function(directory = NULL, files = NULL,
                         cluster_col = "seurat_clusters",
                         sample_col = "Sample_name", aggregate = "mean") {
  pt_load_input(pt_input_config(kind = "h5ad", h5ad_dir = directory,
                                h5ad_paths = .pt_or(files, character()),
                                cluster_col = cluster_col,
                                sample_col = sample_col, aggregate = aggregate))
}


# ---- optional: select samples by a metadata value --------------------------

#' Keep only the samples whose metadata matches
#'
#' Works on a [pt_feature_matrix()] (before enrichment) or a
#' [pt_score_matrix()] (after), e.g. to run the rest of the pipeline on one
#' tumour type.
#'
#' @param obj a feature or score matrix.
#' @param column a metadata column.
#' @param keep values of `column` to retain.
#' @return An object of the same class as `obj`.
#' @export
pt_filter_samples <- function(obj, column, keep) {
  keep <- unique(as.character(keep))
  meta <- obj$metadata$table
  if (!column %in% colnames(meta))
    stop("metadata column '", column, "' not found; available: ",
         paste(colnames(meta), collapse = ", "), call. = FALSE)
  values <- as.character(meta[colnames(obj$data), column])
  cols <- colnames(obj$data)[values %in% keep]
  if (!length(cols))
    stop("no samples have ", column, " in ",
         paste(sort(keep, method = "radix"), collapse = ", "), call. = FALSE)
  sub <- obj$data[, cols, drop = FALSE]
  new_meta <- pt_sample_metadata(meta[cols, , drop = FALSE])
  if (inherits(obj, "pt_feature_matrix")) return(pt_feature_matrix(sub, new_meta))
  pt_score_matrix(sub, new_meta, backend = obj$backend, geneset = obj$geneset,
                  coverage = obj$coverage)
}


# ---- stage 2: enrichment ---------------------------------------------------

#' Run an enrichment backend
#'
#' @param fm a [pt_feature_matrix()].
#' @param backend `"ssgsea"`, `"enrichr"` or `"goslim"`.
#' @param geneset label used for output naming.
#' @param gmt path to a `.gmt` gene-set file.
#' @param cache_dir reuse a previously computed score matrix if present.
#' @param verbose print progress and coverage warnings.
#' @param ... further fields of [pt_enrichment_config()], e.g. `top_n = 200`.
#' @return A [pt_score_matrix()].
#' @export
pt_enrich <- function(fm, backend = "ssgsea", geneset = "geneset", gmt = NULL,
                      cache_dir = NULL, verbose = TRUE, ...) {
  cfg <- pt_enrichment_config(backend = backend, geneset = geneset,
                              gmt_path = gmt, cache_dir = cache_dir, ...)
  pt_run_enrichment(fm, cfg, verbose = verbose)
}


# ---- stage 3: grouping -----------------------------------------------------

#' Resolve a comparison grouping
#'
#' @param sm a [pt_score_matrix()].
#' @param mode `"target"`, `"existing"` or `"auto"`.
#' @param target_col metadata column to compare / colour by.
#' @param scope_col metadata column partitioning samples into independent PCA
#'   runs.
#' @param ... further fields of [pt_grouping_config()], e.g. `n_clusters = 4`.
#' @return A [pt_grouping()].
#' @export
pt_group <- function(sm, mode = "target", target_col = NULL, scope_col = NULL,
                     ...) {
  cfg <- pt_grouping_config(mode = mode, target_col = target_col,
                            scope_col = scope_col, ...)
  pt_resolve_grouping(sm, cfg)
}


# ---- stage 4: PCA ----------------------------------------------------------

#' Run one PCA per scope
#'
#' If `grouping` is omitted, all samples form a single unlabelled scope.
#'
#' @param sm a [pt_score_matrix()].
#' @param grouping an optional [pt_grouping()].
#' @param display_col,secondary_col,size_col optional metadata columns providing
#'   display labels, a secondary colour track and dot sizes.
#' @param ... further fields of [pt_pca_config()], e.g. `k_max = 8`.
#' @return A list of [pt_pca_result()] objects.
#' @export
pt_pca <- function(sm, grouping = NULL, display_col = NULL,
                   secondary_col = NULL, size_col = NULL, ...) {
  if (is.null(grouping)) {
    labels <- stats::setNames(rep("all", ncol(sm$data)), colnames(sm$data))
    grouping <- pt_grouping(labels = labels, label_name = "group")
  }
  pt_run_pca_scopes(sm, grouping, pt_pca_config(...), display_col = display_col,
                    sublabel_col = secondary_col, size_col = size_col)
}


# ---- optional stage: differential pathway analysis -------------------------

#' Compare pathway scores between groups
#'
#' With no `reference` and no explicit `contrasts` this runs one-vs-rest for
#' every group; pass `reference = "X"` for every-group-vs-X or
#' `contrasts = list(c("A", "B"), ...)`.
#'
#' @param sm a [pt_score_matrix()].
#' @param grouping an optional [pt_grouping()] supplying the labels.
#' @param method `"welch"`, `"mannwhitney"` or `"moderated_t"`.
#' @param group_col a metadata column to group by instead of `grouping`.
#' @param reference reference group.
#' @param contrasts explicit `list(c(case, reference), ...)`.
#' @param ... further fields of [pt_diff_config()].
#' @return A [pt_diff_result()].
#' @export
pt_diff <- function(sm, grouping = NULL, method = "welch", group_col = NULL,
                    reference = NULL, contrasts = NULL, ...) {
  cfg <- pt_diff_config(enabled = TRUE, method = method, group_col = group_col,
                        reference = reference, contrasts = contrasts, ...)
  pt_run_diff(sm, cfg, grouping)
}

#' Roll a result table up into broad categories
#'
#' @param table a [pt_diff_result()], its `$table`, or any long data.frame with
#'   a pathway column.
#' @param mapping a term -> category named vector, or a path to a TSV/CSV
#'   (loaded via [pt_load_category_map()]).
#' @param key_col the pathway column.
#' @param value_col the numeric column to aggregate.
#' @param group_cols extra columns to aggregate within; defaults to
#'   `"comparison"` when present.
#' @param stat `"mean"`, `"median"` or `"sum"`.
#' @param significance_col pass `"fdr"` or `"p_value"` to split each category
#'   into significant vs non-significant pathways.
#' @param alpha significance threshold.
#' @return A data.frame, one row per (group..., \[significance,\] category).
#' @export
pt_summarize_categories <- function(table, mapping, key_col = "pathway",
                                    value_col = "effect", group_cols = NULL,
                                    stat = "mean", significance_col = NULL,
                                    alpha = 0.05) {
  if (inherits(table, "pt_diff_result")) table <- table$table
  if (is.character(mapping) && length(mapping) == 1L && file.exists(mapping))
    mapping <- pt_load_category_map(mapping, key_col = key_col)
  if (is.null(group_cols) && "comparison" %in% colnames(table))
    group_cols <- "comparison"
  pt_summarize_by_category(table, mapping, key_col = key_col,
                           value_col = value_col, group_cols = group_cols,
                           stat = stat, significance_col = significance_col,
                           alpha = alpha)
}

#' Write a z-scored pathway x sample QC heatmap
#'
#' @param sm a [pt_score_matrix()].
#' @param out_path output PDF path.
#' @param label_col optional metadata column for the colour strip.
#' @param max_pathways cap on the most-variable pathways.
#' @return The written path, or `NULL`.
#' @export
pt_sanity_heatmap <- function(sm, out_path, label_col = NULL,
                              max_pathways = 200L) {
  pt_render_sanity_heatmap(sm, out_path, label_col = label_col,
                           max_pathways = max_pathways)
}


# ---- stage 5: figures + tables ---------------------------------------------

#' Write the 11 PDFs and 3 TSVs per PCA result
#'
#' Output goes to `<out_dir>/<geneset>/sample_pca/<scope>/`.
#'
#' @param results a [pt_pca_result()] or a list of them.
#' @param out_dir output root.
#' @param geneset label used in the path and file names.
#' @param tables,figures toggle the two output kinds.
#' @return A named list of written paths, keyed by scope.
#' @export
pt_figures <- function(results, out_dir, geneset = "analysis", tables = TRUE,
                       figures = TRUE) {
  if (inherits(results, "pt_pca_result")) results <- list(results)
  viz <- pt_viz_config(make_figures = figures, make_tables = tables)
  pcfg <- pt_pca_config()
  root <- file.path(out_dir, geneset, "sample_pca")
  written <- list()
  for (res in results) {
    d <- file.path(root, res$scope)
    prefix <- paste0(geneset, "_", res$scope)
    title <- paste0(geneset, " / ", res$scope)
    paths <- character()
    if (tables) paths <- c(paths, pt_write_pca_tables(res, d, prefix, pcfg))
    if (figures) paths <- c(paths, pt_render_pca_figures(res, d, prefix, title,
                                                         viz))
    written[[res$scope]] <- paths
  }
  written
}


# ---- all-in-one ------------------------------------------------------------

#' Run the full pipeline
#'
#' @param config a YAML path, a [pt_config()], or `NULL` to build entirely from
#'   `...`.
#' @param verbose print stage-by-stage progress.
#' @param ... dotted overrides, e.g. `"input.matrix_path" = "x.tsv"`.
#' @return A [pt_pipeline_result()].
#' @examples
#' \dontrun{
#' pt_run_pipeline("my_config.yaml")
#' pt_run_pipeline(
#'   "input.matrix_path" = "expr.tsv",
#'   "input.metadata_path" = "meta.tsv",
#'   "enrichment.backend" = "ssgsea",
#'   "enrichment.gmt_path" = "go_bp.gmt",
#'   "grouping.target_col" = "treatment",
#'   output_dir = "results")
#' }
#' @export
pt_run_pipeline <- function(config = NULL, verbose = TRUE, ...) {
  cfg <- if (is.character(config) && length(config) == 1L) {
    pt_config_from_yaml(config)
  } else if (inherits(config, "pt_config")) {
    config
  } else {
    pt_config()
  }
  cfg <- pt_config_set(cfg, ...)
  pt_run(cfg, verbose = verbose)
}
