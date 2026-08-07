## Configuration objects for the pipeline.
##
## Everything the pipeline needs is expressed as a nested list so a run is fully
## described by one `pt_config()` (loadable from / dumpable to YAML).  Section
## names, field names and defaults are identical to the Python package, so the
## same YAML file drives either implementation.

# ---- PCA constants (shared with the Python implementation) -----------------

#' PCA display/rollup constants
#'
#' Upper bound on principal components kept, PCs combined into one observation's
#' signature, pathways per observation in the signature panels, and pathways
#' tabulated per component in the TSV.
#' @export
PT_K_MAX <- 10L
#' @rdname PT_K_MAX
#' @export
PT_TOP_PCS_PER_CLUSTER <- 3L
#' @rdname PT_K_MAX
#' @export
PT_TOP_PATHWAYS_PER_CLUST <- 25L
#' @rdname PT_K_MAX
#' @export
PT_TOP_PATHWAYS_PER_PC <- 50L

#' Blue-white-red colour ramp anchors
#'
#' Matches `colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))`.
#' @export
PT_BWR_COLORS <- c("#0000FF", "#FFFFFF", "#FF0000")


# ---- section constructors --------------------------------------------------

#' Pipeline configuration sections
#'
#' Each function builds one section of a [pt_config()].  Defaults reproduce the
#' Python package's defaults exactly, so a YAML file written by either
#' implementation drives the other.
#'
#' @param kind `"matrix"` or `"h5ad"`.
#' @param matrix_path features x samples table (tsv/csv/parquet-free; see
#'   [pt_load_matrix()]).
#' @param metadata_path samples x attributes table.
#' @param sep field separator for plain-text tables.
#' @param h5ad_dir directory of `.h5ad` files.
#' @param h5ad_paths explicit character vector of `.h5ad` files.
#' @param cluster_col,sample_col obs columns used to build pseudobulk columns.
#' @param aggregate `"mean"` or `"sum"` per-cluster aggregation.
#' @param filter_col,filter_values keep only samples whose
#'   `metadata[[filter_col]]` is in `filter_values`.
#'
#' @return A named list carrying the section's fields.
#' @export
pt_input_config <- function(kind = "matrix",
                            matrix_path = NULL,
                            metadata_path = NULL,
                            sep = "\t",
                            h5ad_dir = NULL,
                            h5ad_paths = character(),
                            cluster_col = "seurat_clusters",
                            sample_col = "Sample_name",
                            aggregate = "mean",
                            filter_col = NULL,
                            filter_values = NULL) {
  list(kind = kind, matrix_path = matrix_path, metadata_path = metadata_path,
       sep = sep, h5ad_dir = h5ad_dir, h5ad_paths = h5ad_paths,
       cluster_col = cluster_col, sample_col = sample_col,
       aggregate = aggregate, filter_col = filter_col,
       filter_values = filter_values)
}

#' @rdname pt_input_config
#' @param backend `"ssgsea"`, `"enrichr"` or `"goslim"`.
#' @param geneset label used for output naming.
#' @param gmt_path a `.gmt` gene-set file.
#' @param min_gene_set_size,max_gene_set_size size bounds on scored gene sets.
#' @param threads worker threads (ssGSEA).
#' @param cache_dir reuse a previously computed score matrix if present.
#' @param missing missing-value policy: `"as_is"`, `"zero"` or `"low"`.
#' @param min_observed_fraction drop features measured in fewer samples than this.
#' @param weight ssGSEA exponent on the rank metric.
#' @param enrichr_library named online Enrichr library (online mode).
#' @param gene_selection `"top_n"` or `"threshold"`.
#' @param top_n,threshold how many / which genes form each sample's gene list.
#' @param enrichr_score `"-log10_padj"` or `"combined_score"`.
#' @param organism organism for the online Enrichr service.
#' @param obo_path,slim_obo_path,gene2go_path,id2sym_path GO-slim ontology inputs.
#' @param goslim_score `"fraction"` or `"count"`.
#' @export
pt_enrichment_config <- function(backend = "ssgsea",
                                 geneset = "GO_BP",
                                 gmt_path = NULL,
                                 min_gene_set_size = 5L,
                                 max_gene_set_size = 5000L,
                                 threads = 8L,
                                 cache_dir = NULL,
                                 missing = "as_is",
                                 min_observed_fraction = 0,
                                 weight = 0.25,
                                 enrichr_library = NULL,
                                 gene_selection = "top_n",
                                 top_n = 200L,
                                 threshold = NULL,
                                 enrichr_score = "-log10_padj",
                                 organism = "human",
                                 obo_path = NULL,
                                 slim_obo_path = NULL,
                                 gene2go_path = NULL,
                                 id2sym_path = NULL,
                                 goslim_score = "fraction") {
  list(backend = backend, geneset = geneset, gmt_path = gmt_path,
       min_gene_set_size = as.integer(min_gene_set_size),
       max_gene_set_size = as.integer(max_gene_set_size),
       threads = as.integer(threads), cache_dir = cache_dir,
       missing = missing, min_observed_fraction = min_observed_fraction,
       weight = weight, enrichr_library = enrichr_library,
       gene_selection = gene_selection, top_n = as.integer(top_n),
       threshold = threshold, enrichr_score = enrichr_score,
       organism = organism, obo_path = obo_path,
       slim_obo_path = slim_obo_path, gene2go_path = gene2go_path,
       id2sym_path = id2sym_path, goslim_score = goslim_score)
}

#' @rdname pt_input_config
#' @param mode `"target"`, `"existing"` or `"auto"`.
#' @param target_col metadata column to compare / colour by.
#' @param scope_col metadata column partitioning samples into independent PCA runs.
#' @param secondary_col optional secondary colouring track.
#' @param size_col metadata column driving dot sizes.
#' @param display_col metadata column for observation display labels.
#' @param method `"kmeans"` or `"hierarchical"` (auto mode).
#' @param n_clusters fixed k; `NULL` searches `k_min:k_max` by silhouette.
#' @param k_min,k_max silhouette search range.
#' @export
pt_grouping_config <- function(mode = "target",
                               target_col = NULL,
                               scope_col = NULL,
                               secondary_col = NULL,
                               size_col = NULL,
                               display_col = NULL,
                               method = "kmeans",
                               n_clusters = NULL,
                               k_min = 2L,
                               k_max = 10L) {
  list(mode = mode, target_col = target_col, scope_col = scope_col,
       secondary_col = secondary_col, size_col = size_col,
       display_col = display_col, method = method, n_clusters = n_clusters,
       k_min = as.integer(k_min), k_max = as.integer(k_max))
}

#' @rdname pt_input_config
#' @param top_pcs_per_cluster PCs combined into one observation's signature.
#' @param top_pathways_per_cluster pathways kept per observation signature.
#' @param top_pathways_per_pc pathways tabulated per component.
#' @param min_observations minimum columns for a meaningful PCA.
#' @param min_features minimum variable features for a meaningful PCA.
#' @param random_state kept for parity with the Python config (unused in R).
#' @export
pt_pca_config <- function(k_max = PT_K_MAX,
                          top_pcs_per_cluster = PT_TOP_PCS_PER_CLUSTER,
                          top_pathways_per_cluster = PT_TOP_PATHWAYS_PER_CLUST,
                          top_pathways_per_pc = PT_TOP_PATHWAYS_PER_PC,
                          min_observations = 3L,
                          min_features = 10L,
                          random_state = 42L) {
  list(k_max = as.integer(k_max),
       top_pcs_per_cluster = as.integer(top_pcs_per_cluster),
       top_pathways_per_cluster = as.integer(top_pathways_per_cluster),
       top_pathways_per_pc = as.integer(top_pathways_per_pc),
       min_observations = as.integer(min_observations),
       min_features = as.integer(min_features),
       random_state = as.integer(random_state))
}

#' @rdname pt_input_config
#' @param enabled run this optional stage.
#' @param group_col metadata column holding the comparison groups.
#' @param reference reference group; `NULL` means one-vs-rest.
#' @param contrasts explicit `list(c(case, reference), ...)`.
#' @param min_group_size minimum samples per side of a comparison.
#' @export
pt_diff_config <- function(enabled = FALSE,
                           group_col = NULL,
                           method = "welch",
                           reference = NULL,
                           contrasts = NULL,
                           min_group_size = 2L,
                           top_n = 25L) {
  list(enabled = enabled, group_col = group_col, method = method,
       reference = reference, contrasts = contrasts,
       min_group_size = as.integer(min_group_size), top_n = as.integer(top_n))
}

#' @rdname pt_input_config
#' @param map_path TSV/CSV term -> category mapping.
#' @param key_col,category_col column names within `map_path`.
#' @param stat `"mean"`, `"median"` or `"sum"`.
#' @param unmapped_label category assigned to terms absent from the map.
#' @param split_significance split each category into significant vs not.
#' @param significance_col `"fdr"` or `"p_value"`.
#' @param alpha significance threshold.
#' @export
pt_category_config <- function(enabled = FALSE,
                               map_path = NULL,
                               key_col = "pathway",
                               category_col = "category",
                               stat = "mean",
                               unmapped_label = "Other",
                               split_significance = TRUE,
                               significance_col = "fdr",
                               alpha = 0.05) {
  list(enabled = enabled, map_path = map_path, key_col = key_col,
       category_col = category_col, stat = stat,
       unmapped_label = unmapped_label,
       split_significance = split_significance,
       significance_col = significance_col, alpha = alpha)
}

#' @rdname pt_input_config
#' @param make_figures,make_tables toggle the two output kinds.
#' @param dpi raster resolution for rasterised heatmap bodies.
#' @param formats output formats (only `"pdf"` is implemented).
#' @param sanity_heatmap draw the full pathway x sample QC heatmap.
#' @param sanity_max_pathways cap on most-variable pathways in that heatmap.
#' @export
pt_viz_config <- function(make_figures = TRUE,
                          make_tables = TRUE,
                          dpi = 120L,
                          formats = "pdf",
                          sanity_heatmap = FALSE,
                          sanity_max_pathways = 200L) {
  list(make_figures = make_figures, make_tables = make_tables,
       dpi = as.integer(dpi), formats = formats,
       sanity_heatmap = sanity_heatmap,
       sanity_max_pathways = as.integer(sanity_max_pathways))
}


# ---- the whole config ------------------------------------------------------

#' Build a full pipeline configuration
#'
#' @param input,enrichment,grouping,pca,diff,categories,viz section lists, each
#'   built by the corresponding `pt_*_config()` constructor.
#' @param output_dir root directory for figures and tables.
#'
#' @return An object of class `pt_config`.
#' @examples
#' cfg <- pt_config(output_dir = "results")
#' cfg <- pt_config_set(cfg, "enrichment.backend" = "ssgsea")
#' cfg$enrichment$backend
#' @export
pt_config <- function(input = pt_input_config(),
                      enrichment = pt_enrichment_config(),
                      grouping = pt_grouping_config(),
                      pca = pt_pca_config(),
                      diff = pt_diff_config(),
                      categories = pt_category_config(),
                      viz = pt_viz_config(),
                      output_dir = "pathwaytheme_output") {
  structure(list(input = input, enrichment = enrichment, grouping = grouping,
                 pca = pca, diff = diff, categories = categories, viz = viz,
                 output_dir = output_dir),
            class = "pt_config")
}

#' @export
print.pt_config <- function(x, ...) {
  cat("<pt_config>\n")
  cat("  input     : kind=", x$input$kind,
      if (!is.null(x$input$matrix_path)) paste0(" matrix=", x$input$matrix_path) else "",
      "\n", sep = "")
  cat("  enrichment: backend=", x$enrichment$backend,
      " geneset=", x$enrichment$geneset, "\n", sep = "")
  cat("  grouping  : mode=", x$grouping$mode,
      " target_col=", .pt_or(x$grouping$target_col, "<none>"), "\n", sep = "")
  cat("  diff      : ", if (isTRUE(x$diff$enabled)) x$diff$method else "disabled",
      "\n", sep = "")
  cat("  output_dir: ", x$output_dir, "\n", sep = "")
  invisible(x)
}

#' Override configuration fields with dotted keys
#'
#' @param config a [pt_config()].
#' @param ... named overrides using dotted section keys, e.g.
#'   `"input.matrix_path" = "expr.tsv"`.  A key without a dot sets a top-level
#'   field such as `output_dir`.
#' @return The modified `pt_config`.
#' @export
pt_config_set <- function(config, ...) {
  overrides <- list(...)
  if (!length(overrides)) return(config)
  nms <- names(overrides)
  if (is.null(nms) || any(!nzchar(nms)))
    stop("pt_config_set() overrides must all be named", call. = FALSE)
  for (i in seq_along(overrides)) {
    parts <- strsplit(nms[[i]], ".", fixed = TRUE)[[1]]
    if (length(parts) == 1L) {
      config[parts[[1]]] <- list(overrides[[i]])
    } else {
      section <- parts[[1]]
      field <- paste(parts[-1], collapse = ".")
      if (is.null(config[[section]]))
        stop("unknown config section '", section, "'", call. = FALSE)
      config[[section]][field] <- list(overrides[[i]])
    }
  }
  config
}

#' Read / write a pipeline configuration as YAML
#'
#' The YAML layout is identical to the Python package's, so a config written by
#' one implementation runs unchanged in the other.
#'
#' @param path YAML file path.
#' @param config a [pt_config()].
#' @return `pt_config_from_yaml()` returns a `pt_config`;
#'   `pt_config_to_yaml()` returns `path` invisibly.
#' @export
pt_config_from_yaml <- function(path) {
  d <- yaml::read_yaml(path)
  if (is.null(d)) d <- list()
  pt_config_from_list(d)
}

#' @rdname pt_config_from_yaml
#' @param d a nested list, e.g. from [yaml::read_yaml()].
#' @export
pt_config_from_list <- function(d) {
  sec <- function(name, ctor) {
    vals <- d[[name]]
    if (is.null(vals)) vals <- list()
    defaults <- ctor()
    unknown <- setdiff(names(vals), names(defaults))
    if (length(unknown))
      stop("unknown key(s) in config section '", name, "': ",
           paste(unknown, collapse = ", "), call. = FALSE)
    for (k in names(vals)) defaults[k] <- list(vals[[k]])
    defaults
  }
  pt_config(
    input = sec("input", pt_input_config),
    enrichment = sec("enrichment", pt_enrichment_config),
    grouping = sec("grouping", pt_grouping_config),
    pca = sec("pca", pt_pca_config),
    diff = sec("diff", pt_diff_config),
    categories = sec("categories", pt_category_config),
    viz = sec("viz", pt_viz_config),
    output_dir = if (is.null(d$output_dir)) "pathwaytheme_output" else d$output_dir
  )
}

#' @rdname pt_config_from_yaml
#' @export
pt_config_to_yaml <- function(config, path) {
  body <- unclass(config)
  # yaml drops NULLs unless they are explicitly represented
  writeLines(yaml::as.yaml(body, handlers = list(
    `NULL` = function(x) structure("~", class = "verbatim")
  )), con = path, useBytes = TRUE)
  invisible(path)
}


# ---- small internal helpers ------------------------------------------------

.pt_or <- function(x, default) if (is.null(x) || !length(x)) default else x

.pt_is_set <- function(x) !is.null(x) && length(x) > 0L && !all(is.na(x)) &&
  !identical(x, "")
