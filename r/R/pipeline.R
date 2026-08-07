## End-to-end orchestrator: io -> enrichment -> grouping -> pca -> viz.
##
## `pt_run(config)` executes all stages and writes, per PCA scope, the 11
## figures and 3 tables into `<output_dir>/<geneset>/sample_pca/<scope>/`.

#' Result of a full pipeline run
#'
#' @param feature_matrix,score_matrix,grouping,pca_results the stage outputs.
#' @param output_dir where results were written.
#' @param written named list of written paths, keyed by scope.
#' @param diff_result,category_summary optional differential-stage outputs.
#' @return An object of class `pt_pipeline_result`.
#' @export
pt_pipeline_result <- function(feature_matrix, score_matrix, grouping,
                               pca_results, output_dir, written = list(),
                               diff_result = NULL, category_summary = NULL) {
  structure(list(feature_matrix = feature_matrix, score_matrix = score_matrix,
                 grouping = grouping, pca_results = pca_results,
                 output_dir = output_dir, written = written,
                 diff_result = diff_result,
                 category_summary = category_summary),
            class = "pt_pipeline_result")
}

#' @export
print.pt_pipeline_result <- function(x, ...) {
  cat("<pt_pipeline_result>\n")
  cat("  score matrix : ", nrow(x$score_matrix$data), " terms x ",
      ncol(x$score_matrix$data), " samples\n", sep = "")
  cat("  PCA results  : ", length(x$pca_results), "\n", sep = "")
  cat("  files written: ", length(unlist(x$written)), "\n", sep = "")
  cat("  output_dir   : ", x$output_dir, "\n", sep = "")
  invisible(x)
}

## Observation display labels: explicit column, else `cl{cluster_id}`, else NULL
.pt_build_display_labels <- function(sm, grouping_cfg) {
  meta <- sm$metadata$table
  if (.pt_is_set(grouping_cfg$display_col) &&
      grouping_cfg$display_col %in% colnames(meta))
    return(stats::setNames(as.character(meta[[grouping_cfg$display_col]]),
                           rownames(meta)))
  if ("cluster_id" %in% colnames(meta))
    return(stats::setNames(paste0("cl", as.character(meta[["cluster_id"]])),
                           rownames(meta)))
  NULL
}

#' Run the whole pipeline from a configuration
#'
#' @param config a [pt_config()].
#' @param verbose print stage-by-stage progress.
#' @return A [pt_pipeline_result()].
#' @export
pt_run <- function(config, verbose = TRUE) {
  log <- function(...) if (verbose) message(...)

  # 1 -- input
  log(sprintf("[1/5] loading input (kind=%s) ...", config$input$kind))
  fm <- pt_load_input(config$input)
  log(sprintf("      feature matrix: %d features x %d samples",
              nrow(fm$data), ncol(fm$data)))
  if (.pt_is_set(config$input$filter_col) &&
      .pt_is_set(config$input$filter_values)) {
    fm <- pt_filter_samples(fm, config$input$filter_col,
                            config$input$filter_values)
    log(sprintf("      filtered to %s in %s: %d samples",
                config$input$filter_col,
                paste(config$input$filter_values, collapse = ", "),
                ncol(fm$data)))
  }

  # 2 -- enrichment
  log(sprintf("[2/5] enrichment (backend=%s, geneset=%s) ...",
              config$enrichment$backend, config$enrichment$geneset))
  sm <- pt_run_enrichment(fm, config$enrichment, verbose = verbose)
  log(sprintf("      score matrix: %d terms x %d samples",
              nrow(sm$data), ncol(sm$data)))
  if (!is.null(sm$coverage)) {
    n_req <- nrow(sm$coverage)
    n_drop <- sum(sm$coverage$status != PT_SCORED)
    log(sprintf("      gene sets: %d of %d scored%s", n_req - n_drop, n_req,
                if (n_drop) sprintf(", %d dropped for insufficient coverage",
                                    n_drop) else ""))
  }

  # 3 -- grouping
  log(sprintf("[3/5] grouping (mode=%s) ...", config$grouping$mode))
  grouping <- pt_resolve_grouping(sm, config$grouping)
  scopes <- pt_scopes(grouping)
  log(sprintf("      %d labels, %d scope(s)",
              length(unique(grouping$labels)), length(scopes)))

  # 4 -- PCA per scope
  log("[4/5] PCA per scope ...")
  display <- .pt_build_display_labels(sm, config$grouping)
  # display labels are keyed by sample id; pt_run_pca_scopes expects a metadata
  # column name, so stash the computed vector on the metadata table
  disp_col <- NULL
  if (!is.null(display)) {
    sm$metadata$table[["__display__"]] <- unname(display[pt_samples(sm)])
    disp_col <- "__display__"
  }
  results <- pt_run_pca_scopes(sm, grouping, config$pca,
                               display_col = disp_col,
                               sublabel_col = config$grouping$secondary_col,
                               size_col = config$grouping$size_col)
  log(sprintf("      %d PCA result(s) (scopes with >= %d observations)",
              length(results), config$pca$min_observations))

  # 5 -- figures + tables
  out_root <- file.path(config$output_dir, config$enrichment$geneset, "sample_pca")
  written <- list()
  if (isTRUE(config$viz$make_figures) || isTRUE(config$viz$make_tables)) {
    log("[5/5] writing figures + tables ...")
    for (res in results) {
      out_dir <- file.path(out_root, res$scope)
      prefix <- paste0(config$enrichment$geneset, "_", res$scope)
      title <- paste0(config$enrichment$geneset, " / ", res$scope)
      paths <- character()
      if (isTRUE(config$viz$make_tables))
        paths <- c(paths, pt_write_pca_tables(res, out_dir, prefix, config$pca))
      if (isTRUE(config$viz$make_figures))
        paths <- c(paths, pt_render_pca_figures(res, out_dir, prefix, title,
                                                config$viz))
      written[[res$scope]] <- paths
      log(sprintf("      %s: %d file(s)", res$scope, length(paths)))
    }
  } else {
    log("[5/5] viz disabled")
  }

  analysis_root <- file.path(config$output_dir, config$enrichment$geneset)
  base_prefix <- config$enrichment$geneset

  # gene-set coverage: which requested sets were scored, and why not
  if (!is.null(sm$coverage) && isTRUE(config$viz$make_tables)) {
    cov_path <- file.path(analysis_root,
                          paste0(base_prefix, "_geneset_coverage.tsv"))
    written[["_coverage"]] <- .pt_write_tsv(sm$coverage, cov_path)
  }

  # optional -- full-matrix sanity heatmap (QC)
  if (isTRUE(config$viz$sanity_heatmap) && isTRUE(config$viz$make_figures)) {
    log("[qc] writing sanity heatmap ...")
    p <- pt_render_sanity_heatmap(
      sm, file.path(analysis_root, paste0(base_prefix, "_sanity_heatmap.pdf")),
      label_col = config$grouping$target_col,
      max_pathways = config$viz$sanity_max_pathways)
    if (!is.null(p)) written[["_qc"]] <- p
  }

  # optional -- differential pathway analysis (parallel to PCA)
  diff_result <- NULL
  category_summary <- NULL
  if (isTRUE(config$diff$enabled)) {
    log(sprintf("[diff] differential analysis (method=%s) ...",
                config$diff$method))
    diff_result <- pt_run_diff(sm, config$diff, grouping)
    n_sig <- sum(diff_result$table$fdr < 0.05, na.rm = TRUE)
    log(sprintf("       %d comparison(s), %d pathway(s) FDR<0.05",
                length(pt_comparisons(diff_result)), n_sig))
    dpaths <- character()
    if (isTRUE(config$viz$make_tables))
      dpaths <- c(dpaths, pt_write_diff_table(diff_result, analysis_root,
                                              base_prefix))
    if (isTRUE(config$viz$make_figures))
      dpaths <- c(dpaths, pt_render_diff_figures(diff_result, analysis_root,
                                                 base_prefix,
                                                 top_n = config$diff$top_n))

    # optional -- roll differential results up into broad categories
    if (isTRUE(config$categories$enabled) &&
        .pt_is_set(config$categories$map_path)) {
      log("[categories] summarizing by category ...")
      mapping <- pt_load_category_map(config$categories$map_path,
                                      key_col = config$categories$key_col,
                                      category_col = config$categories$category_col)
      sig_col <- if (isTRUE(config$categories$split_significance))
        config$categories$significance_col else NULL
      category_summary <- pt_summarize_by_category(
        diff_result$table, mapping, key_col = "pathway", value_col = "effect",
        group_cols = "comparison", stat = config$categories$stat,
        unmapped_label = config$categories$unmapped_label,
        significance_col = sig_col, alpha = config$categories$alpha)
      if (isTRUE(config$viz$make_tables))
        dpaths <- c(dpaths, pt_write_category_table(category_summary,
                                                    analysis_root, base_prefix))
    }
    written[["_diff"]] <- dpaths
  }

  pt_pipeline_result(fm, sm, grouping, results, config$output_dir, written,
                     diff_result = diff_result,
                     category_summary = category_summary)
}
