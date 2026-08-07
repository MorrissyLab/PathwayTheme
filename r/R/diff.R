## Differential pathway analysis over a score matrix + grouping.
##
## Runs a two-group test (see diff_stats.R) for one or more comparisons and
## returns a single long `pt_diff_result` table.
##
## Comparison selection
## --------------------
## * explicit `contrasts = list(c(case, ref), ...)` -- run exactly those, OR
## * a `reference` group                            -- every other group vs it, OR
## * neither                                        -- one-vs-rest for every group.

.pt_diff_labels <- function(sm, config, grouping) {
  if (.pt_is_set(config$group_col)) {
    if (!pt_meta_has(sm$metadata, config$group_col))
      stop("diff$group_col '", config$group_col, "' not in metadata", call. = FALSE)
    labels <- as.character(pt_meta_get(sm$metadata, config$group_col))
    names(labels) <- rownames(sm$metadata$table)
    return(list(labels = labels[pt_samples(sm)], name = config$group_col))
  }
  if (!is.null(grouping))
    return(list(labels = grouping$labels[pt_samples(sm)],
                name = grouping$label_name))
  stop("differential analysis needs diff$group_col or a grouping", call. = FALSE)
}

## [(case, reference_or_NA), ...].  A reference of NA means one-vs-rest.
.pt_plan_comparisons <- function(config, groups) {
  if (.pt_is_set(config$contrasts)) {
    return(lapply(config$contrasts,
                  function(cr) c(as.character(cr[[1L]]), as.character(cr[[2L]]))))
  }
  if (.pt_is_set(config$reference)) {
    ref <- as.character(config$reference)
    return(lapply(setdiff(groups, ref), function(g) c(g, ref)))
  }
  lapply(groups, function(g) c(g, NA_character_))
}

#' Compare pathway scores between groups
#'
#' @param sm a [pt_score_matrix()].
#' @param config a diff section, from [pt_diff_config()].
#' @param grouping an optional [pt_grouping()] providing the labels when
#'   `config$group_col` is not set.
#' @return A [pt_diff_result()] with all comparisons stacked in one long table.
#' @export
pt_run_diff <- function(sm, config, grouping = NULL) {
  resolved <- .pt_diff_labels(sm, config, grouping)
  labels <- resolved$labels
  groups <- sort(unique(labels[!is.na(labels)]), method = "radix")
  plan <- .pt_plan_comparisons(config, groups)

  data <- sm$data                        # pathways x samples
  pathways <- rownames(data)
  frames <- list()

  for (pair in plan) {
    case <- pair[[1L]]
    ref <- pair[[2L]]
    case_cols <- names(labels)[!is.na(labels) & labels == case]
    if (is.na(ref)) {
      ref_cols <- names(labels)[is.na(labels) | labels != case]
      ref_name <- "rest"
    } else {
      ref_cols <- names(labels)[!is.na(labels) & labels == ref]
      ref_name <- as.character(ref)
    }
    case_cols <- case_cols[case_cols %in% colnames(data)]
    ref_cols <- ref_cols[ref_cols %in% colnames(data)]
    if (length(case_cols) < config$min_group_size ||
        length(ref_cols) < config$min_group_size) next

    case_mat <- data[, case_cols, drop = FALSE]
    ref_mat <- data[, ref_cols, drop = FALSE]
    res <- pt_run_test(config$method, case_mat, ref_mat)
    fdr <- pt_benjamini_hochberg(res$p_value)

    frame <- data.frame(
      comparison = paste0(case, "_vs_", ref_name),
      pathway = pathways,
      method = config$method,
      n_case = length(case_cols),
      n_reference = length(ref_cols),
      mean_case = rowMeans(case_mat, na.rm = TRUE),
      mean_reference = rowMeans(ref_mat, na.rm = TRUE),
      effect = res$effect,
      statistic = res$statistic,
      p_value = res$p_value,
      fdr = fdr,
      direction = ifelse(res$effect >= 0, "up", "down"),
      stringsAsFactors = FALSE)
    frame <- frame[order(frame$fdr, frame$p_value, na.last = TRUE,
                         method = "radix"), , drop = FALSE]
    frames[[length(frames) + 1L]] <- frame
  }

  table <- if (length(frames)) do.call(rbind, frames) else .pt_empty_diff_table()
  rownames(table) <- NULL
  pt_diff_result(table = table, method = config$method,
                 group_col = .pt_or(resolved$name, "group"))
}

.pt_empty_diff_table <- function() {
  data.frame(comparison = character(), pathway = character(),
             method = character(), n_case = integer(), n_reference = integer(),
             mean_case = double(), mean_reference = double(), effect = double(),
             statistic = double(), p_value = double(), fdr = double(),
             direction = character(), stringsAsFactors = FALSE)
}
