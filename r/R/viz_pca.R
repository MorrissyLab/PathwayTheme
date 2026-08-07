## PCA figure generation -- the 11 PDFs.
##
## Blocks:
##   1. observation x PC scores heatmap
##   2. PC1-PC2 biplot + pairwise grid  (x2 colourings: primary + secondary label)
##   3. per-observation pathway signatures
##   4. per-PC loading bars + pathway x PC heatmap + clustered heatmap
##   5. union of top +/- loadings per PC  (clustermap + signature view)
##
## All functions consume a `pt_pca_result`.  Colouring uses the result's
## `labels` (primary) and `sublabels` (secondary).

## `cl3 (tumor)` -> `cl3`; otherwise the label unchanged
.pt_short_label <- function(lbl) {
  s <- as.character(lbl)
  ifelse(grepl(" (", s, fixed = TRUE), sub(" \\(.*$", "", s), s)
}

## matplotlib scatter `s` (area in points^2) -> ggplot2 `size` (mm)
.pt_pt2mm <- function(s) sqrt(s) * 25.4 / 72

.pt_blank_theme <- function(base = 9) {
  ggplot2::theme_bw(base_size = base) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   legend.key = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(hjust = 0, size = base + 2))
}


# ---------------------------------------------------------------------------
# Block 1 -- observation x PC scores heatmap
# ---------------------------------------------------------------------------
.pt_draw_scores_heatmap <- function(res, title, out_path) {
  scores <- res$scores
  var_expl <- res$variance_explained
  n_obs <- nrow(scores)
  k <- ncol(scores)
  class_lut <- pt_class_colours(sort(unique(res$labels), method = "radix"))
  row_colors <- pt_label_colours(class_lut, res$labels)

  vmax <- max(1, stats::quantile(abs(scores), 0.98, na.rm = TRUE, names = FALSE))
  body_w <- 0.35 * k + 1.5
  body_h <- 0.32 * n_obs + 1.5

  pt_save_pdf(out_path, width = max(7, body_w + 4), height = max(5, body_h + 2),
              draw = function() {
    pt_draw_clustermap(
      scores, vmin = -vmax, vmax = vmax,
      row_cluster = n_obs > 1L, col_cluster = FALSE,
      row_colors = row_colors,
      xlabels = sprintf("%s\n(%.0f%%)", colnames(scores), var_expl * 100),
      ylabels = rownames(scores),
      title = title, cbar_label = "PC score",
      legend = list(colours = class_lut[sort(names(class_lut), method = "radix")],
                    title = .pt_or(res$labels_name, "group")),
      fontsize_x = 8, fontsize_y = 7)
  })
}


# ---------------------------------------------------------------------------
# Block 2 -- biplot + pairwise grid
# ---------------------------------------------------------------------------
.pt_draw_biplot <- function(scores, loadings, var_expl, point_labels,
                            group_labels, sizes, lut, legend_title, title_str,
                            out_path) {
  cols <- pt_label_colours(lut, group_labels)
  n <- log1p(sizes)
  pt_sizes <- 40 + 320 * (n / max(1, max(n)))

  x <- scores[, "PC1"]
  y <- scores[, "PC2"]
  df <- data.frame(x = x, y = y, lab = point_labels, grp = group_labels,
                   size = .pt_pt2mm(pt_sizes), stringsAsFactors = FALSE)

  lp <- t(loadings[c("PC1", "PC2"), , drop = FALSE])
  mag <- sqrt(lp[, "PC1"]^2 + lp[, "PC2"]^2)
  top10 <- utils::head(order(mag, decreasing = TRUE, method = "radix"), 10L)
  span_x <- .pt_nonzero((max(x) - min(x)) * 0.5)
  span_y <- .pt_nonzero((max(y) - min(y)) * 0.5)
  arrow_scale <- 0.7 * min(span_x, span_y) / max(abs(lp[top10, , drop = FALSE]))
  arr <- data.frame(xend = lp[top10, "PC1"] * arrow_scale,
                    yend = lp[top10, "PC2"] * arrow_scale,
                    lab = pt_truncate(rownames(lp)[top10], 55L),
                    stringsAsFactors = FALSE)

  present <- sort(unique(group_labels[nzchar(group_labels)]), method = "radix")
  breaks <- present
  values <- pt_label_colours(lut, present)
  if (any(!nzchar(group_labels))) {
    breaks <- c(breaks, "")
    values <- c(values, "#dddddd")
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_hline(yintercept = 0, colour = "#bbbbbb", linewidth = 0.25) +
    ggplot2::geom_vline(xintercept = 0, colour = "#bbbbbb", linewidth = 0.25) +
    ggplot2::geom_segment(data = arr,
                          ggplot2::aes(x = 0, y = 0, xend = xend, yend = yend),
                          inherit.aes = FALSE, colour = "#888888",
                          linewidth = 0.3, alpha = 0.7,
                          arrow = grid::arrow(length = grid::unit(2, "mm"),
                                              type = "closed")) +
    ggplot2::geom_point(ggplot2::aes(fill = grp, size = size), shape = 21,
                        colour = "#222222", stroke = 0.3, alpha = 0.9) +
    ggrepel::geom_text_repel(data = arr,
                             ggplot2::aes(x = xend * 1.08, y = yend * 1.08,
                                          label = lab),
                             inherit.aes = FALSE, size = 2.1, colour = "#444444",
                             fontface = "italic", segment.colour = "#cccccc",
                             segment.size = 0.15, max.overlaps = Inf) +
    ggrepel::geom_text_repel(ggplot2::aes(label = lab), size = 2.5,
                             colour = "#222222", fontface = "bold",
                             segment.colour = "#bbbbbb", segment.size = 0.15,
                             max.overlaps = Inf) +
    ggplot2::scale_size_identity() +
    ggplot2::scale_fill_manual(values = stats::setNames(values, breaks),
                               breaks = breaks,
                               labels = ifelse(nzchar(breaks), breaks, "(none)"),
                               name = legend_title) +
    ggplot2::labs(x = sprintf("PC1 (%.1f%% variance)", var_expl[[1L]] * 100),
                  y = sprintf("PC2 (%.1f%% variance)", var_expl[[2L]] * 100),
                  title = title_str) +
    ggplot2::guides(fill = ggplot2::guide_legend(
      override.aes = list(size = 3, shape = 21))) +
    .pt_blank_theme(10)

  pt_save_pdf(out_path, width = 12, height = 9, draw = function() print(p))
}

.pt_draw_pairwise <- function(scores, var_expl, pc_names, point_labels,
                              group_labels, sizes, lut, legend_title,
                              title_str, out_path) {
  k <- ncol(scores)
  n_show <- min(4L, k)
  show_pcs <- pc_names[seq_len(n_show)]
  n <- log1p(sizes)
  pt_sizes <- 30 + 200 * (n / max(1, max(n)))
  cols <- pt_label_colours(lut, group_labels)

  panels <- list()
  for (i in seq_len(n_show)) {
    for (j in seq_len(n_show)) {
      if (i == j) {
        panels[[length(panels) + 1L]] <- .pt_diag_panel(
          show_pcs[[i]], sprintf("%.1f%% var", var_expl[[i]] * 100))
        next
      }
      df <- data.frame(x = scores[, show_pcs[[j]]], y = scores[, show_pcs[[i]]],
                       lab = point_labels, fill = cols,
                       size = .pt_pt2mm(pt_sizes * 0.6), stringsAsFactors = FALSE)
      panels[[length(panels) + 1L]] <-
        ggplot2::ggplot(df, ggplot2::aes(x = x, y = y)) +
        ggplot2::geom_hline(yintercept = 0, colour = "#dddddd", linewidth = 0.2) +
        ggplot2::geom_vline(xintercept = 0, colour = "#dddddd", linewidth = 0.2) +
        ggplot2::geom_point(ggplot2::aes(fill = fill, size = size), shape = 21,
                            colour = "#333333", stroke = 0.2, alpha = 0.85) +
        ggplot2::geom_text(ggplot2::aes(label = lab), size = 1.6,
                           colour = "#333333", hjust = 0, vjust = 0,
                           nudge_x = 0.02) +
        ggplot2::scale_fill_identity() +
        ggplot2::scale_size_identity() +
        ggplot2::labs(x = if (i == n_show) show_pcs[[j]] else NULL,
                      y = if (j == 1L) show_pcs[[i]] else NULL) +
        .pt_blank_theme(7) +
        ggplot2::theme(plot.margin = grid::unit(c(1, 1, 1, 1), "mm"))
    }
  }

  present <- sort(unique(group_labels[nzchar(group_labels)]), method = "radix")
  legend_grob <- .pt_legend_grob(
    stats::setNames(pt_label_colours(lut, present), present),
    legend_title)

  pt_save_pdf(out_path, width = n_show * 3.6 + 2, height = n_show * 3.4,
              draw = function() {
    gridExtra::grid.arrange(
      gridExtra::arrangeGrob(grobs = panels, ncol = n_show),
      legend_grob, ncol = 2L,
      widths = grid::unit.c(grid::unit(1, "null"), grid::unit(1.7, "inches")),
      top = grid::textGrob(title_str, gp = grid::gpar(fontsize = 11),
                           x = 0.02, hjust = 0))
  })
}

.pt_diag_panel <- function(name, sub) {
  df <- data.frame(x = c(0.5, 0.5), y = c(0.55, 0.3),
                   label = c(name, sub), stringsAsFactors = FALSE)
  ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, label = label)) +
    ggplot2::geom_text(size = c(6, 3), fontface = c("bold", "plain"),
                       colour = c("#000000", "#555555")) +
    ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) +
    ggplot2::theme_void()
}

## A standalone legend drawn as a grob, for figures assembled with gridExtra.
.pt_legend_grob <- function(colours, title, fontsize = 8) {
  keys <- names(colours)
  n <- length(keys)
  if (!n) return(grid::nullGrob())
  ys <- grid::unit(0.95, "npc") - grid::unit(seq_len(n) * 4.2, "mm")
  grid::gTree(children = grid::gList(
    grid::textGrob(title, x = grid::unit(2, "mm"), y = grid::unit(0.95, "npc"),
                   just = c("left", "bottom"),
                   gp = grid::gpar(fontsize = fontsize + 1, fontface = "bold")),
    grid::rectGrob(x = grid::unit(2, "mm"), y = ys, width = grid::unit(3, "mm"),
                   height = grid::unit(3, "mm"), just = c("left", "centre"),
                   gp = grid::gpar(fill = unname(colours), col = NA)),
    grid::textGrob(ifelse(nzchar(keys), keys, "(none)"),
                   x = grid::unit(6, "mm"), y = ys, just = c("left", "centre"),
                   gp = grid::gpar(fontsize = fontsize))))
}


# ---------------------------------------------------------------------------
# Block 3 -- per-observation signature panels
# ---------------------------------------------------------------------------
## Rebuild {obs: (top_pcs, signed vector)} from the long signature table.
.pt_per_obs_from_signatures <- function(sig_df) {
  clusters <- unique(as.character(sig_df$cluster))
  out <- list()
  for (cl in clusters) {
    grp <- sig_df[sig_df$cluster == cl, , drop = FALSE]
    out[[cl]] <- list(
      top_pcs = strsplit(as.character(grp$top_pcs)[[1L]], ",", fixed = TRUE)[[1L]],
      values = stats::setNames(as.double(grp$signature_score),
                               as.character(grp$pathway)))
  }
  out
}

## A horizontal bar panel: named values sorted ascending, red up / blue down.
## A pathway can appear twice (the same term among both the largest and the
## smallest loadings of a short list), so bars are positioned numerically and
## the names are attached as axis labels rather than as factor levels.
.pt_bar_panel <- function(values, title, xlab, label_width = 70L,
                          label_size = 1.9, base = 8) {
  ord <- order(values, method = "radix")
  v <- values[ord]
  df <- data.frame(
    y = seq_along(v),
    value = as.double(v),
    fill = ifelse(v > 0, .PT_UP, .PT_DOWN), stringsAsFactors = FALSE)
  ggplot2::ggplot(df, ggplot2::aes(x = value, y = y, fill = fill)) +
    ggplot2::geom_col(width = 0.8, orientation = "y") +
    ggplot2::geom_vline(xintercept = 0, colour = "#333333", linewidth = 0.2) +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_y_continuous(breaks = df$y,
                                labels = pt_truncate(names(v), label_width),
                                expand = ggplot2::expansion(add = 0.6)) +
    ggplot2::labs(title = title, x = xlab, y = NULL) +
    .pt_blank_theme(base) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = label_size * 3),
                   axis.text.x = ggplot2::element_text(size = base - 1),
                   plot.title = ggplot2::element_text(size = base, hjust = 0))
}

.pt_draw_signature_panels <- function(res, title, out_path) {
  per_obs <- .pt_per_obs_from_signatures(res$signatures)
  n_obs <- length(per_obs)
  n_cols <- 4L
  n_rows <- ceiling(n_obs / n_cols)
  panels <- lapply(names(per_obs), function(cl) {
    entry <- per_obs[[cl]]
    .pt_bar_panel(entry$values,
                  sprintf("%s\nsignature PCs: %s", cl,
                          paste(entry$top_pcs, collapse = ", ")),
                  "signature score (red=up, blue=down)")
  })
  pt_save_pdf(out_path, width = n_cols * 5, height = n_rows * 4.2,
              draw = function() {
    gridExtra::grid.arrange(grobs = panels, ncol = n_cols,
                            top = grid::textGrob(title,
                                                 gp = grid::gpar(fontsize = 11)))
  })
}


# ---------------------------------------------------------------------------
# Block 4 -- per-PC loading bars + pathway x PC heatmaps
# ---------------------------------------------------------------------------
.pt_draw_pc_loadings <- function(res, tag_title, out_dir, prefix, top_n = 25L) {
  loadings <- res$loadings
  var_expl <- res$variance_explained
  pc_names <- pt_pc_names(res)
  k <- nrow(loadings)
  written <- character()

  n_cols <- min(3L, k)
  n_rows <- ceiling(k / n_cols)
  panels <- lapply(seq_along(pc_names), function(i) {
    lv <- loadings[i, ]
    ord <- order(lv, decreasing = TRUE, method = "radix")
    pos <- lv[utils::head(ord, top_n)]
    neg <- lv[utils::head(rev(ord), top_n)]
    .pt_bar_panel(c(neg, pos),
                  sprintf("%s (%.1f%% variance)", pc_names[[i]],
                          var_expl[[i]] * 100),
                  "loading (red = positive axis end, blue = negative)",
                  label_width = 65L, label_size = 1.6)
  })
  p <- file.path(out_dir, paste0(prefix, "_pca_pc_loadings.pdf"))
  written <- c(written, pt_save_pdf(p, width = n_cols * 6, height = n_rows * 4.5,
                                    draw = function() {
    gridExtra::grid.arrange(grobs = panels, ncol = n_cols,
      top = grid::textGrob(
        paste0(tag_title,
               " \\u2014 top pathway loadings per principal component"),
        gp = grid::gpar(fontsize = 11)))
  }))

  # pathway x PC heatmap (top 80 by max |loading|)
  max_abs <- apply(abs(loadings), 2L, max)
  top_paths <- utils::head(order(max_abs, decreasing = TRUE, method = "radix"), 80L)
  sub <- t(loadings[, top_paths, drop = FALSE])
  vmax <- stats::quantile(abs(sub), 0.98, na.rm = TRUE, names = FALSE)
  x_labels <- sprintf("%s\n(%.0f%%)", pc_names, var_expl * 100)
  fig_w <- max(7, 0.6 * k + 2)
  fig_h <- max(8, 0.18 * nrow(sub) + 1.5)

  p <- file.path(out_dir, paste0(prefix, "_pca_pathway_pc_heatmap.pdf"))
  written <- c(written, pt_save_pdf(p, width = fig_w + 2, height = fig_h,
                                    draw = function() {
    pt_draw_clustermap(sub, vmin = -vmax, vmax = vmax,
                       row_cluster = FALSE, col_cluster = FALSE,
                       xlabels = x_labels,
                       ylabels = pt_truncate(rownames(sub), 70L),
                       title = sprintf(paste0("%s\nTop %d pathways x top %d PCs ",
                                              "- loadings heatmap (ordered by max |loading|)"),
                                       tag_title, nrow(sub), k),
                       cbar_label = "loading", fontsize_x = 7, fontsize_y = 6)
  }))

  # clustered pathway x PC heatmap
  p <- file.path(out_dir, paste0(prefix, "_pca_pathway_pc_clustermap.pdf"))
  written <- c(written, pt_save_pdf(p, width = max(8, 0.7 * k + 3) + 2,
                                    height = max(9, 0.18 * nrow(sub) + 2),
                                    draw = function() {
    pt_draw_clustermap(sub, vmin = -vmax, vmax = vmax,
                       row_cluster = TRUE, col_cluster = FALSE,
                       xlabels = sprintf("%s (%.0f%%)", pc_names, var_expl * 100),
                       ylabels = pt_truncate(rownames(sub), 70L),
                       title = sprintf(paste0("%s\nTop %d pathways x top %d PCs ",
                                              "- loadings clustered by pathway"),
                                       tag_title, nrow(sub), k),
                       cbar_label = "loading", fontsize_x = 7, fontsize_y = 6)
  }))
  written[!vapply(written, is.null, TRUE)]
}


# ---------------------------------------------------------------------------
# Block 5 -- union of top +/- loadings per PC
# ---------------------------------------------------------------------------
.pt_draw_pc_loading_union <- function(res, tag_title, out_dir, prefix,
                                      top_per_pc_per_sign = 10L) {
  loadings <- res$loadings
  var_expl <- res$variance_explained
  pc_names <- pt_pc_names(res)
  k <- nrow(loadings)
  written <- character()

  union_paths <- character()
  for (i in seq_along(pc_names)) {
    lv <- loadings[i, ]
    ord <- order(lv, decreasing = TRUE, method = "radix")
    union_paths <- c(union_paths,
                     names(lv)[utils::head(ord, top_per_pc_per_sign)],
                     names(lv)[utils::head(rev(ord), top_per_pc_per_sign)])
  }
  union_paths <- unique(union_paths)
  sub_u <- t(loadings[, union_paths, drop = FALSE])
  vmax_u <- stats::quantile(abs(sub_u), 0.98, na.rm = TRUE, names = FALSE)

  x_labels <- sprintf("%s\n%.0f%% var", pc_names, var_expl * 100)
  fig_w <- max(10, 0.8 * k + 5)
  fig_h <- max(11, 0.18 * length(union_paths) + 3)

  # Plot A -- clustered union
  p <- file.path(out_dir, paste0(prefix, "_pca_pc_loading_union_clustermap.pdf"))
  written <- c(written, pt_save_pdf(p, width = fig_w, height = fig_h,
                                    draw = function() {
    pt_draw_clustermap(sub_u, vmin = -vmax_u, vmax = vmax_u,
                       row_cluster = TRUE, col_cluster = FALSE,
                       xlabels = x_labels,
                       ylabels = pt_truncate(rownames(sub_u), 55L),
                       title = sprintf(paste0("%s - union of top +/-%d loadings per ",
                                              "PC (%d pathways) - clustered by pathway"),
                                       tag_title, top_per_pc_per_sign,
                                       length(union_paths)),
                       cbar_label = "loading", fontsize_x = 8, fontsize_y = 6,
                       xlabel_rot = 35)
  }))

  # Plot B -- signature view (rows grouped by home PC + sign)
  home_idx <- apply(abs(sub_u), 1L, which.max)
  signed_home <- sub_u[cbind(seq_len(nrow(sub_u)), home_idx)]
  ord <- order(home_idx, as.integer(signed_home < 0), -abs(signed_home),
               method = "radix")
  sub_ord <- sub_u[ord, , drop = FALSE]
  pc_color <- stats::setNames(
    .PT_TAB10[((seq_along(pc_names) - 1L) %% 10L) + 1L], pc_names)
  row_colors <- unname(pc_color[pc_names[home_idx[ord]]])

  p <- file.path(out_dir, paste0(prefix, "_pca_pc_loading_union_signature.pdf"))
  written <- c(written, pt_save_pdf(p, width = fig_w, height = fig_h,
                                    draw = function() {
    pt_draw_clustermap(sub_ord, vmin = -vmax_u, vmax = vmax_u,
                       row_cluster = FALSE, col_cluster = FALSE,
                       row_colors = row_colors,
                       xlabels = x_labels,
                       ylabels = pt_truncate(rownames(sub_ord), 55L),
                       title = sprintf(paste0("%s - union of top +/-%d loadings per ",
                                              "PC (%d pathways) - grouped by top PC + sign"),
                                       tag_title, top_per_pc_per_sign,
                                       length(union_paths)),
                       cbar_label = "loading",
                       legend = list(colours = pc_color,
                                     title = "top PC (max |loading|)"),
                       fontsize_x = 8, fontsize_y = 6, xlabel_rot = 35)
  }))
  written[!vapply(written, is.null, TRUE)]
}


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------
#' Render all 11 PCA PDFs for one PCA result
#'
#' @param res a [pt_pca_result()].
#' @param out_dir directory to write into; created if needed.
#' @param prefix file-name prefix.
#' @param title figure title.
#' @param viz a viz section, from [pt_viz_config()].
#' @return The written paths.
#' @export
pt_render_pca_figures <- function(res, out_dir, prefix, title,
                                  viz = pt_viz_config()) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  scores <- res$scores
  loadings <- res$loadings
  var_expl <- res$variance_explained
  pc_names <- pt_pc_names(res)
  group_labels <- unname(res$labels)
  sub_labels <- if (length(res$sublabels)) unname(res$sublabels) else group_labels
  point_labels <- .pt_short_label(rownames(scores))
  sizes <- if (!is.null(res$sizes)) res$sizes else rep(1, nrow(scores))

  class_lut <- pt_class_colours(sort(unique(group_labels), method = "radix"))
  ct_lut <- pt_celltype_colours(unique(sub_labels))
  primary_name <- .pt_or(res$labels_name, "group")
  secondary_name <- .pt_or(res$sublabels_name, "secondary")

  # Block 1
  .pt_draw_scores_heatmap(
    res, sprintf("%s\nPCA cluster scores (top %d PCs)", title, pt_k(res)),
    file.path(out_dir, paste0(prefix, "_pca_cluster_scores_heatmap.pdf")))

  # Block 2 -- biplots + pairwise, primary colouring
  .pt_draw_biplot(scores, loadings, var_expl, point_labels, group_labels, sizes,
                  class_lut, primary_name,
                  sprintf("%s\nPCA biplot (colour = %s)", title, primary_name),
                  file.path(out_dir, paste0(prefix, "_pca_scatter_pc1_pc2_by_class.pdf")))
  .pt_draw_pairwise(scores, var_expl, pc_names, point_labels, group_labels, sizes,
                    class_lut, primary_name,
                    sprintf("%s - pairwise PCA scatter (colour = %s)", title,
                            primary_name),
                    file.path(out_dir, paste0(prefix, "_pca_scatter_pairwise_by_class.pdf")))
  # Block 2 -- secondary colouring
  .pt_draw_biplot(scores, loadings, var_expl, point_labels, sub_labels, sizes,
                  ct_lut, secondary_name,
                  sprintf("%s\nPCA biplot (colour = %s)", title, secondary_name),
                  file.path(out_dir, paste0(prefix, "_pca_scatter_pc1_pc2_by_celltype.pdf")))
  .pt_draw_pairwise(scores, var_expl, pc_names, point_labels, sub_labels, sizes,
                    ct_lut, secondary_name,
                    sprintf("%s - pairwise PCA scatter (colour = %s)", title,
                            secondary_name),
                    file.path(out_dir, paste0(prefix, "_pca_scatter_pairwise_by_celltype.pdf")))
  # Block 3
  .pt_draw_signature_panels(
    res, sprintf("%s - per-cluster pathway signatures (PCA)", title),
    file.path(out_dir, paste0(prefix, "_pca_per_cluster_signatures.pdf")))
  # Block 4
  .pt_draw_pc_loadings(res, title, out_dir, prefix)
  # Block 5
  .pt_draw_pc_loading_union(res, title, out_dir, prefix)

  # prefixes contain characters that are regex metacharacters (`|` in pseudobulk
  # scope names), so match on plain string boundaries rather than a pattern
  files <- list.files(out_dir)
  hits <- startsWith(files, paste0(prefix, "_pca_")) &
    endsWith(files, ".pdf")
  sort(file.path(out_dir, files[hits]), method = "radix")
}

.pt_nonzero <- function(v) if (!is.finite(v) || v == 0) 1 else v
