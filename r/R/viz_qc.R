## Full-matrix sanity heatmap -- a QC view of the whole score matrix.
##
## Z-scores each pathway across samples (population sd) and draws a clustered
## heatmap, optionally with a colour strip for one metadata label.  Very large
## matrices are capped to the most-variable pathways so the figure stays
## legible, and the cap is stated in the title.

#' Write a z-scored pathway x sample QC heatmap
#'
#' @param sm a [pt_score_matrix()].
#' @param out_path output PDF path.
#' @param label_col optional metadata column driving the column colour strip.
#' @param max_pathways cap on the most-variable pathways (0 keeps all).
#' @param title heatmap title.
#' @return The written path, or `NULL` when the matrix is too small to draw.
#' @export
pt_render_sanity_heatmap <- function(sm, out_path, label_col = NULL,
                                     max_pathways = 200L,
                                     title = "Pathway score matrix (z-scored)") {
  x <- sm$data
  x[is.na(x)] <- 0
  sds <- sqrt(rowMeans((x - rowMeans(x))^2))
  x <- x[sds > 1e-12, , drop = FALSE]
  if (nrow(x) < 2L || ncol(x) < 2L) return(NULL)
  n_variable <- nrow(x)

  if (max_pathways > 0L && nrow(x) > max_pathways) {
    v <- sqrt(rowMeans((x - rowMeans(x))^2))
    keep <- utils::head(order(v, decreasing = TRUE, method = "radix"), max_pathways)
    x <- x[keep, , drop = FALSE]
    title <- sprintf("%s \\u2014 top %d most-variable of %d", title,
                     max_pathways,
                     n_variable)
  }

  mu <- rowMeans(x)
  sdv <- sqrt(rowMeans((x - mu)^2))
  z <- (x - mu) / sdv
  rownames(z) <- pt_truncate(rownames(z), 60L)

  col_colors <- NULL
  lut <- NULL
  if (!is.null(label_col) && pt_meta_has(sm$metadata, label_col)) {
    labs <- as.character(pt_meta_get(sm$metadata, label_col))[
      match(pt_samples(sm), rownames(sm$metadata$table))]
    lut <- pt_class_colours(sort(unique(labs), method = "radix"))
    col_colors <- pt_label_colours(lut, labs)
  }

  h <- min(40, max(4, 0.16 * nrow(z) + 2))
  w <- min(30, max(6, 0.18 * ncol(z) + 3))
  show_y <- if (nrow(z) <= 120L) rownames(z) else NA

  pt_save_pdf(out_path, width = w + 1.5, height = h, draw = function() {
    pt_draw_clustermap(
      z, vmin = -2, vmax = 2, row_cluster = TRUE, col_cluster = TRUE,
      col_colors = col_colors, xlabels = colnames(z), ylabels = show_y,
      title = title, cbar_label = "z-score",
      legend = if (is.null(lut)) NULL else
        list(colours = lut[sort(names(lut), method = "radix")],
             title = label_col),
      fontsize_x = 6, fontsize_y = 5, xlabel_rot = 90)
  })
}
