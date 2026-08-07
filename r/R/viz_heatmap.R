## A grid-based clustered-heatmap engine.
##
## Six of the package's figures are the same picture with different inputs: a
## rasterised matrix body, optional row/column dendrograms, optional annotation
## colour strips, tick labels, a colour bar and a discrete legend.  Rather than
## depend on a heatmap package, the layout is built once here and reused, which
## keeps the figures consistent and the dependency list short.

#' Open a PDF device, run a drawing function, close the device
#'
#' Retries briefly on the transient permission errors Windows raises when a
#' viewer still holds the previous file open.
#'
#' @param path output file.
#' @param width,height device size in inches.
#' @param draw a function of no arguments that draws onto the open device.
#' @return `path` if the file was written, otherwise `NULL`.
#' @export
pt_save_pdf <- function(path, width, height, draw) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  for (attempt in seq_len(8L)) {
    ok <- tryCatch({
      grDevices::pdf(path, width = width, height = height, onefile = TRUE)
      on.exit(grDevices::dev.off(), add = FALSE)
      draw()
      grDevices::dev.off()
      on.exit()
      TRUE
    }, error = function(e) {
      try(grDevices::dev.off(), silent = TRUE)
      FALSE
    })
    if (ok) return(path)
    Sys.sleep(2)
  }
  message("  SAVE FAILED: ", path)
  NULL
}

## ---- dendrograms -----------------------------------------------------------

## Cluster the rows (or columns) of a matrix the way seaborn's clustermap does:
## euclidean distance, average linkage.
.pt_cluster <- function(m, method = "average") {
  if (nrow(m) < 2L) return(NULL)
  d <- stats::dist(m, method = "euclidean")
  if (!all(is.finite(d))) return(NULL)
  stats::hclust(d, method = method)
}

## Line segments of a dendrogram in (leaf position, height) coordinates.
.pt_dendro_segments <- function(hc) {
  n <- length(hc$order)
  pos <- numeric(n)
  pos[hc$order] <- seq_len(n)
  xs <- numeric(nrow(hc$merge))
  ys <- hc$height
  segs <- list()
  child_xy <- function(k) {
    if (k < 0) c(pos[[-k]], 0) else c(xs[[k]], ys[[k]])
  }
  for (i in seq_len(nrow(hc$merge))) {
    a <- child_xy(hc$merge[i, 1L])
    b <- child_xy(hc$merge[i, 2L])
    xs[[i]] <- (a[[1L]] + b[[1L]]) / 2
    y <- ys[[i]]
    segs[[length(segs) + 1L]] <- c(a[[1L]], a[[2L]], a[[1L]], y)
    segs[[length(segs) + 1L]] <- c(b[[1L]], b[[2L]], b[[1L]], y)
    segs[[length(segs) + 1L]] <- c(a[[1L]], y, b[[1L]], y)
  }
  m <- do.call(rbind, segs)
  colnames(m) <- c("x0", "y0", "x1", "y1")
  as.data.frame(m)
}

## Draw a dendrogram inside the current viewport.  `side` is "top" (leaves along
## x, height growing upward) or "left" (leaves along y, height growing leftward).
.pt_draw_dendro <- function(hc, n, side = "top") {
  if (is.null(hc)) return(invisible(NULL))
  segs <- .pt_dendro_segments(hc)
  hmax <- max(hc$height, na.rm = TRUE)
  if (!is.finite(hmax) || hmax <= 0) hmax <- 1
  leaf <- function(v) (v - 0.5) / n
  hgt <- function(v) v / hmax
  if (identical(side, "top")) {
    grid::grid.segments(leaf(segs$x0), hgt(segs$y0), leaf(segs$x1), hgt(segs$y1),
                        gp = grid::gpar(col = "#555555", lwd = 0.6))
  } else {
    # rows are drawn top-to-bottom, so leaf 1 sits at the top
    grid::grid.segments(1 - hgt(segs$y0), 1 - leaf(segs$x0),
                        1 - hgt(segs$y1), 1 - leaf(segs$x1),
                        gp = grid::gpar(col = "#555555", lwd = 0.6))
  }
}

## ---- the heatmap itself ----------------------------------------------------

#' Draw a clustered heatmap
#'
#' @param mat numeric matrix; rows are drawn top to bottom.
#' @param vmin,vmax colour-scale limits.  Defaults to a symmetric scale at the
#'   98th percentile of `abs(mat)`.
#' @param row_cluster,col_cluster cluster and reorder rows / columns.
#' @param method linkage method passed to [stats::hclust()].
#' @param row_colors,col_colors optional annotation colour vectors, one entry
#'   per row / column *before* clustering reorders them.
#' @param xlabels,ylabels tick labels; `NULL` uses the dimnames, `NA` hides them.
#' @param title figure title, drawn top-left.
#' @param cbar_label label under the colour bar.
#' @param legend optional `list(colours = <named vector>, title = )` drawn as
#'   swatches at the top left.
#' @param fontsize_x,fontsize_y,fontsize_title point sizes.
#' @param colours the colour ramp.
#' @param xlabel_rot rotation for the column labels, in degrees.
#' @return The row and column orders, invisibly.
#' @export
pt_draw_clustermap <- function(mat, vmin = NULL, vmax = NULL,
                               row_cluster = TRUE, col_cluster = FALSE,
                               method = "average",
                               row_colors = NULL, col_colors = NULL,
                               xlabels = NULL, ylabels = NULL,
                               title = NULL, cbar_label = NULL,
                               legend = NULL,
                               fontsize_x = 7, fontsize_y = 6,
                               fontsize_title = 10,
                               colours = pt_bwr(), xlabel_rot = 0) {
  n_row <- nrow(mat)
  n_col <- ncol(mat)
  if (is.null(vmax)) {
    vmax <- stats::quantile(abs(mat), 0.98, na.rm = TRUE, names = FALSE)
    if (!is.finite(vmax) || vmax <= 0) vmax <- 1
  }
  if (is.null(vmin)) vmin <- -vmax

  row_hc <- if (isTRUE(row_cluster)) .pt_cluster(mat, method) else NULL
  col_hc <- if (isTRUE(col_cluster)) .pt_cluster(t(mat), method) else NULL
  row_ord <- if (is.null(row_hc)) seq_len(n_row) else row_hc$order
  col_ord <- if (is.null(col_hc)) seq_len(n_col) else col_hc$order

  body <- mat[row_ord, col_ord, drop = FALSE]
  if (is.null(xlabels)) xlabels <- colnames(mat)
  if (is.null(ylabels)) ylabels <- rownames(mat)
  show_y <- !identical(ylabels, NA) && length(ylabels) == n_row
  show_x <- !identical(xlabels, NA) && length(xlabels) == n_col
  xlab_ord <- if (show_x) xlabels[col_ord] else NULL
  ylab_ord <- if (show_y) ylabels[row_ord] else NULL

  # ---- layout ----
  left_pad <- grid::unit(if (is.null(legend) && is.null(cbar_label)) 0.4 else 1.5,
                         "inches")
  row_dend_w <- grid::unit(if (is.null(row_hc)) 0 else 0.5, "inches")
  row_ann_w <- grid::unit(if (is.null(row_colors)) 0 else 0.14, "inches")
  y_lab_w <- if (show_y)
    grid::unit(1, "strwidth", list(.pt_widest(ylab_ord))) + grid::unit(2, "mm")
  else grid::unit(0.15, "inches")

  title_h <- grid::unit(if (is.null(title)) 0.15 else
    0.28 * (1 + .pt_n_lines(title)), "inches")
  col_dend_h <- grid::unit(if (is.null(col_hc)) 0 else 0.45, "inches")
  col_ann_h <- grid::unit(if (is.null(col_colors)) 0 else 0.14, "inches")
  x_lab_h <- if (show_x)
    grid::unit(0.13 * max(1, .pt_n_lines(xlab_ord)) +
                 (if (xlabel_rot != 0) 0.55 else 0.05), "inches")
  else grid::unit(0.1, "inches")

  grid::grid.newpage()
  lay <- grid::grid.layout(
    nrow = 5L, ncol = 5L,
    widths = grid::unit.c(left_pad, row_dend_w, row_ann_w,
                          grid::unit(1, "null"), y_lab_w),
    heights = grid::unit.c(title_h, col_dend_h, col_ann_h,
                           grid::unit(1, "null"), x_lab_h))
  grid::pushViewport(grid::viewport(layout = lay,
                                    width = grid::unit(1, "npc") - grid::unit(4, "mm"),
                                    height = grid::unit(1, "npc") - grid::unit(4, "mm")))

  vp <- function(r, c) grid::viewport(layout.pos.row = r, layout.pos.col = c)

  # title
  if (!is.null(title)) {
    grid::pushViewport(vp(1L, 4L))
    grid::grid.text(title, x = 0, y = 0.5, just = c("left", "centre"),
                    gp = grid::gpar(fontsize = fontsize_title))
    grid::popViewport()
  }

  # dendrograms
  if (!is.null(col_hc)) {
    grid::pushViewport(vp(2L, 4L))
    .pt_draw_dendro(col_hc, n_col, side = "top")
    grid::popViewport()
  }
  if (!is.null(row_hc)) {
    grid::pushViewport(vp(4L, 2L))
    .pt_draw_dendro(row_hc, n_row, side = "left")
    grid::popViewport()
  }

  # annotation strips
  if (!is.null(col_colors)) {
    grid::pushViewport(vp(3L, 4L))
    grid::grid.raster(matrix(col_colors[col_ord], nrow = 1L),
                      width = 1, height = 1, interpolate = FALSE)
    grid::popViewport()
  }
  if (!is.null(row_colors)) {
    grid::pushViewport(vp(4L, 3L))
    grid::grid.raster(matrix(row_colors[row_ord], ncol = 1L),
                      width = 1, height = 1, interpolate = FALSE)
    grid::popViewport()
  }

  # body
  grid::pushViewport(vp(4L, 4L))
  grid::grid.raster(pt_map_colours(body, vmin, vmax, colours),
                    width = 1, height = 1, interpolate = FALSE)
  grid::grid.rect(gp = grid::gpar(col = "#888888", fill = NA, lwd = 0.5))
  grid::popViewport()

  # tick labels
  if (show_y) {
    grid::pushViewport(vp(4L, 5L))
    grid::grid.text(ylab_ord, x = grid::unit(1, "mm"),
                    y = grid::unit(1 - (seq_len(n_row) - 0.5) / n_row, "npc"),
                    just = c("left", "centre"),
                    gp = grid::gpar(fontsize = fontsize_y))
    grid::popViewport()
  }
  if (show_x) {
    grid::pushViewport(vp(5L, 4L))
    grid::grid.text(xlab_ord,
                    x = grid::unit((seq_len(n_col) - 0.5) / n_col, "npc"),
                    y = grid::unit(1, "npc") - grid::unit(1, "mm"),
                    rot = xlabel_rot,
                    just = if (xlabel_rot != 0) c("right", "centre") else c("centre", "top"),
                    gp = grid::gpar(fontsize = fontsize_x))
    grid::popViewport()
  }

  # colour bar + legend live in the left gutter
  grid::pushViewport(vp(4L, 1L))
  .pt_draw_colourbar(vmin, vmax, colours, cbar_label)
  grid::popViewport()
  if (!is.null(legend)) {
    grid::pushViewport(grid::viewport(layout.pos.row = 4L, layout.pos.col = 1L))
    .pt_draw_legend(legend$colours, legend$title, y_top = 0.62)
    grid::popViewport()
  }

  grid::popViewport()
  invisible(list(row_order = row_ord, col_order = col_ord))
}

.pt_draw_colourbar <- function(vmin, vmax, colours, label, y_top = 0.98,
                               height = 0.16) {
  grid::pushViewport(grid::viewport(x = grid::unit(2, "mm"),
                                    y = grid::unit(y_top, "npc"),
                                    width = grid::unit(4, "mm"),
                                    height = grid::unit(height, "npc"),
                                    just = c("left", "top")))
  grid::grid.raster(matrix(rev(colours), ncol = 1L), width = 1, height = 1,
                    interpolate = TRUE)
  grid::grid.rect(gp = grid::gpar(col = "#666666", fill = NA, lwd = 0.4))
  grid::grid.text(.pt_fmt(vmax), x = grid::unit(1, "npc") + grid::unit(1, "mm"),
                  y = 1, just = c("left", "centre"), gp = grid::gpar(fontsize = 6))
  grid::grid.text(.pt_fmt(vmin), x = grid::unit(1, "npc") + grid::unit(1, "mm"),
                  y = 0, just = c("left", "centre"), gp = grid::gpar(fontsize = 6))
  if (!is.null(label))
    grid::grid.text(label, x = 0, y = grid::unit(-2, "mm"),
                    just = c("left", "top"), gp = grid::gpar(fontsize = 7))
  grid::popViewport()
}

.pt_draw_legend <- function(colours, title = NULL, y_top = 0.6,
                            fontsize = 7, max_items = 24L) {
  if (!length(colours)) return(invisible(NULL))
  keys <- names(colours)
  if (length(keys) > max_items) {
    colours <- colours[seq_len(max_items)]
    keys <- names(colours)
  }
  step <- grid::unit(3.6, "mm")
  y <- grid::unit(y_top, "npc")
  if (!is.null(title)) {
    grid::grid.text(title, x = grid::unit(2, "mm"), y = y,
                    just = c("left", "centre"),
                    gp = grid::gpar(fontsize = fontsize + 1, fontface = "bold"))
    y <- y - step
  }
  for (i in seq_along(colours)) {
    grid::grid.rect(x = grid::unit(2, "mm"), y = y, width = grid::unit(2.6, "mm"),
                    height = grid::unit(2.6, "mm"), just = c("left", "centre"),
                    gp = grid::gpar(fill = colours[[i]], col = NA))
    grid::grid.text(if (nzchar(keys[[i]])) keys[[i]] else "(none)",
                    x = grid::unit(5.6, "mm"), y = y, just = c("left", "centre"),
                    gp = grid::gpar(fontsize = fontsize))
    y <- y - step
  }
  invisible(NULL)
}

## ---- small helpers ---------------------------------------------------------

.pt_widest <- function(x) {
  if (!length(x)) return("")
  x[[which.max(nchar(as.character(x)))]]
}

.pt_n_lines <- function(x) {
  max(vapply(strsplit(as.character(x), "\n", fixed = TRUE), length, 1L), 1L)
}

.pt_fmt <- function(v) formatC(v, digits = 2L, format = "g")
