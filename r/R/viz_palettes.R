## Deterministic colour palettes + shared drawing constants, so the same label
## always maps to the same colour (and to the same colour the Python
## implementation picks).

.PT_EXTRAS <- c("#1f77b4", "#2ca02c", "#9467bd", "#8c564b", "#e377c2", "#17becf",
                "#ff7f0e", "#bcbd22", "#7f6d6d", "#3d8b8b", "#c2a44b", "#a4a4a4")

## Known biological cell-type colours (kept for parity; harmless elsewhere)
.PT_CELL_TYPE_COLORS <- c(
  tumor = "#d62728", tumor_proliferating = "#a02020", proliferating = "#ff9896",
  myeloid_proliferating = "#fb7858", low_quality = "#7f7f7f",
  low_numbers = "#bdbdbd", fibroblast = "#8c564b", CAF = "#a0522d",
  myofibroblast = "#c49c94", endothelial = "#1f77b4",
  endothelial_lymphatic = "#4ba6d5", vascular = "#155b8f", pericyte = "#9edae5",
  `endo-peri` = "#17becf", smooth_muscle_cell = "#5b8db5",
  macrophage = "#2ca02c", monocyte = "#5fbb5f", myeloid = "#1f8a1f",
  dendritic_cell = "#8de08d", pDC = "#b9ddb9", mononuclear_phagocyte = "#3d9b3d",
  `T-cell` = "#9467bd", `CD8-T` = "#7e4faf", `B-cell` = "#c5b0d5",
  NK_cell = "#5e3a8a", Treg = "#a986d3", MAIT = "#dac4ec",
  lymphocyte = "#6f3a93", leukocyte = "#b48fd1", skeletal_muscle = "#e377c2",
  muscle_cell = "#f5b6e0", neuronal_glial = "#ff7f0e", neuron = "#ffae5d",
  Schwann_cells = "#d65a05", adipocyte = "#bcbd22", erythrocyte = "#a3070b",
  pneumocyte = "#6a8eb0")

#' Truncate a long term name for display
#'
#' @param term the label.
#' @param n maximum characters, ellipsis included.
#' @return A character vector no longer than `n`.
#' @export
pt_truncate <- function(term, n = 80L) {
  t <- as.character(term)
  ifelse(nchar(t) <= n, t, paste0(substr(t, 1L, n - 3L), "..."))
}

#' Deterministic colour lookup tables for labels
#'
#' `pt_class_colours()` is the palette for the primary comparison labels
#' (tumor red, low_quality grey, everything else from a fixed rotation).
#' `pt_celltype_colours()` is the palette for a secondary label track, reusing
#' known biological colours where the label matches one.
#'
#' Empty labels are left out of the table entirely -- R cannot look a name up by
#' the empty string -- and [pt_label_colours()] resolves them to the grey
#' fallback instead.
#'
#' @param labels character vector of labels.
#' @return A named character vector of hex colours.
#' @export
pt_class_colours <- function(labels) {
  .pt_build_lut(labels, c(tumor = "#d62728", low_quality = "#7f7f7f"))
}

#' @rdname pt_class_colours
#' @export
pt_celltype_colours <- function(labels) {
  .pt_build_lut(labels, .PT_CELL_TYPE_COLORS)
}

.pt_build_lut <- function(labels, known) {
  out <- character()
  k <- 0L
  for (lbl in as.character(labels)) {
    if (is.na(lbl) || !nzchar(lbl) || lbl %in% names(out)) next
    if (lbl %in% names(known)) {
      out[[lbl]] <- unname(known[[lbl]])
    } else {
      out[[lbl]] <- .PT_EXTRAS[[(k %% length(.PT_EXTRAS)) + 1L]]
      k <- k + 1L
    }
  }
  out
}

#' Resolve labels to colours, falling back for unknown and empty labels
#'
#' @param lut a lookup table from [pt_class_colours()] or
#'   [pt_celltype_colours()].
#' @param labels the labels to colour.
#' @param default colour used for labels absent from `lut`, including `NA` and
#'   the empty string.
#' @return An unnamed character vector as long as `labels`.
#' @export
pt_label_colours <- function(lut, labels, default = "#dddddd") {
  out <- unname(lut[as.character(labels)])
  out[is.na(out)] <- default
  out
}

#' Blue-white-red diverging colour ramp
#'
#' The R equivalent of the Python `R_BWR` colormap: linear RGB interpolation
#' between blue, white and red.
#'
#' @param n number of colours to generate.
#' @return A character vector of `n` hex colours.
#' @export
pt_bwr <- function(n = 256L) {
  ramp <- grDevices::colorRamp(PT_BWR_COLORS, space = "rgb")
  rgb <- ramp(seq(0, 1, length.out = n))
  grDevices::rgb(rgb[, 1L], rgb[, 2L], rgb[, 3L], maxColorValue = 255)
}

#' Map values onto a diverging colour scale
#'
#' Values outside `[vmin, vmax]` are clamped to the end colours, matching
#' matplotlib's default `imshow` behaviour.  `NA` becomes light grey.
#'
#' @param x a numeric vector or matrix.
#' @param vmin,vmax the ends of the scale.
#' @param colours the ramp, e.g. from [pt_bwr()].
#' @return A character object shaped like `x`.
#' @export
pt_map_colours <- function(x, vmin, vmax, colours = pt_bwr()) {
  n <- length(colours)
  span <- vmax - vmin
  if (!is.finite(span) || span <= 0) span <- 1
  idx <- floor((as.numeric(x) - vmin) / span * (n - 1L)) + 1L
  idx <- pmin(pmax(idx, 1L), n)
  out <- colours[idx]
  out[is.na(as.numeric(x))] <- "#eeeeee"
  if (is.matrix(x)) out <- matrix(out, nrow = nrow(x), dimnames = dimnames(x))
  out
}

## Matplotlib's tab10, used to colour "home PC" row annotations.
.PT_TAB10 <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
               "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf")

## The two directional colours used throughout the bar figures.
.PT_UP <- "#d62728"
.PT_DOWN <- "#1f77b4"
