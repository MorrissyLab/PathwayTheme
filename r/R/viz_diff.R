## Figures for the differential pathway stage: one volcano + one top-pathway
## bar panel per comparison.

.PT_NEGLOG_CAP <- 320   # -log10(p) cap so p == 0 rows stay on-plot

.pt_neglog10 <- function(p) {
  v <- -log10(as.double(p))
  v[!is.finite(v)] <- .PT_NEGLOG_CAP
  pmin(pmax(v, 0), .PT_NEGLOG_CAP)
}

.pt_draw_volcano <- function(df, comparison, out_path, alpha = 0.05) {
  d <- df[!is.na(df$effect) & !is.na(df$fdr), , drop = FALSE]
  if (!nrow(d)) return(NULL)
  d$neglog <- .pt_neglog10(d$fdr)
  sig <- d$fdr < alpha
  d$sig <- ifelse(!sig, sprintf("FDR \\u2265 %g", alpha),
                  ifelse(d$effect >= 0, "up", "down"))
  d$sig <- factor(d$sig, levels = c(sprintf("FDR \\u2265 %g", alpha), "up", "down"))
  pal <- stats::setNames(c("#bbbbbb", .PT_UP, .PT_DOWN), levels(d$sig))

  lab <- d[sig, , drop = FALSE]
  lab <- lab[order(lab$fdr, method = "radix"), , drop = FALSE]
  lab <- utils::head(lab, 8L)
  lab$label <- pt_truncate(lab$pathway, 40L)

  p <- ggplot2::ggplot(d, ggplot2::aes(x = effect, y = neglog, colour = sig)) +
    ggplot2::geom_vline(xintercept = 0, colour = "#bbbbbb", linewidth = 0.25) +
    ggplot2::geom_hline(yintercept = -log10(alpha), colour = "#888888",
                        linewidth = 0.3, linetype = "dashed") +
    ggplot2::geom_point(size = 0.9, alpha = 0.75) +
    ggplot2::scale_colour_manual(values = pal, name = NULL, drop = FALSE) +
    ggplot2::labs(x = "effect (mean_case \\u2212 mean_reference)",
                  y = "\\u2212log10(FDR)", title = comparison) +
    .pt_blank_theme(10) +
    ggplot2::theme(legend.position = c(0.99, 0.99),
                   legend.justification = c(1, 1),
                   legend.background = ggplot2::element_blank(),
                   legend.text = ggplot2::element_text(size = 7))
  if (nrow(lab))
    p <- p + ggrepel::geom_text_repel(
      data = lab, ggplot2::aes(label = label), inherit.aes = TRUE,
      size = 2, colour = "#333333", show.legend = FALSE, max.overlaps = Inf,
      segment.colour = "#cccccc", segment.size = 0.15)

  pt_save_pdf(out_path, width = 7, height = 6, draw = function() print(p))
}

.pt_draw_top_bars <- function(df, comparison, out_path, top_n = 25L) {
  d <- df[!is.na(df$effect), , drop = FALSE]
  if (!nrow(d)) return(NULL)
  d <- d[order(d$fdr, na.last = TRUE, method = "radix"), , drop = FALSE]
  d <- utils::head(d, top_n)
  d <- d[order(d$effect, method = "radix"), , drop = FALSE]
  d$row_i <- seq_len(nrow(d))
  d$fill <- ifelse(d$effect >= 0, .PT_UP, .PT_DOWN)

  p <- ggplot2::ggplot(d, ggplot2::aes(x = effect, y = row_i, fill = fill)) +
    ggplot2::geom_col(width = 0.8, orientation = "y") +
    ggplot2::geom_vline(xintercept = 0, colour = "#333333", linewidth = 0.25) +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_y_continuous(breaks = d$row_i,
                                labels = pt_truncate(d$pathway, 60L),
                                expand = ggplot2::expansion(add = 0.6)) +
    ggplot2::labs(x = "effect (mean_case \\u2212 mean_reference)", y = NULL,
                  title = sprintf("%s \\u2014 top %d by FDR", comparison,
                                  nrow(d))) +
    .pt_blank_theme(10) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 6))

  pt_save_pdf(out_path, width = 9, height = max(3, 0.32 * nrow(d) + 1),
              draw = function() print(p))
}

#' Volcano and top-pathway figures for a differential result
#'
#' One volcano plus one top-pathway bar panel per comparison.
#'
#' @param result a [pt_diff_result()].
#' @param out_dir directory to write into.
#' @param prefix file-name prefix.
#' @param top_n pathways shown in the bar panel.
#' @return The written paths.
#' @export
pt_render_diff_figures <- function(result, out_dir, prefix, top_n = 25L) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  written <- character()
  for (comp in pt_comparisons(result)) {
    sub <- result$table[result$table$comparison == comp, , drop = FALSE]
    safe <- gsub("[^A-Za-z0-9._-]", "_", comp)
    vp <- file.path(out_dir, paste0(prefix, "_diff_", safe, "_volcano.pdf"))
    bp <- file.path(out_dir, paste0(prefix, "_diff_", safe, "_top.pdf"))
    .pt_draw_volcano(sub, comp, vp)
    .pt_draw_top_bars(sub, comp, bp, top_n = top_n)
    written <- c(written, Filter(file.exists, c(vp, bp)))
  }
  written
}
