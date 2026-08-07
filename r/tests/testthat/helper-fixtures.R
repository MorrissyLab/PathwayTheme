## Shared synthetic fixtures -- no external data or network needed.

## A tiny GO-style .gmt: 15 terms, each 8-18 genes drawn from 200 genes.
synthetic_gmt <- function(genes = sprintf("G%03d", 0:199), seed = 0L) {
  set.seed(seed)
  path <- tempfile(fileext = ".gmt")
  lines <- vapply(0:14, function(t) {
    members <- sample(genes, sample(8:17, 1L))
    paste(c(sprintf("TERM_%02d", t), "desc", members), collapse = "\t")
  }, "")
  writeLines(lines, path)
  list(path = path, genes = genes)
}

## 200 genes x 12 samples, with two latent groups baked in for structure.
synthetic_feature_matrix <- function(genes = sprintf("G%03d", 0:199),
                                     seed = 1L) {
  set.seed(seed)
  cols <- sprintf("S%02d", 0:11)
  x <- matrix(rnorm(length(genes) * 12L), nrow = length(genes),
              dimnames = list(genes, cols))
  x[1:40, 1:6] <- x[1:40, 1:6] + 2
  meta <- data.frame(
    sample_id = rep(c("A", "B"), each = 6L),
    cluster_id = as.character(rep(0:5, 2L)),
    grp = rep(c("hi", "lo"), each = 6L),
    cell_type = rep(c("tumor", "macrophage", "T-cell"), 4L),
    n_cells = sample(20:500, 12L),
    row.names = cols, stringsAsFactors = FALSE)
  pt_feature_matrix(x, pt_sample_metadata(meta))
}

## 40 terms x 12 samples continuous score matrix + metadata.
synthetic_score_matrix <- function(seed = 2L) {
  set.seed(seed)
  cols <- sprintf("S%02d", 0:11)
  x <- matrix(rnorm(40L * 12L), nrow = 40L,
              dimnames = list(sprintf("TERM_%02d", 0:39), cols))
  x[1:10, 1:6] <- x[1:10, 1:6] + 3
  meta <- data.frame(
    sample_id = rep(c("A", "B"), each = 6L),
    cluster_id = as.character(rep(0:5, 2L)),
    grp = rep(c("hi", "lo"), each = 6L),
    cell_type = rep(c("tumor", "macrophage", "T-cell"), 4L),
    n_cells = sample(20:500, 12L),
    row.names = cols, stringsAsFactors = FALSE)
  pt_score_matrix(x, pt_sample_metadata(meta), backend = "synthetic",
                  geneset = "SYN")
}
