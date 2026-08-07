## Cross-language parity check: run every stage in R on the same inputs the
## Python implementation used, and report the agreement stage by stage.
##
## Run r/parity/run_python.py first, then from the repo root:
##
##     Rscript r/parity/compare.R

suppressMessages(pkgload::load_all(file.path("r"), quiet = TRUE))

HERE <- file.path("r", "parity")
DATA <- file.path(HERE, "data")
PY <- file.path(HERE, "python")

read_tsv <- function(name, rownames = FALSE) {
  utils::read.table(file.path(PY, name), sep = "\t", header = TRUE,
                    row.names = if (rownames) 1L else NULL,
                    check.names = FALSE, quote = "\"", comment.char = "",
                    stringsAsFactors = FALSE)
}
read_mat <- function(name) as.matrix(read_tsv(name, rownames = TRUE))

results <- list()
report <- function(stage, detail, max_abs, rel = NA_real_, tol = 1e-9) {
  ok <- is.finite(max_abs) && max_abs <= tol
  results[[length(results) + 1L]] <<- data.frame(
    stage = stage, detail = detail, max_abs_diff = max_abs,
    max_rel_diff = rel, tol = tol, pass = ok, stringsAsFactors = FALSE)
  cat(sprintf("%-4s %-26s %-34s max|diff| = %.3e\n", if (ok) "PASS" else "FAIL",
              stage, detail, max_abs))
  invisible(ok)
}

## Compare two numeric matrices after aligning names.
cmp_mat <- function(stage, detail, r_mat, py_mat, tol = 1e-9) {
  common_r <- intersect(rownames(r_mat), rownames(py_mat))
  common_c <- intersect(colnames(r_mat), colnames(py_mat))
  if (length(common_r) != nrow(py_mat) || length(common_c) != ncol(py_mat)) {
    cat(sprintf("FAIL %-26s %-34s shape/name mismatch: R %dx%d vs py %dx%d\n",
                stage, detail, nrow(r_mat), ncol(r_mat), nrow(py_mat),
                ncol(py_mat)))
    results[[length(results) + 1L]] <<- data.frame(
      stage = stage, detail = paste(detail, "(name mismatch)"),
      max_abs_diff = NA_real_, max_rel_diff = NA_real_, tol = tol,
      pass = FALSE, stringsAsFactors = FALSE)
    return(invisible(FALSE))
  }
  a <- r_mat[common_r, common_c, drop = FALSE]
  b <- py_mat[common_r, common_c, drop = FALSE]
  d <- abs(a - b)
  scale <- pmax(abs(b), 1e-12)
  report(stage, detail, max(d, na.rm = TRUE),
         max(d / scale, na.rm = TRUE), tol)
}

cmp_vec <- function(stage, detail, r_v, py_v, tol = 1e-9) {
  d <- abs(as.double(r_v) - as.double(py_v))
  report(stage, detail, max(d, na.rm = TRUE), NA_real_, tol)
}

## --------------------------------------------------------------------------
cat("\n=== stage 1-2: input + enrichment ===\n")
gmt <- file.path(DATA, "genesets.gmt")
fm <- pt_load_matrix(file.path(DATA, "expr.tsv"),
                     metadata = file.path(DATA, "meta.tsv"))
cat(sprintf("     R feature matrix: %d x %d\n", nrow(fm$data), ncol(fm$data)))

sm <- suppressWarnings(pt_enrich(fm, backend = "ssgsea", gmt = gmt,
                                 geneset = "SYN", verbose = FALSE))
cmp_mat("ssgsea", "NES, continuous matrix", sm$data, read_mat("ssgsea.tsv"),
        tol = 1e-9)

fm_ties <- pt_load_matrix(file.path(DATA, "expr_ties.tsv"),
                          metadata = file.path(DATA, "meta.tsv"))
sm_ties <- suppressWarnings(pt_enrich(fm_ties, backend = "ssgsea", gmt = gmt,
                                      geneset = "SYN", verbose = FALSE))
cmp_mat("ssgsea", "NES, 45% tied (zero) entries", sm_ties$data,
        read_mat("ssgsea_ties.tsv"), tol = 1e-9)

cov_r <- pt_geneset_coverage(pt_features(fm), pt_parse_gmt(gmt), 5L, 5000L)
cov_py <- read_tsv("coverage.tsv")
cov_r <- cov_r[order(cov_r$pathway, method = "radix"), ]
cov_py <- cov_py[order(cov_py$pathway, method = "radix"), ]
same_status <- identical(cov_r$status, cov_py$status) &&
  identical(cov_r$pathway, cov_py$pathway)
report("coverage", "status + n_present per set",
       if (same_status) max(abs(cov_r$n_present - cov_py$n_present)) else Inf)

ora <- suppressWarnings(pt_enrich(fm, backend = "enrichr", gmt = gmt,
                                  geneset = "SYN", top_n = 120L,
                                  verbose = FALSE))
cmp_mat("enrichr", "-log10(adjusted p) per sample", ora$data, read_mat("ora.tsv"),
        tol = 1e-9)

goslim <- suppressWarnings(pt_enrich(fm, backend = "goslim", gmt = gmt,
                                     geneset = "SYN", top_n = 120L,
                                     verbose = FALSE))
cmp_mat("goslim", "membership fraction per sample", goslim$data,
        read_mat("goslim.tsv"), tol = 1e-12)

## --------------------------------------------------------------------------
cat("\n=== stage 3-4: grouping + PCA ===\n")
groups <- pt_group(sm, mode = "target", target_col = "grp")
res <- pt_pca(sm, groups, size_col = "n_cells")[[1L]]

cmp_vec("pca", "variance explained", res$variance_explained,
        read_tsv("pca_variance.tsv")$variance_explained, tol = 1e-10)
cmp_mat("pca", "scores (sign-matched)", res$scores, read_mat("pca_scores.tsv"),
        tol = 1e-9)
cmp_mat("pca", "loadings (sign-matched)", res$loadings,
        read_mat("pca_loadings.tsv"), tol = 1e-9)

sig_r <- res$signatures
sig_py <- read_tsv("pca_signatures.tsv")
key <- function(d) paste(d$cluster, d$rank, d$pathway, sep = "|")
aligned <- identical(key(sig_r), key(sig_py))
report("pca", "signatures (order + scores)",
       if (aligned) max(abs(sig_r$signature_score - sig_py$signature_score))
       else Inf, tol = 1e-9)

groups_scoped <- pt_group(sm, mode = "target", target_col = "grp",
                          scope_col = "sample_id")
scoped <- pt_pca(sm, groups_scoped)
scoped_r <- do.call(rbind, lapply(scoped, function(r)
  data.frame(scope = r$scope, PC = pt_pc_names(r),
             variance_explained = r$variance_explained,
             stringsAsFactors = FALSE)))
scoped_py <- read_tsv("pca_scoped_variance.tsv")
report("pca", "scoped variance explained",
       if (identical(paste(scoped_r$scope, scoped_r$PC),
                     paste(scoped_py$scope, scoped_py$PC)))
         max(abs(scoped_r$variance_explained - scoped_py$variance_explained))
       else Inf, tol = 1e-10)

## --------------------------------------------------------------------------
cat("\n=== stage 5: differential ===\n")
cmp_diff <- function(label, r_res, py_file, tol = 1e-9) {
  a <- r_res$table
  b <- read_tsv(py_file)
  k <- function(d) paste(d$comparison, d$pathway, sep = "|")
  a <- a[match(k(b), k(a)), , drop = FALSE]
  for (col in c("effect", "statistic", "p_value", "fdr")) {
    cmp_vec("diff", sprintf("%s: %s", label, col), a[[col]], b[[col]], tol)
  }
  ok_dir <- identical(as.character(a$direction), as.character(b$direction))
  report("diff", sprintf("%s: direction labels", label),
         if (ok_dir) 0 else Inf)
}

for (m in c("welch", "mannwhitney", "moderated_t")) {
  cmp_diff(m, pt_diff(sm, groups, method = m), sprintf("diff_%s.tsv", m))
}
cmp_diff("welch vs reference group",
         pt_diff(sm, method = "welch", group_col = "arm", reference = "ctrl"),
         "diff_reference.tsv")
cmp_diff("moderated_t explicit contrast",
         pt_diff(sm, method = "moderated_t", group_col = "cell_type",
                 contrasts = list(c("tumor", "T-cell"))),
         "diff_contrast.tsv")

## limma helpers probed directly
lv <- read_tsv("limma_var.tsv")$var
helpers_py <- read_tsv("limma_helpers.tsv")
prior <- pt_fit_fdist(lv, rep(4, length(lv)))
tri <- vapply(c(0.01, 0.1, 0.5, 1.0, 5.0), pt_trigamma_inverse, 1)
cmp_vec("diff", "limma fitFDist + trigammaInverse",
        c(prior$d0, prior$s0_sq, tri), helpers_py$value, tol = 1e-9)

## --------------------------------------------------------------------------
cat("\n=== category roll-up ===\n")
d <- pt_diff(sm, groups, method = "moderated_t")
summ_r <- pt_summarize_categories(d, file.path(DATA, "categories.tsv"),
                                  significance_col = "fdr")
summ_py <- read_tsv("categories.tsv")
k <- function(x) paste(x$comparison, x$significance, x$category, sep = "|")
aligned <- identical(k(summ_r), k(summ_py))
report("categories", "grouped mean effect",
       if (aligned) max(abs(summ_r$mean_effect - summ_py$mean_effect)) else Inf)
report("categories", "counts (n_pathways/up/down)",
       if (aligned) max(abs(c(summ_r$n_pathways - summ_py$n_pathways,
                              summ_r$n_up - summ_py$n_up,
                              summ_r$n_down - summ_py$n_down))) else Inf)

## --------------------------------------------------------------------------
tbl <- do.call(rbind, results)
cat("\n=== summary ===\n")
cat(sprintf("%d of %d checks pass\n", sum(tbl$pass), nrow(tbl)))
failed <- tbl[!tbl$pass, , drop = FALSE]
if (nrow(failed)) {
  cat("\nfailures:\n")
  print(failed[, c("stage", "detail", "max_abs_diff", "tol")], row.names = FALSE)
}
dir.create(file.path(HERE, "r"), recursive = TRUE, showWarnings = FALSE)
utils::write.table(tbl, file.path(HERE, "r", "parity_report.tsv"), sep = "\t",
                   quote = FALSE, row.names = FALSE)
cat("\nreport written to", file.path(HERE, "r", "parity_report.tsv"), "\n")
if (nrow(failed)) quit(status = 1L)
