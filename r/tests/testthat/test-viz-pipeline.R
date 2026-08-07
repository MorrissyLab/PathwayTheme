test_that("the same label always gets the same colour", {
  lut <- pt_class_colours(c("tumor", "alpha", "beta", "low_quality"))
  expect_equal(unname(lut[["tumor"]]), "#d62728")
  expect_equal(unname(lut[["low_quality"]]), "#7f7f7f")
  # unknown labels take the rotation in order of first appearance
  expect_equal(unname(lut[["alpha"]]), "#1f77b4")
  expect_equal(unname(lut[["beta"]]), "#2ca02c")
  expect_equal(pt_class_colours(c("alpha", "beta")), lut[c("alpha", "beta")])

  ct <- pt_celltype_colours(c("macrophage", "unknown_type", ""))
  expect_equal(unname(ct[["macrophage"]]), "#2ca02c")
  # the empty label is never a key -- R cannot look one up by "" -- and is
  # resolved to the fallback at draw time instead
  expect_false("" %in% names(ct))
  expect_equal(pt_label_colours(ct, c("macrophage", "", NA, "nope")),
               c("#2ca02c", "#dddddd", "#dddddd", "#dddddd"))
})

test_that("truncation keeps short labels and ellipsises long ones", {
  expect_equal(pt_truncate("short", 10L), "short")
  expect_equal(pt_truncate(strrep("x", 20L), 10L), paste0(strrep("x", 7L), "..."))
  expect_length(pt_truncate(c("a", strrep("b", 30L)), 10L), 2L)
})

test_that("the colour ramp runs blue to white to red and clamps", {
  ramp <- pt_bwr(3L)
  expect_equal(toupper(substr(ramp, 1L, 7L)), c("#0000FF", "#FFFFFF", "#FF0000"))
  cols <- pt_map_colours(c(-99, 0, 99, NA), vmin = -1, vmax = 1, colours = ramp)
  expect_equal(toupper(cols[1:3]), c("#0000FF", "#FFFFFF", "#FF0000"))
  expect_equal(cols[[4L]], "#eeeeee")
})

test_that("PCA output is the documented 11 figures and 3 tables", {
  res <- pt_pca(synthetic_score_matrix(), k_max = 4L)[[1L]]
  out <- file.path(tempdir(), "pt_viz_test")
  unlink(out, recursive = TRUE)

  written <- pt_figures(res, out, geneset = "SYN")[["__all__"]]
  pdfs <- grep("\\.pdf$", written, value = TRUE)
  tsvs <- grep("\\.tsv$", written, value = TRUE)
  expect_length(pdfs, 11L)
  expect_length(tsvs, 3L)
  expect_true(all(file.exists(written)))
  expect_true(all(file.size(pdfs) > 1000))

  expect_setequal(basename(tsvs), paste0("SYN___all___pca_", c(
    "variance_explained.tsv", "top_pathways_per_pc.tsv",
    "per_cluster_signatures.tsv")))
  expect_setequal(basename(pdfs), paste0("SYN___all___pca_", c(
    "cluster_scores_heatmap.pdf", "scatter_pc1_pc2_by_class.pdf",
    "scatter_pairwise_by_class.pdf", "scatter_pc1_pc2_by_celltype.pdf",
    "scatter_pairwise_by_celltype.pdf", "per_cluster_signatures.pdf",
    "pc_loadings.pdf", "pathway_pc_heatmap.pdf", "pathway_pc_clustermap.pdf",
    "pc_loading_union_clustermap.pdf", "pc_loading_union_signature.pdf")))
})

test_that("the variance table carries a cumulative column", {
  res <- pt_pca(synthetic_score_matrix(), k_max = 4L)[[1L]]
  out <- file.path(tempdir(), "pt_tab_test")
  unlink(out, recursive = TRUE)
  paths <- pt_write_pca_tables(res, out, "SYN", pt_pca_config())
  tab <- utils::read.table(grep("variance", paths, value = TRUE), sep = "\t",
                           header = TRUE)
  expect_equal(tab$variance_explained, res$variance_explained, tolerance = 1e-12)
  expect_equal(tab$cumulative_variance, cumsum(res$variance_explained),
               tolerance = 1e-12)
})

test_that("differential figures are written per comparison", {
  sm <- synthetic_score_matrix()
  g <- pt_group(sm, mode = "target", target_col = "grp")
  d <- pt_diff(sm, g, method = "welch")
  out <- file.path(tempdir(), "pt_diff_viz")
  unlink(out, recursive = TRUE)
  paths <- pt_render_diff_figures(d, out, "SYN")
  expect_length(paths, 2L * length(pt_comparisons(d)))
  expect_true(all(file.exists(paths)))
})

test_that("the QC heatmap is written, and skipped when too small", {
  sm <- synthetic_score_matrix()
  p <- pt_sanity_heatmap(sm, file.path(tempdir(), "pt_qc.pdf"), label_col = "grp")
  expect_true(file.exists(p))
  tiny <- pt_subset_samples(sm, pt_samples(sm)[[1L]])
  expect_null(pt_sanity_heatmap(tiny, file.path(tempdir(), "pt_qc2.pdf")))
})

test_that("the full pipeline runs end to end from a config", {
  fx <- synthetic_gmt()
  fm <- synthetic_feature_matrix()
  expr_path <- tempfile(fileext = ".tsv")
  meta_path <- tempfile(fileext = ".tsv")
  utils::write.table(fm$data, expr_path, sep = "\t", quote = FALSE, col.names = NA)
  utils::write.table(fm$metadata$table, meta_path, sep = "\t", quote = FALSE,
                     col.names = NA)
  out <- file.path(tempdir(), "pt_pipeline_test")
  unlink(out, recursive = TRUE)

  res <- suppressWarnings(suppressMessages(pt_run_pipeline(
    "input.matrix_path" = expr_path,
    "input.metadata_path" = meta_path,
    "enrichment.backend" = "ssgsea",
    "enrichment.gmt_path" = fx$path,
    "enrichment.geneset" = "SYN",
    "grouping.mode" = "target",
    "grouping.target_col" = "grp",
    "grouping.secondary_col" = "cell_type",
    "grouping.size_col" = "n_cells",
    "diff.enabled" = TRUE,
    "diff.method" = "moderated_t",
    "viz.sanity_heatmap" = TRUE,
    output_dir = out, verbose = FALSE)))

  expect_s3_class(res, "pt_pipeline_result")
  expect_s3_class(res$score_matrix, "pt_score_matrix")
  expect_length(res$pca_results, 1L)
  expect_s3_class(res$diff_result, "pt_diff_result")

  scope_dir <- file.path(out, "SYN", "sample_pca", "__all__")
  expect_length(list.files(scope_dir, pattern = "\\.pdf$"), 11L)
  expect_length(list.files(scope_dir, pattern = "\\.tsv$"), 3L)
  expect_true(file.exists(file.path(out, "SYN", "SYN_geneset_coverage.tsv")))
  expect_true(file.exists(file.path(out, "SYN", "SYN_sanity_heatmap.pdf")))
  expect_true(file.exists(file.path(out, "SYN",
                                    "SYN_differential_moderated_t.tsv")))
  # display labels default to cl<cluster_id> when the column is present
  expect_true(all(grepl("^cl", rownames(res$pca_results[[1L]]$scores))))
})

test_that("the CLI writes a config that loads back", {
  path <- tempfile(fileext = ".yaml")
  expect_output(pt_main(c("init", "-o", path)), "wrote default config")
  expect_s3_class(pt_config_from_yaml(path), "pt_config")
  expect_output(pt_main(character()), "usage:")
  expect_output(pt_main("nosuch"), "usage:")
})
