test_that("defaults match the documented pipeline settings", {
  cfg <- pt_config()
  expect_s3_class(cfg, "pt_config")
  expect_equal(cfg$enrichment$backend, "ssgsea")
  expect_equal(cfg$enrichment$weight, 0.25)
  expect_equal(cfg$enrichment$missing, "as_is")
  expect_equal(cfg$pca$k_max, 10L)
  expect_equal(cfg$pca$top_pcs_per_cluster, 3L)
  expect_false(cfg$diff$enabled)
  expect_true(cfg$viz$make_figures)
  expect_equal(cfg$output_dir, "pathwaytheme_output")
})

test_that("dotted overrides reach nested sections and the top level", {
  cfg <- pt_config_set(pt_config(),
                       "input.matrix_path" = "expr.tsv",
                       "enrichment.top_n" = 50L,
                       output_dir = "results")
  expect_equal(cfg$input$matrix_path, "expr.tsv")
  expect_equal(cfg$enrichment$top_n, 50L)
  expect_equal(cfg$output_dir, "results")
  expect_error(pt_config_set(cfg, "nosuch.field" = 1), "unknown config section")
})

test_that("YAML round-trips, NULLs included", {
  path <- tempfile(fileext = ".yaml")
  cfg <- pt_config_set(pt_config(), "grouping.target_col" = "treatment",
                       "diff.enabled" = TRUE)
  pt_config_to_yaml(cfg, path)
  back <- pt_config_from_yaml(path)
  expect_equal(back$grouping$target_col, "treatment")
  expect_true(back$diff$enabled)
  expect_null(back$grouping$scope_col)
  expect_equal(back$enrichment$weight, cfg$enrichment$weight)
})

test_that("an unknown YAML key is an error, not a silent no-op", {
  path <- tempfile(fileext = ".yaml")
  writeLines(c("enrichment:", "  backend: ssgsea", "  bakcend: typo"), path)
  expect_error(pt_config_from_yaml(path), "unknown key")
})

test_that("a config written by the Python package loads unchanged", {
  path <- tempfile(fileext = ".yaml")
  writeLines(c(
    "input:", "  kind: matrix", "  matrix_path: expr.tsv", "  sep: \"\\t\"",
    "enrichment:", "  backend: ssgsea", "  geneset: GO_BP",
    "  gmt_path: go_bp.gmt", "  min_gene_set_size: 5",
    "grouping:", "  mode: target", "  target_col: treatment",
    "diff:", "  enabled: true", "  method: moderated_t",
    "output_dir: results"), path)
  cfg <- pt_config_from_yaml(path)
  expect_equal(cfg$input$matrix_path, "expr.tsv")
  expect_equal(cfg$enrichment$gmt_path, "go_bp.gmt")
  expect_equal(cfg$diff$method, "moderated_t")
  expect_equal(cfg$output_dir, "results")
  # untouched sections keep their defaults
  expect_equal(cfg$pca$k_max, 10L)
})
