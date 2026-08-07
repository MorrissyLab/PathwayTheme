test_that("target mode reads labels straight from the metadata column", {
  sm <- synthetic_score_matrix()
  g <- pt_group(sm, mode = "target", target_col = "grp")
  expect_identical(unname(g$labels), sm$metadata$table$grp)
  expect_identical(names(g$labels), pt_samples(sm))
  expect_identical(g$label_name, "grp")
  expect_null(g$scope)
})

test_that("a missing target column is an error naming the alternatives", {
  sm <- synthetic_score_matrix()
  expect_error(pt_group(sm, mode = "target", target_col = "nope"), "not in metadata")
  expect_error(pt_group(sm, mode = "target"), "requires grouping\\$target_col")
  expect_error(pt_group(sm, mode = "sideways", target_col = "grp"),
               "Unknown grouping\\$mode")
  expect_error(pt_group(sm, mode = "target", target_col = "grp",
                        scope_col = "nope"), "scope_col")
})

test_that("auto clustering recovers a planted two-group split", {
  sm <- synthetic_score_matrix()   # first 6 samples are shifted
  g <- pt_group(sm, mode = "auto", n_clusters = 2L)
  labels <- g$labels
  expect_length(unique(labels), 2L)
  # the planted groups must not be split across clusters
  expect_length(unique(labels[1:6]), 1L)
  expect_length(unique(labels[7:12]), 1L)
  expect_false(labels[[1L]] == labels[[7L]])
})

test_that("hierarchical clustering is available and also finds the split", {
  sm <- synthetic_score_matrix()
  g <- pt_group(sm, mode = "auto", method = "hierarchical", n_clusters = 2L)
  expect_length(unique(g$labels[1:6]), 1L)
  expect_length(unique(g$labels[7:12]), 1L)
})

test_that("silhouette selection picks a k without being told one", {
  sm <- synthetic_score_matrix()
  g <- pt_group(sm, mode = "auto", k_min = 2L, k_max = 5L)
  expect_gte(length(unique(g$labels)), 2L)
})

test_that("silhouette width behaves on an obvious two-cluster layout", {
  x <- rbind(matrix(rnorm(20, 0), ncol = 2), matrix(rnorm(20, 20), ncol = 2))
  good <- pt_silhouette(x, rep(c("a", "b"), each = 10L))
  bad <- pt_silhouette(x, rep(c("a", "b"), times = 10L))
  expect_gt(good, 0.9)
  expect_lt(bad, good)
  expect_true(is.na(pt_silhouette(x, rep("a", 20L))))
})

test_that("filtering keeps the chosen samples and their metadata", {
  fm <- synthetic_feature_matrix()
  kept <- pt_filter_samples(fm, "grp", "hi")
  expect_equal(ncol(kept$data), 6L)
  expect_true(all(kept$metadata$table$grp == "hi"))
  expect_error(pt_filter_samples(fm, "grp", "nope"), "no samples have")
  expect_error(pt_filter_samples(fm, "nope", "hi"), "not found")
})


# ---- category roll-up ------------------------------------------------------

test_that("categories aggregate the requested statistic with counts", {
  tab <- data.frame(
    pathway = c("P1", "P2", "P3", "P4"),
    effect = c(2, 4, -1, -3),
    fdr = c(0.01, 0.2, 0.001, 0.9),
    stringsAsFactors = FALSE)
  map <- c(P1 = "alpha", P2 = "alpha", P3 = "beta")

  summ <- pt_summarize_by_category(tab, map)
  expect_setequal(summ$category, c("alpha", "beta", "Other"))
  alpha <- summ[summ$category == "alpha", ]
  expect_equal(alpha$mean_effect, 3)
  expect_equal(alpha$n_pathways, 2L)
  expect_equal(alpha$n_up, 2L)
  expect_equal(alpha$n_down, 0L)
  # P4 is unmapped and lands in the fallback category
  expect_equal(summ$n_pathways[summ$category == "Other"], 1L)
})

test_that("the significance split uses alpha and treats NA as non-significant", {
  tab <- data.frame(pathway = c("P1", "P2", "P3"), effect = c(1, 2, 3),
                    fdr = c(0.01, 0.2, NA_real_), stringsAsFactors = FALSE)
  map <- c(P1 = "alpha", P2 = "alpha", P3 = "alpha")
  summ <- pt_summarize_by_category(tab, map, significance_col = "fdr",
                                   alpha = 0.05)
  expect_setequal(summ$significance, c("significant", "non_significant"))
  expect_equal(summ$n_pathways[summ$significance == "significant"], 1L)
  expect_equal(summ$n_pathways[summ$significance == "non_significant"], 2L)
})

test_that("median and sum are honoured, and bad arguments are rejected", {
  tab <- data.frame(pathway = c("P1", "P2", "P3"), effect = c(1, 2, 30),
                    stringsAsFactors = FALSE)
  map <- c(P1 = "a", P2 = "a", P3 = "a")
  expect_equal(pt_summarize_by_category(tab, map, stat = "median")$median_effect, 2)
  expect_equal(pt_summarize_by_category(tab, map, stat = "sum")$sum_effect, 33)
  expect_error(pt_summarize_by_category(tab, map, stat = "mode"), "mean\\|median\\|sum")
  expect_error(pt_summarize_by_category(tab, map, key_col = "nope"), "key_col")
  expect_error(pt_summarize_by_category(tab, map, value_col = "nope"), "value_col")
  expect_error(pt_summarize_by_category(tab, map, significance_col = "nope"),
               "significance_col")
})

test_that("a category map is read from disk, first duplicate winning", {
  path <- tempfile(fileext = ".tsv")
  writeLines(c("pathway\tcategory", "P1\talpha", "P2\tbeta", "P1\tGAMMA"), path)
  map <- pt_load_category_map(path)
  expect_equal(unname(map[["P1"]]), "alpha")
  expect_length(map, 2L)
  expect_error(pt_load_category_map(path, key_col = "term"), "column 'term'")
})

test_that("a differential result can be summarised directly", {
  sm <- synthetic_score_matrix()
  g <- pt_group(sm, mode = "target", target_col = "grp")
  d <- pt_diff(sm, g, method = "welch")
  map <- stats::setNames(rep(c("a", "b"), length.out = 40L), rownames(sm$data))
  summ <- pt_summarize_categories(d, map, significance_col = "fdr")
  # grouped by comparison automatically, since the column is present
  expect_true("comparison" %in% names(summ))
  expect_setequal(unique(summ$comparison), pt_comparisons(d))
})
