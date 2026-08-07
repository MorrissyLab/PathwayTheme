test_that("matrix contracts coerce names to character and keep order", {
  x <- matrix(1:6, nrow = 3, dimnames = list(1:3, 10:11))
  fm <- pt_feature_matrix(x, data.frame(a = 1:2, row.names = c("10", "11")))
  expect_identical(pt_features(fm), c("1", "2", "3"))
  expect_identical(pt_samples(fm), c("10", "11"))
  expect_type(fm$data, "double")
})

test_that("metadata alignment reorders and reports what is missing", {
  meta <- pt_sample_metadata(data.frame(grp = c("a", "b", "c"),
                                        row.names = c("S1", "S2", "S3")))
  aligned <- pt_align_to(meta, c("S3", "S1"))
  expect_identical(rownames(aligned$table), c("S3", "S1"))
  expect_identical(unname(pt_meta_get(aligned, "grp")), c("c", "a"))
  expect_error(pt_align_to(meta, c("S1", "S9")), "missing from metadata")
  expect_error(pt_meta_get(meta, "nope"), "not found")
  expect_false(pt_meta_has(meta, "nope"))
})

test_that("scopes fall back to a single bucket and keep first-appearance order", {
  labels <- stats::setNames(c("a", "b", "a"), c("S1", "S2", "S3"))
  expect_identical(names(pt_scopes(pt_grouping(labels))), "__all__")

  scope <- stats::setNames(c("z", "y", "z"), c("S1", "S2", "S3"))
  g <- pt_grouping(labels, scope)
  sc <- pt_scopes(g)
  expect_identical(names(sc), c("z", "y"))     # not alphabetical
  expect_identical(sc$z, c("S1", "S3"))
})

test_that("score matrices subset samples with their metadata", {
  sm <- synthetic_score_matrix()
  sub <- pt_subset_samples(sm, c("S03", "S00"))
  expect_identical(pt_samples(sub), c("S03", "S00"))
  expect_identical(rownames(sub$metadata$table), c("S03", "S00"))
  expect_equal(sub$data[, "S00"], sm$data[, "S00"])
  expect_identical(sub$backend, sm$backend)
})

test_that("differential results expose comparisons and significant rows", {
  tab <- data.frame(comparison = c("a_vs_rest", "a_vs_rest", "b_vs_rest"),
                    pathway = c("P1", "P2", "P1"),
                    fdr = c(0.01, 0.4, NA_real_), stringsAsFactors = FALSE)
  res <- pt_diff_result(tab, method = "welch")
  expect_identical(pt_comparisons(res), c("a_vs_rest", "b_vs_rest"))
  sig <- pt_significant(res, 0.05)
  expect_equal(nrow(sig), 1L)          # NA fdr is not significant
  expect_identical(sig$pathway, "P1")
})

test_that("descending order breaks ties the way NumPy and pandas do", {
  x <- c(a = 1, b = 1, c = 2)
  # a stable ascending order reversed: c, then b, then a
  expect_identical(pathwaytheme:::.pt_order_desc(x), c(3L, 2L, 1L))
})
