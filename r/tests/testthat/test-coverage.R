test_that("every requested gene set gets a recorded fate", {
  features <- c("A", "B", "C", "D", "E")
  gs <- list(fine = c("A", "B", "C"),
             small = c("A", "ZZ"),
             absent = c("X", "Y"),
             huge = c("A", "B", "C", "D", "E"))
  cov <- pt_geneset_coverage(features, gs, min_size = 3L, max_size = 4L)

  expect_setequal(cov$pathway, names(gs))
  status <- stats::setNames(cov$status, cov$pathway)
  expect_equal(status[["fine"]], PT_SCORED)
  expect_equal(status[["small"]], PT_TOO_SMALL)
  expect_equal(status[["absent"]], PT_NO_OVERLAP)
  expect_equal(status[["huge"]], PT_TOO_LARGE)

  n_present <- stats::setNames(cov$n_present, cov$pathway)
  expect_equal(n_present[["small"]], 1L)
  expect_equal(cov$fraction_present[cov$pathway == "small"], 0.5)
})

test_that("coverage matches on upper-cased symbols, as the backend does", {
  cov <- pt_geneset_coverage(c("brca1", "tp53"), list(T1 = c("BRCA1", "TP53")),
                             min_size = 1L, max_size = 10L)
  expect_equal(cov$n_present, 2L)
  expect_equal(cov$status, PT_SCORED)
})

test_that("dropped sets raise a warning naming them", {
  cov <- pt_geneset_coverage(c("A"), list(gone = c("X", "Y")), 2L, 10L)
  expect_warning(pt_warn_on_dropped(cov, geneset = "SYN"), "will not be scored")
  expect_warning(pt_warn_on_dropped(cov), "gone")
  expect_silent(pt_warn_on_dropped(cov, verbose = FALSE))
})

test_that("the default missing policy changes nothing but says so", {
  m <- matrix(c(1, NA, 3, 4), nrow = 2, dimnames = list(c("a", "b"), c("s1", "s2")))
  expect_warning(out <- pt_apply_missing_policy(m, "as_is"), "as_is")
  expect_identical(out, m)
})

test_that("zero and low policies fill explicitly", {
  m <- matrix(c(1, NA, 3, 4), nrow = 2, dimnames = list(c("a", "b"), c("s1", "s2")))
  z <- suppressMessages(pt_apply_missing_policy(m, "zero"))
  expect_equal(z[["b", "s1"]], 0)

  low <- suppressMessages(pt_apply_missing_policy(m, "low"))
  # one unit below sample s1's own observed minimum (1)
  expect_equal(low[["b", "s1"]], 0)

  m2 <- matrix(c(5, NA, 3, 4), nrow = 2, dimnames = list(c("a", "b"), c("s1", "s2")))
  low2 <- suppressMessages(pt_apply_missing_policy(m2, "low"))
  expect_equal(low2[["b", "s1"]], 4)
})

test_that("under-observed features are dropped before filling", {
  m <- matrix(c(1, NA, NA, 2, NA, NA), nrow = 3,
              dimnames = list(c("a", "b", "c"), c("s1", "s2")))
  out <- suppressMessages(suppressWarnings(
    pt_apply_missing_policy(m, "as_is", min_observed_fraction = 0.9)))
  expect_identical(rownames(out), "a")
})

test_that("an unknown policy is rejected", {
  expect_error(pt_apply_missing_policy(matrix(1), "guess"), "must be one of")
})
