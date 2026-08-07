test_that("Welch's test matches stats::t.test row by row", {
  set.seed(21)
  case <- matrix(rnorm(5 * 6, mean = 1), nrow = 5)
  ref <- matrix(rnorm(5 * 7), nrow = 5)
  got <- pt_welch(case, ref)
  for (i in 1:5) {
    tt <- stats::t.test(case[i, ], ref[i, ], var.equal = FALSE)
    expect_equal(got$statistic[[i]], unname(tt$statistic))
    expect_equal(got$p_value[[i]], tt$p.value)
  }
  expect_equal(got$effect, rowMeans(case) - rowMeans(ref))
})

test_that("Mann-Whitney matches stats::wilcox.test", {
  set.seed(22)
  case <- matrix(rnorm(4 * 9, mean = 0.8), nrow = 4)
  ref <- matrix(rnorm(4 * 9), nrow = 4)
  got <- pt_mannwhitney(case, ref)
  for (i in 1:4) {
    w <- stats::wilcox.test(case[i, ], ref[i, ], exact = FALSE, correct = TRUE)
    expect_equal(got$statistic[[i]], unname(w$statistic))
    expect_equal(got$p_value[[i]], w$p.value)
  }
})

test_that("trigammaInverse really inverts trigamma", {
  for (y in c(0.3, 1, 4, 25)) {
    expect_equal(pt_trigamma_inverse(trigamma(y)), y, tolerance = 1e-6)
  }
  expect_equal(pt_trigamma_inverse(1e8), 1 / sqrt(1e8))
  expect_equal(pt_trigamma_inverse(1e-8), 1e8)
})

test_that("fitFDist recovers the prior it was given", {
  set.seed(23)
  d0_true <- 8
  s0_true <- 2.5
  df <- 6
  # variances drawn from the scaled-F prior the model assumes
  s2 <- s0_true * stats::rchisq(4000, d0_true) / d0_true *
    stats::rchisq(4000, df) / df
  prior <- pt_fit_fdist(s2, rep(df, length(s2)))
  expect_equal(prior$d0, d0_true, tolerance = 0.25 * d0_true)
  expect_equal(prior$s0_sq, s0_true, tolerance = 0.15 * s0_true)
})

test_that("the moderated t-test agrees with limma", {
  skip_if_not_installed("limma")
  set.seed(24)
  n1 <- 6L
  n2 <- 7L
  y <- matrix(rnorm(60 * (n1 + n2)), nrow = 60)
  y[1:10, seq_len(n1)] <- y[1:10, seq_len(n1)] + 1.5
  case <- y[, seq_len(n1), drop = FALSE]
  ref <- y[, n1 + seq_len(n2), drop = FALSE]

  got <- pt_moderated_t(case, ref)
  design <- cbind(Intercept = 1, Group = c(rep(1, n1), rep(0, n2)))
  fit <- limma::eBayes(limma::lmFit(y, design))
  expect_equal(got$statistic, unname(fit$t[, "Group"]), tolerance = 1e-8)
  expect_equal(got$p_value, unname(fit$p.value[, "Group"]), tolerance = 1e-8)
  expect_equal(got$effect, unname(fit$coefficients[, "Group"]), tolerance = 1e-8)
})

test_that("BH adjustment matches p.adjust and passes NAs through", {
  set.seed(25)
  p <- runif(50)
  expect_equal(pt_benjamini_hochberg(p), stats::p.adjust(p, "BH"))

  p2 <- c(p, NA_real_)
  got <- pt_benjamini_hochberg(p2)
  expect_true(is.na(got[[51L]]))
  # the NA is excluded from the ranking, so the rest are unchanged
  expect_equal(got[1:50], stats::p.adjust(p, "BH"))
  expect_true(all(is.na(pt_benjamini_hochberg(c(NA_real_, NA_real_)))))
})

test_that("an unknown method names the ones that exist", {
  expect_error(pt_run_test("ttest", matrix(1), matrix(1)), "Unknown diff method")
  expect_error(pt_run_test("ttest", matrix(1), matrix(1)), "moderated_t")
})

test_that("one-vs-rest runs a comparison per group", {
  sm <- synthetic_score_matrix()
  g <- pt_group(sm, mode = "target", target_col = "grp")
  d <- pt_diff(sm, g, method = "welch")
  expect_setequal(pt_comparisons(d), c("hi_vs_rest", "lo_vs_rest"))
  expect_equal(nrow(d$table), 2L * nrow(sm$data))
  expect_equal(unique(d$table$n_case), 6L)
  expect_equal(unique(d$table$n_reference), 6L)
})

test_that("a reference group compares everything else against it", {
  sm <- synthetic_score_matrix()
  d <- pt_diff(sm, group_col = "cell_type", reference = "tumor")
  expect_setequal(pt_comparisons(d),
                  c("T-cell_vs_tumor", "macrophage_vs_tumor"))
})

test_that("explicit contrasts run exactly what was asked for", {
  sm <- synthetic_score_matrix()
  d <- pt_diff(sm, group_col = "cell_type",
               contrasts = list(c("tumor", "T-cell")))
  expect_identical(pt_comparisons(d), "tumor_vs_T-cell")
})

test_that("groups smaller than min_group_size are skipped", {
  sm <- synthetic_score_matrix()
  sm$metadata$table$tiny <- c("solo", rep("rest", 11L))
  d <- pt_diff(sm, group_col = "tiny", min_group_size = 2L)
  expect_false("solo_vs_rest" %in% pt_comparisons(d))
})

test_that("rows are sorted by FDR then p, and direction follows the effect", {
  sm <- synthetic_score_matrix()
  g <- pt_group(sm, mode = "target", target_col = "grp")
  d <- pt_diff(sm, g, method = "moderated_t")
  one <- d$table[d$table$comparison == "hi_vs_rest", ]
  expect_true(all(diff(one$fdr) >= -1e-12))
  expect_identical(one$direction, ifelse(one$effect >= 0, "up", "down"))
  expect_equal(one$effect, one$mean_case - one$mean_reference, tolerance = 1e-12)
})

test_that("differential analysis without labels is an error", {
  expect_error(pt_run_diff(synthetic_score_matrix(), pt_diff_config()),
               "needs diff\\$group_col or a grouping")
  expect_error(pt_diff(synthetic_score_matrix(), group_col = "nope"),
               "not in metadata")
})
