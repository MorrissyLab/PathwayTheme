test_that("PCA reproduces the score matrix it decomposed", {
  sm <- synthetic_score_matrix()
  # keep every component the data supports, so the reconstruction is exact
  res <- pt_pca(sm, k_max = ncol(sm$data) - 1L)[[1L]]

  x <- sm$data
  mu <- rowMeans(x)
  sds <- sqrt(rowMeans((x - mu)^2))
  xn <- (x - mu) / sds
  # scores %*% loadings rebuilds the centred, scaled matrix within the kept rank
  rebuilt <- t(res$scores %*% res$loadings)
  expect_equal(rebuilt, xn[rownames(rebuilt), ], tolerance = 1e-8)
})

test_that("variance explained is a decreasing fraction summing to at most one", {
  res <- pt_pca(synthetic_score_matrix())[[1L]]
  v <- res$variance_explained
  expect_true(all(diff(v) <= 1e-12))
  expect_true(all(v >= 0))
  expect_lte(sum(v), 1 + 1e-12)
})

test_that("component signs are fixed by the largest-magnitude loading", {
  res <- pt_pca(synthetic_score_matrix())[[1L]]
  for (k in seq_len(pt_k(res))) {
    lv <- res$loadings[k, ]
    expect_gt(lv[[which.max(abs(lv))]], 0)
  }
})

test_that("the same input gives the same signs every time", {
  sm <- synthetic_score_matrix()
  a <- pt_pca(sm)[[1L]]
  b <- pt_pca(sm)[[1L]]
  expect_equal(a$loadings, b$loadings)
  expect_equal(a$scores, b$scores)
})

test_that("K is capped by both k_max and the observation count", {
  sm <- synthetic_score_matrix()                 # 12 samples
  expect_equal(pt_k(pt_pca(sm)[[1L]]), 10L)      # k_max default
  expect_equal(pt_k(pt_pca(sm, k_max = 4L)[[1L]]), 4L)
  small <- pt_subset_samples(sm, pt_samples(sm)[1:5])
  expect_equal(pt_k(pt_pca(small)[[1L]]), 4L)    # n_obs - 1
})

test_that("scopes too small to decompose are skipped, not faked", {
  sm <- synthetic_score_matrix()
  expect_null(pt_run_pca(sm$data[, 1:2, drop = FALSE],
                         stats::setNames(c("a", "b"), pt_samples(sm)[1:2]),
                         pt_pca_config()))
  expect_null(pt_run_pca(sm$data[1:3, , drop = FALSE],
                         stats::setNames(rep("a", 12L), pt_samples(sm)),
                         pt_pca_config()))
})

test_that("constant features are dropped before the decomposition", {
  sm <- synthetic_score_matrix()
  sm$data["TERM_00", ] <- 5                       # zero variance
  res <- pt_pca(sm)[[1L]]
  expect_false("TERM_00" %in% colnames(res$loadings))
  expect_equal(ncol(res$loadings), nrow(sm$data) - 1L)
})

test_that("missing scores are filled with zero, loudly", {
  sm <- synthetic_score_matrix()
  sm$data[1L, 1L] <- NA_real_
  expect_warning(pt_pca(sm), "missing score")
})

test_that("a scope column splits the run and keeps the scope labels", {
  sm <- synthetic_score_matrix()
  g <- pt_group(sm, mode = "target", target_col = "grp", scope_col = "sample_id")
  res <- pt_pca(sm, g)
  expect_length(res, 2L)
  expect_setequal(vapply(res, function(r) r$scope, ""), c("A", "B"))
  expect_equal(nrow(res[[1L]]$scores), 6L)
})

test_that("signatures combine the top PCs and rank by absolute score", {
  res <- pt_pca(synthetic_score_matrix(), top_pcs_per_cluster = 2L,
                top_pathways_per_cluster = 5L)[[1L]]
  sig <- res$signatures
  expect_equal(nrow(sig), nrow(res$scores) * 5L)
  expect_setequal(unique(sig$rank), 1:5)

  first <- sig[sig$cluster == sig$cluster[[1L]], ]
  expect_true(all(diff(abs(first$signature_score)) <= 1e-12))
  expect_identical(first$direction, ifelse(first$signature_score > 0, "up", "down"))

  # the score really is the sum over the top PCs of score * loading
  pcs <- strsplit(first$top_pcs[[1L]], ",", fixed = TRUE)[[1L]]
  expect_length(pcs, 2L)
  obs <- first$cluster[[1L]]
  manual <- res$scores[obs, pcs] %*% res$loadings[pcs, first$pathway[[1L]]]
  expect_equal(as.numeric(manual), first$signature_score[[1L]])
})

test_that("the per-PC table is ordered by absolute loading", {
  res <- pt_pca(synthetic_score_matrix())[[1L]]
  tab <- pt_top_pathways_per_pc(res$loadings, res$variance_explained, 5L)
  expect_equal(nrow(tab), pt_k(res) * 5L)
  pc1 <- tab[tab$PC == "PC1", ]
  expect_true(all(diff(abs(pc1$loading)) <= 1e-12))
  expect_equal(pc1$variance_explained[[1L]], res$variance_explained[[1L]])
})

test_that("display labels, sublabels and sizes ride along without changing maths", {
  sm <- synthetic_score_matrix()
  g <- pt_group(sm, mode = "target", target_col = "grp")
  plain <- pt_pca(sm, g)[[1L]]
  annotated <- pt_pca(sm, g, display_col = "cluster_id",
                      secondary_col = "cell_type", size_col = "n_cells")[[1L]]
  expect_equal(unname(annotated$scores), unname(plain$scores))
  expect_identical(rownames(annotated$scores), sm$metadata$table$cluster_id)
  expect_setequal(unique(annotated$sublabels), unique(sm$metadata$table$cell_type))
  expect_equal(annotated$sizes, as.double(sm$metadata$table$n_cells))
})
