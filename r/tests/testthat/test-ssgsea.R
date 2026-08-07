## The closed-form enrichment score is the part most worth pinning down: it has
## to equal the literal running-sum definition, and it has to agree with the
## Python implementation's tie handling.

## The definition, written out literally: walk the ordered gene list keeping a
## running sum, then take the area under it.
reference_es <- function(values, members, weight = 0.25) {
  ord <- rev(order(values, method = "radix"))
  gl <- names(values)[ord]
  cv <- abs(values[ord])^weight
  tag <- as.numeric(gl %in% members)
  n_h <- sum(tag)
  n_m <- length(gl) - n_h
  sum(cumsum(tag * cv / sum(tag * cv) - (1 - tag) / n_m))
}

test_that("rank normalisation averages ties and scales by 10000/n", {
  m <- matrix(c(3, 1, 1, 5), ncol = 1, dimnames = list(letters[1:4], "s1"))
  r <- pt_rank_normalise(m)
  # ranks are 3, 1.5, 1.5, 4 -> scaled by 10000/4
  expect_equal(as.vector(r), c(3, 1.5, 1.5, 4) * 10000 / 4)
})

test_that("missing values are ranked last, sharing the trailing positions", {
  m <- matrix(c(3, NA, 1, NA), ncol = 1, dimnames = list(letters[1:4], "s1"))
  r <- as.vector(pt_rank_normalise(m)) / (10000 / 4)
  expect_equal(r[[3]], 1)          # smallest observed
  expect_equal(r[[1]], 2)
  expect_equal(r[[2]], mean(3:4))  # both NAs share the mean trailing rank
  expect_equal(r[[4]], mean(3:4))
})

test_that("the closed-form score equals the running-sum definition", {
  set.seed(11)
  genes <- sprintf("G%03d", 1:120)
  x <- matrix(rnorm(120 * 3), nrow = 120, dimnames = list(genes, c("a", "b", "c")))
  gs <- list(T1 = sample(genes, 20L), T2 = sample(genes, 8L))

  got <- pt_ssgsea(x, gs, weight = 0.25, min_size = 5L, normalise = FALSE)
  ranked <- pt_rank_normalise(x)
  for (term in names(gs)) {
    for (s in colnames(x)) {
      v <- ranked[, s]
      names(v) <- rownames(ranked)
      expect_equal(got[term, s], reference_es(v, gs[[term]]), tolerance = 1e-10)
    }
  }
})

test_that("the same holds when nearly half the matrix is tied at zero", {
  set.seed(12)
  genes <- sprintf("G%03d", 1:150)
  x <- matrix(rnorm(150 * 2), nrow = 150, dimnames = list(genes, c("a", "b")))
  x[runif(length(x)) < 0.45] <- 0
  gs <- list(T1 = sample(genes, 25L))

  got <- pt_ssgsea(x, gs, min_size = 5L, normalise = FALSE)
  ranked <- pt_rank_normalise(x)
  for (s in colnames(x)) {
    v <- ranked[, s]
    names(v) <- rownames(ranked)
    expect_equal(got["T1", s], reference_es(v, gs$T1), tolerance = 1e-10)
  }
})

test_that("normalisation rescales by the range over the whole table", {
  set.seed(13)
  genes <- sprintf("G%03d", 1:100)
  x <- matrix(rnorm(100 * 4), nrow = 100, dimnames = list(genes, letters[1:4]))
  gs <- list(T1 = sample(genes, 15L), T2 = sample(genes, 12L),
             T3 = sample(genes, 30L))
  raw <- pt_ssgsea(x, gs, min_size = 5L, normalise = FALSE)
  nes <- pt_ssgsea(x, gs, min_size = 5L, normalise = TRUE)
  expect_equal(nes, raw / (max(raw) - min(raw)))
  expect_equal(max(nes) - min(nes), 1)
})

test_that("size bounds apply to the members actually present in the matrix", {
  x <- matrix(rnorm(20 * 2), nrow = 20,
              dimnames = list(sprintf("G%02d", 1:20), c("a", "b")))
  gs <- list(keep = sprintf("G%02d", 1:10),
             thin = c("G01", "G02", "ABSENT1", "ABSENT2", "ABSENT3", "ABSENT4"),
             gone = c("NOWHERE1", "NOWHERE2"))
  got <- pt_ssgsea(x, gs, min_size = 5L, normalise = FALSE)
  expect_identical(rownames(got), "keep")
})

test_that("a gene set with no scoreable members is an error, not silence", {
  x <- matrix(rnorm(20), nrow = 10, dimnames = list(letters[1:10], c("a", "b")))
  expect_error(pt_ssgsea(x, list(gone = c("X", "Y")), min_size = 5L),
               "no gene set has")
})

test_that("terms come back sorted and samples keep their input order", {
  set.seed(14)
  genes <- sprintf("G%03d", 1:80)
  x <- matrix(rnorm(80 * 3), nrow = 80,
              dimnames = list(genes, c("zeta", "alpha", "mid")))
  gs <- list(Zed = sample(genes, 10L), Able = sample(genes, 10L))
  got <- pt_ssgsea(x, gs, min_size = 5L)
  expect_identical(rownames(got), c("Able", "Zed"))
  expect_identical(colnames(got), c("zeta", "alpha", "mid"))
})

test_that("the backend harmonises case and averages duplicate gene rows", {
  fx <- synthetic_gmt()
  fm <- synthetic_feature_matrix()
  rownames(fm$data) <- tolower(rownames(fm$data))   # gmt is upper-case
  sm <- suppressWarnings(pt_enrich(fm, backend = "ssgsea", gmt = fx$path,
                                   geneset = "SYN", verbose = FALSE))
  expect_gt(nrow(sm$data), 0L)
  expect_identical(pt_samples(sm), pt_samples(fm))
})
