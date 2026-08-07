test_that("over-representation matches a hand-computed hypergeometric test", {
  background <- sprintf("G%02d", 1:100)
  gs <- list(hit = sprintf("G%02d", 1:20), miss = sprintf("G%02d", 81:100))
  query <- sprintf("G%02d", 1:10)

  res <- pt_ora(query, gs, background)
  hit <- res[res$Term == "hit", ]
  # 10 of the query's 10 genes are in a 20-gene set out of 100
  expect_equal(hit$P.value, stats::phyper(9, 20, 80, 10, lower.tail = FALSE))
  expect_equal(hit$Overlap, "10/20")
  # Haldane-Anscombe corrected odds ratio
  expect_equal(hit$Odds.Ratio, ((10 + 0.5) * (100 - 20 - 10 + 10 + 0.5)) /
                 ((20 - 10 + 0.5) * (10 - 10 + 0.5)))
  expect_equal(hit$Combined.Score, -log(hit$P.value) * hit$Odds.Ratio)
  # a set with no overlap is omitted entirely
  expect_false("miss" %in% res$Term)
})

test_that("the query and the sets are both restricted to the background", {
  res <- pt_ora(c("A", "B", "OFFPIECE"), list(T1 = c("A", "B", "OFFPIECE")),
                background = c("A", "B", "C", "D"))
  expect_equal(res$Overlap, "2/2")
})

test_that("gene-list selection honours top_n and threshold", {
  col <- stats::setNames(c(5, 1, 4, NA, 3), letters[1:5])
  cfg <- pt_enrichment_config(gene_selection = "top_n", top_n = 2L)
  expect_identical(pt_select_genes(col, cfg), c("a", "c"))

  cfg <- pt_enrichment_config(gene_selection = "threshold", threshold = 3.5)
  expect_setequal(pt_select_genes(col, cfg), c("a", "c"))

  cfg <- pt_enrichment_config(gene_selection = "threshold")
  expect_error(pt_select_genes(col, cfg), "requires enrichment\\$threshold")
})

test_that("goslim scores the fraction or the count of a sample's genes", {
  fx <- synthetic_gmt()
  fm <- synthetic_feature_matrix()
  frac <- pt_enrich(fm, backend = "goslim", gmt = fx$path, geneset = "SYN",
                    top_n = 50L, verbose = FALSE)
  cnt <- pt_enrich(fm, backend = "goslim", gmt = fx$path, geneset = "SYN",
                   top_n = 50L, goslim_score = "count", verbose = FALSE)
  common <- intersect(rownames(frac$data), rownames(cnt$data))
  expect_equal(frac$data[common, ], cnt$data[common, ] / 50)
  expect_true(all(frac$data >= 0 & frac$data <= 1))
})

test_that("an unknown backend names the ones that exist", {
  expect_error(pt_get_backend("nosuch"), "Unknown enrichment backend")
  expect_error(pt_get_backend("nosuch"), "ssgsea")
})

test_that("backends can be registered and are looked up case-insensitively", {
  pt_register_backend("Constant", function(fm, config, verbose) {
    matrix(1, nrow = 2L, ncol = ncol(fm$data),
           dimnames = list(c("T1", "T2"), pt_samples(fm)))
  })
  fm <- synthetic_feature_matrix()
  sm <- pt_enrich(fm, backend = "constant", geneset = "SYN", verbose = FALSE)
  expect_equal(dim(sm$data), c(2L, 12L))
  expect_true(all(sm$data == 1))
})

test_that("the cache is reused, and extended one sample at a time", {
  fx <- synthetic_gmt()
  fm <- synthetic_feature_matrix()
  cache_dir <- tempfile()

  first <- suppressWarnings(pt_enrich(pt_filter_samples(fm, "grp", "hi"),
                                      backend = "ssgsea", gmt = fx$path,
                                      geneset = "SYN", cache_dir = cache_dir,
                                      verbose = FALSE))
  cache_file <- file.path(cache_dir, "SYN_ssgsea_scores.tsv")
  expect_true(file.exists(cache_file))
  expect_equal(ncol(first$data), 6L)

  # asking for all 12 must score the 6 that are missing, not drop them
  full <- suppressWarnings(suppressMessages(
    pt_enrich(fm, backend = "ssgsea", gmt = fx$path, geneset = "SYN",
              cache_dir = cache_dir, verbose = FALSE)))
  expect_equal(ncol(full$data), 12L)
  # the cached six are returned unchanged
  expect_equal(full$data[rownames(first$data), colnames(first$data)],
               first$data, tolerance = 1e-12)
})

test_that("the missing-value policy is part of the cache key", {
  cfg <- pt_enrichment_config(geneset = "SYN", cache_dir = "cache")
  expect_match(pathwaytheme:::.pt_cache_path(cfg), "SYN_ssgsea_scores.tsv$")
  cfg$missing <- "low"
  expect_match(pathwaytheme:::.pt_cache_path(cfg), "miss-low")
})

test_that("coverage travels with the score matrix", {
  fx <- synthetic_gmt()
  fm <- synthetic_feature_matrix()
  sm <- suppressWarnings(pt_enrich(fm, backend = "ssgsea", gmt = fx$path,
                                   geneset = "SYN", verbose = FALSE))
  expect_s3_class(sm$coverage, "data.frame")
  expect_true(all(c("pathway", "n_present", "status") %in% names(sm$coverage)))
})
