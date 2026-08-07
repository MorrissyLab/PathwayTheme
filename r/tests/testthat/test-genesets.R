test_that("gmt parsing keeps terms, drops the description, ignores blanks", {
  path <- tempfile(fileext = ".gmt")
  writeLines(c("T1\tsome description\tA\tB\tC",
               "T2\tdesc\tB\tD\t",          # trailing empty field
               "TOO_SHORT\tdesc"), path)    # no members at all
  gs <- pt_parse_gmt(path)
  expect_named(gs, c("T1", "T2"))
  expect_identical(gs$T1, c("A", "B", "C"))
  expect_identical(gs$T2, c("B", "D"))
  expect_setequal(pt_gene_universe(gs), c("A", "B", "C", "D"))
})

test_that("an unparseable gmt is an error rather than an empty result", {
  path <- tempfile(fileext = ".gmt")
  writeLines(c("nothing here", "still nothing"), path)
  expect_error(pt_parse_gmt(path), "No gene sets parsed")
})

test_that("case harmonisation only fires when it improves overlap", {
  gs <- list(T1 = c("BRCA1", "TP53"))

  # lower-case features: upper-casing lifts overlap from 0 to 2
  h <- pt_harmonize_case(c("brca1", "tp53"), gs)
  expect_identical(h$transform("brca1"), "BRCA1")
  expect_identical(h$gene_sets$T1, c("BRCA1", "TP53"))

  # already matching: leave both sides alone
  h <- pt_harmonize_case(c("BRCA1", "TP53"), gs)
  expect_identical(h$transform("Brca1"), "Brca1")
  expect_identical(h$gene_sets, gs)
})
