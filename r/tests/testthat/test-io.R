write_matrix <- function(x, path, sep = "\t") {
  utils::write.table(x, path, sep = sep, quote = FALSE, col.names = NA)
  path
}

test_that("a tab-separated matrix loads with metadata attached", {
  x <- matrix(1:6, nrow = 3, dimnames = list(c("g1", "g2", "g3"), c("s1", "s2")))
  mp <- write_matrix(x, tempfile(fileext = ".tsv"))
  meta_path <- tempfile(fileext = ".tsv")
  writeLines(c("\tgrp", "s1\thi", "s2\tlo"), meta_path)

  fm <- pt_load_matrix(mp, metadata = meta_path)
  expect_equal(dim(fm$data), c(3L, 2L))
  expect_identical(pt_samples(fm), c("s1", "s2"))
  expect_identical(unname(pt_meta_get(fm$metadata, "grp")), c("hi", "lo"))
})

test_that("without metadata each sample gets a minimal record", {
  x <- matrix(1:4, nrow = 2, dimnames = list(c("g1", "g2"), c("s1", "s2")))
  fm <- pt_load_matrix(write_matrix(x, tempfile(fileext = ".tsv")))
  expect_identical(unname(pt_meta_get(fm$metadata, "sample")), c("s1", "s2"))
})

test_that("a sample missing from the metadata is an error, not a silent drop", {
  x <- matrix(1:4, nrow = 2, dimnames = list(c("g1", "g2"), c("s1", "s2")))
  mp <- write_matrix(x, tempfile(fileext = ".tsv"))
  meta_path <- tempfile(fileext = ".tsv")
  writeLines(c("\tgrp", "s1\thi"), meta_path)
  expect_error(pt_load_matrix(mp, metadata = meta_path),
               "absent from metadata")
})

test_that("csv is detected by extension and non-numeric cells become NA", {
  path <- tempfile(fileext = ".csv")
  writeLines(c(",s1,s2", "g1,1,2", "g2,broken,4"), path)
  fm <- pt_load_matrix(path)
  expect_true(is.na(fm$data[["g2", "s1"]]))
  expect_equal(fm$data[["g2", "s2"]], 4)
})

test_that("unsupported formats say so rather than failing obscurely", {
  expect_error(pt_load_matrix(tempfile(fileext = ".parquet")), "not supported")
  expect_error(pt_load_input(pt_input_config(kind = "sideways")),
               "Unknown input kind")
  expect_error(pt_load_input(pt_input_config(kind = "matrix")),
               "matrix_path is required")
})

test_that("h5ad input reports the missing dependency clearly", {
  skip_if(requireNamespace("hdf5r", quietly = TRUE),
          "hdf5r is installed, so the fallback path is not exercised")
  expect_error(pt_load_h5ad(directory = tempdir()), "hdf5r")
})
