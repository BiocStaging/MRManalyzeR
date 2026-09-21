# imputeMissing(): fills below each compound's lowest measurement.
#
# Eight injections, S8 a blank. With scalar = 0.2 and no batches:
#   A  non-blank values 10, 20, 40, 50 -> ceiling 2; S3, S4, S7 missing; the
#      blank's own NA must stay NA.
#   B  S2 is a zero, which counts as missing; minimum 5 -> ceiling 1.
#   C  non-blank minimum 2 -> ceiling 0.4. The blank's 0.5 is lower, and would
#      give 0.1 if blanks counted.
# Split by Batch (S1-S4 = b1, S5-S8 = b2): A's ceiling is 2 in b1 (min 10) and
# 8 in b2 (min 40).

mk_imp <- function() {
  ids <- paste0("S", 1:8)
  struct::DatasetExperiment(
    data = data.frame(
      A = c(10, 20, NA, NA, 40, 50, NA, NA),
      B = c(5, 0, 5, 5, 5, 5, 5, 5),
      C = c(NA, 2, 3, 4, 100, NA, 200, 0.5),
      row.names = ids),
    sample_meta = data.frame(
      Sample_type = c(rep("Sample", 7), "Blank"),
      Batch       = rep(c("b1", "b2"), each = 4),
      row.names   = ids),
    variable_meta = data.frame(Compound = c("A", "B", "C"),
                               row.names = c("A", "B", "C")))
}

test_that("random fills lie between 0 and the ceiling and are not all equal", {
  d <- as.data.frame(imputeMissing(mk_imp(), 0.2)$data)
  a <- d$A[c(3, 4, 7)]
  expect_true(all(a > 0 & a < 2))
  expect_gt(length(unique(a)), 1)
  expect_true(is.na(d$A[8]))   # the blank is left alone
})

test_that("the same seed gives the same fills, another seed does not", {
  a1 <- as.data.frame(imputeMissing(mk_imp(), 0.2)$data)$A
  a2 <- as.data.frame(imputeMissing(mk_imp(), 0.2)$data)$A
  a3 <- as.data.frame(imputeMissing(mk_imp(), 0.2, seed = 7)$data)$A
  expect_identical(a1, a2)
  expect_false(identical(a1, a3))
})

test_that("constant fills every cell with the ceiling, as before", {
  out <- imputeMissing(mk_imp(), 0.2, method = "constant")
  d  <- as.data.frame(out$data)
  vm <- as.data.frame(out$variable_meta)
  expect_equal(d$A[c(3, 4, 7)], c(2, 2, 2))
  expect_equal(d$B[2], 1)             # a zero counts as missing
  expect_equal(d$C[c(1, 6)], c(0.4, 0.4))   # the blank's 0.5 does not count
  expect_equal(vm["A", "impute_max"], 2)
  expect_equal(vm["A", "n_imputed"], 3)
  expect_equal(vm["A", "frac_imputed"], round(3 / 7, 3))
})

test_that("by gives each batch its own ceiling and keeps the record", {
  out <- imputeMissing(mk_imp(), 0.2, method = "constant", by = "Batch")
  d  <- as.data.frame(out$data)
  vm <- as.data.frame(out$variable_meta)
  expect_equal(d$A[c(3, 4, 7)], c(2, 2, 8))
  expect_equal(unlist(vm["A", c("impute_max_b1", "impute_max_b2")],
                      use.names = FALSE), c(2, 8))
  expect_equal(vm["A", "n_imputed"], 3)
  expect_error(imputeMissing(mk_imp(), 0.2, by = "Nope"),
               "not a sample_meta column")
})
