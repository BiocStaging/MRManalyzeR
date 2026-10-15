# filterBlanks(): mask a sample value that is not far enough above the blanks.
#
# Two blanks and three samples. With blank_filter = 3 the threshold is three
# times the blank mean:
#   A  blanks 10, 10 -> mean 10 -> threshold 30; samples 20, 100, 40
#      so only the 20 is masked.
#   B  blanks 0,  0  -> threshold 0; nothing is below it, nothing masked.
#   C  blanks NA, NA -> no usable blank, so the compound is left alone.

.de_blank <- function() {
  ids <- c("B1", "B2", "S1", "S2", "S3")
  struct::DatasetExperiment(
    data = data.frame(
      A = c(10, 10, 20, 100, 40),
      B = c(0, 0, 5, 5, 5),
      C = c(NA, NA, 7, 8, 9),
      row.names = ids),
    sample_meta = data.frame(
      Sample_type = c("Blank", "Blank", "Sample", "Sample", "Sample"),
      row.names   = ids),
    variable_meta = data.frame(Compound = c("A", "B", "C"),
                               row.names = c("A", "B", "C")))
}

test_that("values below n x the blank mean are masked, and counted", {
  out <- filterBlanks(.de_blank(), blank_filter = 3)
  d   <- as.data.frame(out$data)
  vm  <- as.data.frame(out$variable_meta)

  expect_true(is.na(d["S1", "A"]))          # 20 < 30
  expect_equal(d[c("S2", "S3"), "A"], c(100, 40))
  expect_equal(vm["A", "blank_threshold"], 30)

  # The count is of every cell masked, the blank injections included: both
  # blanks sit at 10 and the threshold they set is 30, so they go too.
  expect_equal(vm["A", "n_masked_blank"], 3)
  expect_equal(sum(is.na(d$A)), 3)
  expect_true(all(is.na(d[c("B1", "B2"), "A"])))
})

test_that("a compound with no usable blank is left alone", {
  out <- filterBlanks(.de_blank(), blank_filter = 3)
  d   <- as.data.frame(out$data)
  vm  <- as.data.frame(out$variable_meta)

  expect_equal(d[c("S1", "S2", "S3"), "C"], c(7, 8, 9))
  expect_true(is.na(vm["C", "blank_threshold"]))
  expect_equal(vm["C", "n_masked_blank"], 0)
  expect_equal(vm["B", "n_masked_blank"], 0)   # threshold of zero masks nothing
})

test_that("median summarises the blanks when asked", {
  # Three blanks, one of them high: mean 30, median 10. The threshold is
  # three times each, so an outlying blank moves the mean and not the median.
  ids <- c("B1", "B2", "B3", "S1", "S2")
  de  <- struct::DatasetExperiment(
    data = data.frame(A = c(10, 10, 70, 50, 100), row.names = ids),
    sample_meta = data.frame(
      Sample_type = c(rep("Blank", 3), "Sample", "Sample"),
      row.names   = ids),
    variable_meta = data.frame(Compound = "A", row.names = "A"))

  mean_thr <- as.data.frame(
    filterBlanks(de, 3)$variable_meta)["A", "blank_threshold"]
  median_thr <- as.data.frame(
    filterBlanks(de, 3, summary = "median")$variable_meta)["A", "blank_threshold"]
  expect_equal(mean_thr, 90)       # (10 + 10 + 70)/3 = 30, x3
  expect_equal(median_thr, 30)     # median 10, x3

  # and that changes what survives: 50 is above 30 but below 90.
  expect_true(is.na(as.data.frame(filterBlanks(de, 3)$data)["S1", "A"]))
  expect_equal(as.data.frame(
    filterBlanks(de, 3, summary = "median")$data)["S1", "A"], 50)
})

test_that("too few blanks is a warning, an error or a skip, as asked", {
  de <- .de_blank()
  sm <- as.data.frame(de$sample_meta)
  sm$Sample_type <- c("Blank", "Sample", "Sample", "Sample", "Sample")
  de$sample_meta <- sm

  expect_warning(filterBlanks(de, 3), "1 blank")
  expect_error(filterBlanks(de, 3, missing_blanks = "error"), "1 blank")
  expect_message(out <- filterBlanks(de, 3, missing_blanks = "skip"),
                 "Skipping")
  expect_equal(as.data.frame(out$data), as.data.frame(de$data))
})

test_that("the blank column has to exist", {
  expect_error(filterBlanks(.de_blank(), 3, blank_head = "Nope"),
               "not a sample_meta column")
  expect_error(filterBlanks(matrix(1:4, 2), 3), "DatasetExperiment")
})
