# addCVMetrics(): per-compound %RSD for the QCs and for the study samples.
#
# The numbers are constructed so each CV is obvious:
#   A  QCs 100, 100, 100      -> sd 0,    CV 0
#      samples 80, 100, 120   -> sd 20,   CV 20
#   B  QCs 90, 100, 110       -> sd 10,   CV 10
#      samples 50, 100, 150   -> sd 50,   CV 50
# so B's sample-to-QC ratio is 5 and A's is Inf (a QC CV of zero).

.de_cv <- function() {
  ids <- paste0("S", 1:6)
  struct::DatasetExperiment(
    data = data.frame(A = c(100, 100, 100, 80, 100, 120),
                      B = c(90, 100, 110, 50, 100, 150),
                      row.names = ids),
    sample_meta = data.frame(
      Sample_type = rep(c("QC", "Sample"), each = 3),
      row.names   = ids),
    variable_meta = data.frame(Compound = c("A", "B"),
                               row.names = c("A", "B")))
}

test_that("CVs are the per-compound %RSD of each sample type", {
  vm <- as.data.frame(addCVMetrics(.de_cv())$variable_meta)

  expect_equal(vm["A", "CV_QC"], 0)
  expect_equal(vm["A", "CV_sample"], 20)
  expect_equal(vm["B", "CV_QC"], 10)
  expect_equal(vm["B", "CV_sample"], 50)
  expect_equal(vm["B", "CV_sample_vs_QC"], 5)
})

test_that("one QC injection gives NA rather than a CV of zero", {
  de <- .de_cv()
  sm <- as.data.frame(de$sample_meta)
  sm$Sample_type <- c("QC", "Sample", "Sample", "Sample", "Sample", "Sample")
  de$sample_meta <- sm

  vm <- as.data.frame(addCVMetrics(de)$variable_meta)
  expect_true(all(is.na(vm$CV_QC)))
  expect_false(any(is.na(vm$CV_sample)))
})

test_that("a second run says it is replacing the first run's CVs", {
  once <- addCVMetrics(.de_cv())
  expect_warning(addCVMetrics(once), "overwriting")
})

test_that("the sample-type column has to exist, and de has to be a dataset", {
  # A missing column is not an error: the dataset comes back untouched, so a
  # run without QCs still produces a report.
  de <- .de_cv()
  expect_silent(out <- addCVMetrics(de, sample_type_head = "Nope"))
  expect_false("CV_QC" %in% colnames(as.data.frame(out$variable_meta)))

  expect_error(addCVMetrics(data.frame(x = 1)), "DatasetExperiment")
})
