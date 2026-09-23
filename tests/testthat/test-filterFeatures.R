# Tests for filterFeatures() / filterSamples().
#
# The arithmetic is deliberately native rather than delegated to pmp: the
# calculation is a column mean of !is.na, MRManalyzeR's matrices are
# samples-in-rows where pmp wants features-in-rows, and a runtime dependency
# for that is not worth the transpose. pmp is used here instead as an ORACLE -
# the last test asserts the two agree, so the claim "matches
# pmp::filter_peaks_by_fraction" is checked rather than asserted.

# npOnly is the case the whole design exists for: absent from Control, present
# throughout NP. 'within' must keep it and 'across' must not.
mk_de <- function() {
  struct::DatasetExperiment(
    data = data.frame(
      good   = c(10, 12, 11, 9, 10, 11, 50, 51),
      npOnly = c(NA, NA, NA, 20, 22, 21, 30, NA),
      sparse = c(5, NA, NA, NA, NA, 6, NA, NA),
      qcOnly = c(NA, NA, NA, NA, NA, NA, 99, 98),
      row.names = c(paste0("S", 1:6), "QC1", "QC2")),
    sample_meta = data.frame(
      Sample_type = c(rep("Sample", 6), "QC", "QC"),
      Treatment   = c(rep("Control", 3), rep("NP", 3), "QC", "QC"),
      row.names   = c(paste0("S", 1:6), "QC1", "QC2")),
    variable_meta = data.frame(
      Compound  = c("good", "npOnly", "sparse", "qcOnly"),
      row.names = c("good", "npOnly", "sparse", "qcOnly")))
}

test_that("'within' keeps a feature that is real in only one group", {
  out <- suppressMessages(filterFeatures(mk_de(), min_frac = 0.5,
                                         method = "within",
                                         group_col = "Treatment"))
  expect_true("npOnly" %in% colnames(out$data))   # 3/3 within NP
  expect_false("sparse" %in% colnames(out$data))  # 1/3 in both groups
})

test_that("'across' deletes the same feature - the trap 'within' avoids", {
  # npOnly is detected in 3 of the 6 study samples, so its across-fraction is
  # exactly 0.5. Threshold deliberately above that: at 0.5 the >= rule would
  # keep it, which is correct but tests nothing.
  out <- suppressMessages(filterFeatures(mk_de(), min_frac = 0.6,
                                         method = "across"))
  expect_false("npOnly" %in% colnames(out$data))
  expect_true("good" %in% colnames(out$data))
})

test_that("min_frac is inclusive - a feature exactly on the threshold is kept", {
  # Pins the boundary rule the documentation states, so a future change to
  # > vs >= fails here rather than silently shifting every result.
  keep <- suppressMessages(filterFeatures(mk_de(), min_frac = 0.5,
                                          method = "across"))
  drop <- suppressMessages(filterFeatures(mk_de(), min_frac = 0.5 + 1e-9,
                                          method = "across"))
  expect_true("npOnly" %in% colnames(keep$data))   # 0.5 >= 0.5
  expect_false("npOnly" %in% colnames(drop$data))
})

test_that("QC and blank injections never count towards the fraction", {
  # qcOnly is present in both QCs and no study sample. If QCs counted it would
  # score 2/8; excluded, it scores 0 and must go.
  out <- suppressMessages(filterFeatures(mk_de(), min_frac = 0.25,
                                         method = "across"))
  expect_false("qcOnly" %in% colnames(out$data))
  # ...but the QC ROWS are still in the dataset - correctBatch() needs them
  expect_true(all(c("QC1", "QC2") %in% rownames(out$data)))
})

test_that("method 'QC' falls back with a warning when there are no QCs", {
  de <- mk_de()
  sm <- as.data.frame(de$sample_meta)
  sm$Sample_type <- "Sample"
  de$sample_meta <- sm
  expect_warning(
    out <- suppressMessages(filterFeatures(de, min_frac = 0.5, method = "QC",
                                           group_col = "Treatment")),
    "no QC injections")
  expect_s4_class(out, "DatasetExperiment")
})

test_that("detected_frac is recorded for every retained feature", {
  out <- suppressMessages(filterFeatures(mk_de(), min_frac = 0,
                                         method = "within",
                                         group_col = "Treatment"))
  vm <- as.data.frame(out$variable_meta)
  expect_true("detected_frac" %in% colnames(vm))
  expect_equal(vm["good", "detected_frac"], 1)
  expect_equal(vm["npOnly", "detected_frac"], 1)   # best group, not overall
})

test_that("filterSamples drops failing study rows but only reports failing QCs", {
  de <- struct::DatasetExperiment(
    data = data.frame(
      f1 = c(10, 11, NA, NA), f2 = c(20, 21, NA, NA), f3 = c(30, 31, NA, 33),
      row.names = c("S1", "S2", "S3_bad", "QC_bad")),
    sample_meta = data.frame(
      Sample_type = c("Sample", "Sample", "Sample", "QC"),
      row.names   = c("S1", "S2", "S3_bad", "QC_bad")),
    variable_meta = data.frame(Compound = c("f1", "f2", "f3"),
                               row.names = c("f1", "f2", "f3")))
  expect_warning(
    out <- suppressMessages(filterSamples(de, max_na = 0.5)),
    "were KEPT")
  expect_false("S3_bad" %in% rownames(out$data))
  expect_true("QC_bad" %in% rownames(out$data))
  expect_true("na_frac" %in% colnames(as.data.frame(out$sample_meta)))
})

test_that("our fractions agree with pmp::filter_peaks_by_fraction", {
  skip_if_not_installed("pmp")
  de <- mk_de()
  X  <- as.data.frame(de$data)
  sm <- as.data.frame(de$sample_meta)
  study <- sm$Sample_type == "Sample"

  # pmp wants features in rows, samples in columns
  peak_tbl <- t(as.matrix(X[study, , drop = FALSE]))
  classes  <- as.character(sm$Treatment[study])

  for (mf in c(0.34, 0.5, 0.67, 1)) {
    ours <- suppressMessages(filterFeatures(de, min_frac = mf,
                                            method = "within",
                                            group_col = "Treatment"))
    theirs <- suppressWarnings(suppressMessages(
      pmp::filter_peaks_by_fraction(peak_tbl, min_frac = mf,
                                    classes = classes, method = "within")))
    # rownames() works on a matrix and on a SummarizedExperiment alike, so
    # the test does not need SummarizedExperiment declared.
    kept_pmp <- rownames(theirs)
    expect_setequal(colnames(ours$data), kept_pmp)
  }
})
