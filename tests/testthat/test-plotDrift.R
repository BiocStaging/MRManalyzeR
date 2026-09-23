test_that("plotDrift refuses a non-numeric injection order", {
  ids <- paste0("S", 1:3)
  de <- struct::DatasetExperiment(
    data = data.frame(A = c(1, 2, 3), row.names = ids),
    sample_meta = data.frame(Sample_type = "QC",
                             Injection_order = c("_1", "_2", "_3"),
                             row.names = ids),
    variable_meta = data.frame(Compound = "A", row.names = "A"))
  expect_error(plotDrift(de, "A"), "must be numeric")
})
