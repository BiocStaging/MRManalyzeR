# End-to-end smoke test for mrmRunExample(). Renders both reports and
# checks the four expected output artefacts exist. Slow (~30-60s) and
# pulls in plotly + rmarkdown — skipped on CRAN and skipped unless the
# env var MRMANALYZER_RUN_E2E is set, so the default `devtools::test()`
# stays fast. Set the env var to force a full smoke check:
#   Sys.setenv(MRMANALYZER_RUN_E2E = "true"); devtools::test()

test_that("mrmRunExample produces both reports + xlsx + RDS", {
  skip_on_cran()
  skip_if_not(nzchar(Sys.getenv("MRMANALYZER_RUN_E2E")),
              "Set MRMANALYZER_RUN_E2E=true to run the e2e smoke test")
  skip_if_not_installed("rmarkdown")
  skip_if_not_installed("plotly")

  td <- withr::local_tempdir()
  res <- mrmRunExample(results_dir = td, open = FALSE)

  out_dir <- res$output_directory
  expect_equal(normalizePath(out_dir), normalizePath(td))

  files <- list.files(out_dir, recursive = TRUE)
  # Two HTML reports
  htmls <- grep("\\.html$", files, value = TRUE)
  expect_length(htmls, 2)
  # Processed-matrix RDS + xlsx; stats xlsx
  expect_true(any(grepl("\\.RDS$",   files)))
  expect_true(any(grepl("\\.xlsx$",  files)))
  expect_true(any(grepl("_stats\\.xlsx$", files)))
})

# The same workflow with the report renders off: no pandoc, no plotly, so it
# runs everywhere the package is checked. This is the one test that drives
# runMRManalyzeR() end to end - config in, workbook read, matrix built,
# outputs written - rather than a step in isolation.
test_that("mrmRunExample builds the outputs without rendering reports", {
  td  <- withr::local_tempdir()
  res <- suppressWarnings(suppressMessages(
    mrmRunExample(results_dir = td, open = FALSE, render = FALSE)))

  expect_equal(normalizePath(res$output_directory), normalizePath(td))

  files <- list.files(td, recursive = TRUE)
  expect_true(any(grepl("[.]RDS$", files)))
  expect_true(any(grepl("[.]xlsx$", files)))
  expect_length(grep("[.]html$", files), 0)   # nothing rendered

  # The config it ran from is written beside the outputs, pointing at them.
  expect_true("example_config.yml" %in% files)
  cfg <- loadConfig(file.path(td, "example_config.yml"))
  expect_equal(cfg$project$paths$result_dir, td)
  expect_false(cfg$project$stats_report$execute)

  # And the dataset came back, not just the files.
  de <- res$datasetExperiment
  expect_s4_class(de, "DatasetExperiment")
  expect_gt(ncol(as.data.frame(de$data)), 0)
  expect_gt(nrow(as.data.frame(de$data)), 0)
  expect_true("Compound" %in% colnames(as.data.frame(de$variable_meta)))
})
