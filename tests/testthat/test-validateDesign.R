# A stats report that will not run must not be validated: its comparisons are
# never tested, so an error there must not stop a processing run, and warnings
# about levels it will never test bury the warnings that matter. The runner
# renders the report only when execute is True, so that is the condition.

mk_design_de <- function() {
  struct::DatasetExperiment(
    data = data.frame(f1 = c(1, 2, 3, 4, 5, 6), f2 = c(2, 3, 4, 5, 6, 7),
                      row.names = paste0("S", 1:6)),
    sample_meta = data.frame(Sample_type = "Sample",
                             Group = rep(c("a", "b"), each = 3),
                             row.names = paste0("S", 1:6)),
    variable_meta = data.frame(Compound = c("f1", "f2"),
                               row.names = c("f1", "f2")))
}

# One comparison asking for levels the data does not have.
ghost_cfg <- function(execute) {
  list(project = list(
    PeakMatrixProcessing = list(tl_data = list(enabled = TRUE)),
    stats_report = list(
      execute = execute,
      comparisons = list(enabled = TRUE, entries = list(
        list(name = "ghost", method = "welch",
             compare = list(factor = "Group", levels = c("x", "y"))))))))
}

test_that("validateDesign skips comparisons when execute is False", {
  v <- validateDesign(mk_design_de(), ghost_cfg(FALSE))
  expect_false(any(grepl("not present in 'Group'", v$errors)))
})

test_that("validateDesign still checks them when the report runs", {
  v <- validateDesign(mk_design_de(), ghost_cfg(TRUE))
  expect_true(any(grepl("not present in 'Group'", v$errors)))
})

test_that("validateConfig skips comparisons when execute is False", {
  cfg <- ghost_cfg(FALSE)
  # anova is not a valid method for a 2-level comparison
  cfg$project$stats_report$comparisons$entries[[1]]$method <- "anova"
  expect_length(validateConfig(cfg)$errors, 0)

  cfg$project$stats_report$execute <- TRUE
  expect_true(any(grepl("not valid for 2 level", validateConfig(cfg)$errors)))
})

test_that("a comparison switched off with enabled: False is left out", {
  sec <- list(enabled = TRUE, entries = list(
    list(name = "a"), list(name = "b", enabled = FALSE),
    list(name = "c", enabled = TRUE)))
  kept <- .section_entries(sec, "comparisons")
  expect_equal(vapply(kept, function(e) e$name, ""), c("a", "c"))

  cfg <- ghost_cfg(TRUE)
  cfg$project$stats_report$comparisons$entries[[1]]$enabled <- FALSE
  v <- validateDesign(mk_design_de(), cfg)
  expect_false(any(grepl("not present in 'Group'", v$errors)))
})

test_that("a config needs only the data-source block it uses", {
  # No skyline_data: or matrix_data: block at all - absent means off.
  cfg <- list(project = list(
    PeakMatrixProcessing = list(tl_data = list(enabled = TRUE))))
  expect_length(validateConfig(cfg)$errors, 0)

  cfg$project$PeakMatrixProcessing$skyline_data <- list(enabled = TRUE)
  expect_match(validateConfig(cfg)$errors, "exactly one", all = FALSE)
})

test_that("validateDesign does not count an empty batch label as a batch", {
  ids <- paste0("S", 1:6)
  de <- struct::DatasetExperiment(
    data = data.frame(f1 = c(1, 2, 3, 4, 5, 6), row.names = ids),
    sample_meta = data.frame(
      Sample_type = c("QC", "QC", "Sample", "QC", "QC", "Sample"),
      Chrom_Batch = c("b1", "b1", "b1", "b2", "b2", NA),
      row.names = ids),
    variable_meta = data.frame(Compound = "f1", row.names = "f1"))
  v <- validateDesign(de)
  expect_false(any(grepl("'NA'", v$warnings)))
})
