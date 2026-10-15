test_that(".stats_key_df() is well formed", {
  key <- .stats_key_df()

  expect_named(key, c("Sheet", "Column", "Meaning"))
  expect_gt(nrow(key), 30)
  expect_false(any(is.na(unlist(key))))
  expect_false(any(!nzchar(unlist(key))))
  # One definition per (sheet, column) - a duplicate means two different
  # explanations of the same heading, and only one of them can be right.
  expect_equal(anyDuplicated(paste(key$Sheet, key$Column)), 0L)
})

# The key is only useful while it matches what runStats() actually writes.
# These two directions catch the drift: a column defined but never written,
# and a column written but never defined.
test_that(".stats_key_df() matches the schemas runStats() returns", {

  schemas <- list(
    stats         = names(.empty_stats_df()),
    correlations  = names(.empty_corr_df()),
    linear_models = names(.empty_lm_df()),
    ion_ratios    = c("ratio", "numerator", "denominator", "comparison",
                      "sample", "group", "value", "method", "statistic",
                      "p_value", "p_adj"),
    summary       = c("Compound", "group", "n", "mean", "sd", "se"))

  key    <- .stats_key_df()
  generic <- key$Column[key$Sheet == "(all sheets)"]

  for(sheet in names(schemas)){
    defined <- key$Column[key$Sheet == sheet]

    # No phantom definitions.
    expect_equal(setdiff(defined, schemas[[sheet]]), character(0))
    # No undefined headings, once the shared ones are allowed for.
    expect_equal(setdiff(schemas[[sheet]], c(defined, generic)),
                 character(0))
  }
})

test_that("the key sheet is written first and covers only the sheets present", {
  skip_if_not_installed("openxlsx")

  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "stats")
  openxlsx::writeData(wb, "stats", .empty_stats_df())

  .write_key_sheet(wb)

  tmp <- tempfile(fileext = ".xlsx")
  openxlsx::saveWorkbook(wb, tmp, overwrite = TRUE)
  on.exit(unlink(tmp), add = TRUE)

  # worksheetOrder is applied when the workbook is written; the in-memory
  # sheets() listing stays in creation order, so the tab order has to be read
  # back from the file - which is the thing that matters anyway.
  expect_equal(openxlsx::getSheetNames(tmp)[1], "key")

  got <- openxlsx::read.xlsx(tmp, sheet = "key")
  expect_true(all(c("stats", "(all sheets)") %in% got$Sheet))
  # correlations was never written, so its definitions must not appear.
  expect_false("correlations" %in% got$Sheet)
})

test_that("append_stats_xlsx() skips empty tables but still writes the key", {
  skip_if_not_installed("openxlsx")

  tmp <- tempfile(fileext = ".xlsx")
  on.exit(unlink(tmp), add = TRUE)

  append_stats_xlsx(tmp, list(stats         = .empty_stats_df(),
                              correlations  = .empty_corr_df(),
                              linear_models = .empty_lm_df()),
                    group_summary = data.frame(Compound = "A", group = "g",
                                               n = 3L, mean = 1, sd = 0.1,
                                               se = 0.06))

  sh <- openxlsx::sheets(openxlsx::loadWorkbook(tmp))
  expect_equal(sh[1], "key")
  expect_true("summary" %in% sh)
  # Zero-row inputs produce no headers-only tabs.
  expect_false(any(c("stats", "correlations", "linear_models") %in% sh))
})
