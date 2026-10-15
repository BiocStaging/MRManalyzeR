# readPeakMatrix(): the vendor-neutral reader. Takes any sample x analyte
# table, either orientation, from a data frame or an xlsx sheet, and returns
# samples as rows with the compound names kept exactly as written.

.pm_df <- function() {
  data.frame(
    Sample = c("S1", "S2", "S3"),
    `12,13-DiHOME` = c(10, 20, 30),
    `PGE2 (free)`  = c(1, 2, 3),
    check.names = FALSE, stringsAsFactors = FALSE)
}

test_that("a samples-in-rows table comes back as it went in", {
  m <- readPeakMatrix(.pm_df())

  expect_equal(rownames(m), c("S1", "S2", "S3"))
  # Brackets and commas survive: make.names() here would break the match
  # against feature_metadata later.
  expect_equal(colnames(m), c("12,13-DiHOME", "PGE2 (free)"))
  expect_equal(m[["12,13-DiHOME"]], c(10, 20, 30))
  expect_true(is.numeric(m[[1]]))
})

test_that("samples in columns are transposed", {
  wide <- data.frame(
    Compound = c("12,13-DiHOME", "PGE2 (free)"),
    S1 = c(10, 1), S2 = c(20, 2), S3 = c(30, 3),
    check.names = FALSE, stringsAsFactors = FALSE)

  m <- readPeakMatrix(wide, orientation = "samples_cols")
  expect_equal(rownames(m), c("S1", "S2", "S3"))
  expect_equal(colnames(m), c("12,13-DiHOME", "PGE2 (free)"))
  expect_equal(m[["PGE2 (free)"]], c(1, 2, 3))
})

test_that("the identifier column can be named, and must exist", {
  df <- .pm_df()
  df$run_id <- df$Sample
  df$Sample <- NULL                       # identifier is no longer first

  m <- readPeakMatrix(df, id_col = "run_id")
  expect_equal(rownames(m), c("S1", "S2", "S3"))
  expect_error(readPeakMatrix(df, id_col = "nope"), "not a column")
})

test_that("non-detect tokens become NA instead of text", {
  df <- .pm_df()
  df[["12,13-DiHOME"]] <- c("#N/A", "20", "n.d.")

  m <- readPeakMatrix(df)
  expect_true(is.numeric(m[["12,13-DiHOME"]]))
  expect_equal(m[["12,13-DiHOME"]], c(NA, 20, NA))
})

test_that("duplicate samples and single-column tables are refused", {
  df <- .pm_df()
  df$Sample <- c("S1", "S1", "S3")
  expect_error(readPeakMatrix(df), "duplicate sample identifier")
  expect_error(readPeakMatrix(data.frame(Sample = "S1")),
               "at least one measurement column")
})

test_that("it reads the sheet of a workbook", {
  skip_if_not_installed("openxlsx")
  td   <- withr::local_tempdir()
  xlsx <- file.path(td, "matrix.xlsx")
  openxlsx::write.xlsx(list(matrix_data = .pm_df()), xlsx)

  # openxlsx rewrites spaces in column names on read, so the values are
  # checked by position rather than by name here.
  m <- readPeakMatrix(xlsx)
  expect_equal(rownames(m), c("S1", "S2", "S3"))
  expect_equal(unname(m[[1]]), c(10, 20, 30))
  expect_equal(unname(m[[2]]), c(1, 2, 3))
})
