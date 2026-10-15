# The two columns the output sheets lead with, after the compound name:
#   feature_metadata$Lowest_measured  - the floor the instrument reported
#   removed_features$Failed_step      - which step took the compound out
#
# A small Skyline workbook is built so each exclusion route is exercised:
#   KEEP_A  detected everywhere
#   KEEP_B  detected everywhere, lower values
#   USER_X  curated out in the workbook (Report = NO)
#   SPARSE  measured in one sample of four, so the detection filter takes it
#   GONE    no values at all

.wb_steps <- function(dir) {
  xlsx <- file.path(dir, "steps.xlsx")
  samp <- paste0("S", 1:4)
  sheet <- data.frame(
    Molecule = c("KEEP_A", "KEEP_B", "USER_X", "SPARSE", "GONE", "LONG_DP"),
    rbind(c(100, 110, 120, 130),
          c(5,   6,   7,   8),
          c(50,  50,  50,  50),
          c(9,   NA,  NA,  NA),
          c(NA,  NA,  NA,  NA),
          # the floor here is 1807.4204, and 0.00123456 once scaled down
          c(1807.4204, 9000.5, 12345.678, 2000)),
    check.names = FALSE)
  colnames(sheet) <- c("Molecule", samp)
  openxlsx::write.xlsx(list(skyline_data = sheet), xlsx)
  xlsx
}

.run_steps <- function(xlsx, ...) {
  fdata <- data.frame(
    Processing_name = c("KEEP_A", "KEEP_B", "USER_X", "SPARSE", "GONE",
                        "LONG_DP"),
    Compound        = c("KEEP_A", "KEEP_B", "USER_X", "SPARSE", "GONE",
                        "LONG_DP"),
    Report          = c("YES", "YES", "NO", "YES", "YES", "YES"),
    Comment         = "",
    stringsAsFactors = FALSE)
  meta <- data.frame(Name = paste0("S", 1:4), Include = "YES",
                     Sample_type = "Sample", Treatment = c("A", "A", "B", "B"),
                     stringsAsFactors = FALSE)
  suppressWarnings(suppressMessages(
    processDataset(fdata, meta, xlsx_path = xlsx, data_source = "skyline",
                   signal_filter = FALSE, blank_filter = FALSE,
                   processing_batch = FALSE, ...)))
}

test_that("feature_metadata leads with the lowest value measured", {
  td  <- withr::local_tempdir()
  out <- .run_steps(.wb_steps(td))
  vm  <- as.data.frame(out$dataset$variable_meta)

  expect_equal(colnames(vm)[1:2], c("Compound", "Lowest_measured"))
  expect_equal(vm["KEEP_A", "Lowest_measured"], 100)
  expect_equal(vm["KEEP_B", "Lowest_measured"], 5)

  # Rounded to the nearest whole number: the column is read by a person, and
  # a peak area means nothing at the fourth decimal place. The floor here is
  # 1807.4204.
  expect_equal(vm["LONG_DP", "Lowest_measured"], 1807)
  expect_true(all(vm$Lowest_measured[!is.na(vm$Lowest_measured)] %%
                  1 == 0))

  # Measured, not imputed: filling the gaps must not lower the floor.
  filled <- .run_steps(.wb_steps(td), replace_MVs = 0.5)
  vmf    <- as.data.frame(filled$dataset$variable_meta)
  expect_equal(vmf["KEEP_B", "Lowest_measured"], 5)
  expect_equal(as.data.frame(filled$measured$variable_meta)[["Lowest_measured"]],
               vmf[["Lowest_measured"]])
})

test_that("removed_features leads with the step that took each compound", {
  td  <- withr::local_tempdir()
  out <- .run_steps(.wb_steps(td),
                    filter_features = list(min_frac = 0.5, method = "across",
                                           sample_col = "Sample_type",
                                           sample_labels = "Sample"))
  rf <- out$excluded_features

  expect_equal(colnames(rf)[1:2], c("Compound", "Failed_step"))
  step <- stats::setNames(rf$Failed_step, rf$Compound)

  expect_equal(unname(step["USER_X"]), "User exclude")
  expect_equal(unname(step["GONE"]), "Never measured")
  expect_equal(unname(step["SPARSE"]), "Detection filter")
  expect_false(any(c("KEEP_A", "KEEP_B") %in% rf$Compound))
})

test_that("the signal filter names itself as the step", {
  td <- withr::local_tempdir()
  xlsx <- .wb_steps(td)

  fdata <- data.frame(
    Processing_name = c("KEEP_A", "KEEP_B", "USER_X", "SPARSE", "GONE",
                        "LONG_DP"),
    Compound        = c("KEEP_A", "KEEP_B", "USER_X", "SPARSE", "GONE",
                        "LONG_DP"),
    Report          = "YES", Comment = "",
    LOD             = c(NA, 1000, NA, NA, NA, NA), # KEEP_B is all below its LOD
    stringsAsFactors = FALSE)
  meta <- data.frame(Name = paste0("S", 1:4), Include = "YES",
                     Sample_type = "Sample", stringsAsFactors = FALSE)

  out <- suppressWarnings(suppressMessages(
    processDataset(fdata, meta, xlsx_path = xlsx, data_source = "skyline",
                   signal_filter = "LOD", blank_filter = FALSE,
                   processing_batch = FALSE)))
  rf <- out$excluded_features
  expect_match(rf$Failed_step[rf$Compound == "KEEP_B"], "Signal filter \\(LOD\\)")
})

test_that("nothing removed still gives the column", {
  td <- withr::local_tempdir()
  xlsx <- file.path(td, "clean.xlsx")
  samp <- paste0("S", 1:3)
  sheet <- data.frame(Molecule = c("A", "B"),
                      rbind(c(10, 11, 12), c(20, 21, 22)), check.names = FALSE)
  colnames(sheet) <- c("Molecule", samp)
  openxlsx::write.xlsx(list(skyline_data = sheet), xlsx)

  fdata <- data.frame(Processing_name = c("A", "B"), Compound = c("A", "B"),
                      Report = "YES", Comment = "", stringsAsFactors = FALSE)
  meta  <- data.frame(Name = samp, Include = "YES", Sample_type = "Sample",
                      stringsAsFactors = FALSE)
  out <- suppressWarnings(suppressMessages(
    processDataset(fdata, meta, xlsx_path = xlsx, data_source = "skyline",
                   signal_filter = FALSE, blank_filter = FALSE,
                   processing_batch = FALSE)))

  expect_equal(nrow(out$excluded_features), 0)
  expect_true("Failed_step" %in% colnames(out$excluded_features))
})

test_that("the floor is read off the study samples, not the QCs or blanks", {
  td   <- withr::local_tempdir()
  xlsx <- file.path(td, "qc.xlsx")
  ids  <- c("S1", "S2", "QC1", "B1")
  # Both compounds are lower in the QC and the blank than in any sample, so a
  # floor taken over everything would report those instead.
  sheet <- data.frame(Molecule = c("A", "B"),
                      rbind(c(100, 110, 5, 1),
                            c(50,  60,  7, 8)),
                      check.names = FALSE)
  colnames(sheet) <- c("Molecule", ids)
  openxlsx::write.xlsx(list(skyline_data = sheet), xlsx)

  fdata <- data.frame(Processing_name = c("A", "B"), Compound = c("A", "B"),
                      Report = "YES", Comment = "", stringsAsFactors = FALSE)
  meta  <- data.frame(Name = ids, Include = "YES",
                      Sample_type = c("Sample", "Sample", "QC", "Blank"),
                      stringsAsFactors = FALSE)

  out <- suppressWarnings(suppressMessages(
    processDataset(fdata, meta, xlsx_path = xlsx, data_source = "skyline",
                   signal_filter = FALSE, blank_filter = FALSE,
                   processing_batch = FALSE)))
  vm <- as.data.frame(out$dataset$variable_meta)

  expect_equal(vm["A", "Lowest_measured"], 100)
  expect_equal(vm["B", "Lowest_measured"], 50)

  # The QC and blank values are still in the matrix - they are just not what
  # the floor is read from.
  d <- as.data.frame(out$dataset$data)
  expect_equal(min(d[["A"]], na.rm = TRUE), 1)
})
