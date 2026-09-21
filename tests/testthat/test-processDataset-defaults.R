# Defaults processDataset() must honour when a config leaves a key out.
#
# runMRManalyzeR() used to pass `data_tab_names %||% "skyline_data"` whatever
# the source, and processDataset() defaulted to "skyline_data" too, so a
# TargetLynx config that omitted the key - documented as auto-matching the
# 'lcms_data*' sheets - failed looking for a sheet called skyline_data. It also
# passed an absent `normalize:` through as NULL, which reached
# normaliseMatrix() and stopped the run with "argument is of length zero".

tl_run <- function(xlsx, ...) {
  fdata <- openxlsx::read.xlsx(xlsx, sheet = "feature_metadata")
  meta  <- openxlsx::read.xlsx(xlsx, sheet = "sample_metadata")
  suppressWarnings(suppressMessages(
    processDataset(fdata, meta, xlsx_path = xlsx, data_source = "targetlynx",
                   datatype = "Area", snr = 3, bc_header = "Chrom_Batch",
                   blank_head = "Sample_type", ...)))
}

test_that("TargetLynx without data_tab_names reads the lcms_data sheets", {
  xlsx <- system.file("extdata", "example_data.xlsx", package = "MRManalyzeR")
  skip_if(!nzchar(xlsx), "Bundled example_data.xlsx not installed")

  implicit <- tl_run(xlsx, blank_filter = FALSE)
  explicit <- tl_run(xlsx, blank_filter = FALSE, data_tab_names = "lcms_data")
  expect_equal(as.data.frame(implicit$dataset$data),
               as.data.frame(explicit$dataset$data))
})

test_that("NULL normalize, blank_filter and replace_MVs mean off, like FALSE", {
  xlsx <- system.file("extdata", "example_data.xlsx", package = "MRManalyzeR")
  skip_if(!nzchar(xlsx), "Bundled example_data.xlsx not installed")

  off   <- tl_run(xlsx, blank_filter = FALSE, normalize = FALSE,
                  replace_MVs = FALSE)
  nulls <- tl_run(xlsx, blank_filter = NULL, normalize = NULL,
                  replace_MVs = NULL)
  expect_equal(as.data.frame(nulls$dataset$data),
               as.data.frame(off$dataset$data))
})

test_that("Skyline without data_tab_names still reads the skyline_data sheet", {
  td   <- withr::local_tempdir()
  xlsx <- file.path(td, "skyline.xlsx")
  samp <- paste0("S", 1:4)
  sheet <- data.frame(Molecule = c("PGE2", "PGD2", "TXB2"),
                      matrix(c(10, 11, 12, 13,
                               20, 21, 22, 23,
                               30, 31, 32, 33), nrow = 3, byrow = TRUE,
                             dimnames = list(NULL, samp)),
                      check.names = FALSE)
  openxlsx::write.xlsx(list(skyline_data = sheet), xlsx)

  fdata <- data.frame(Processing_name = sheet$Molecule,
                      Compound = sheet$Molecule,
                      Report = "YES", Comment = "")
  meta  <- data.frame(Name = samp, Include = "YES", Sample_type = "Sample")

  run <- function(...) suppressWarnings(suppressMessages(
    processDataset(fdata, meta, xlsx_path = xlsx, data_source = "skyline",
                   signal_filter = FALSE, blank_filter = FALSE,
                   processing_batch = FALSE, ...)))

  expect_equal(as.data.frame(run()$dataset$data),
               as.data.frame(run(data_tab_names = "skyline_data")$dataset$data))
})

test_that("a dropped feature says whether it was never measured or filtered", {
  td   <- withr::local_tempdir()
  xlsx <- file.path(td, "skyline.xlsx")
  samp <- paste0("S", 1:3)
  sheet <- data.frame(
    Molecule = c("PGE2", "PGD2", "TXB2"),
    rbind(c("10", "12", "11"),
          c("#N/A", "#N/A", "#N/A"),
          c("5", "6", "7")),
    check.names = FALSE)
  colnames(sheet) <- c("Molecule", samp)
  openxlsx::write.xlsx(list(skyline_data = sheet), xlsx)

  fdata <- data.frame(Processing_name = sheet$Molecule,
                      Compound = sheet$Molecule, Report = "YES",
                      Comment = "", LOD = c(NA, NA, 100))
  meta  <- data.frame(Name = samp, Include = "YES", Sample_type = "Sample")

  out <- suppressWarnings(suppressMessages(
    processDataset(fdata, meta, xlsx_path = xlsx, data_source = "skyline",
                   signal_filter = "LOD", blank_filter = FALSE,
                   processing_batch = FALSE)))
  ex <- out$excluded_features
  expect_equal(ex$Comment[ex$Compound == "PGD2"],
               "No values in the included samples")
  expect_equal(ex$Comment[ex$Compound == "TXB2"],
               "All values below LOD threshold")
})

test_that("a zero is read as a non-detect", {
  td   <- withr::local_tempdir()
  xlsx <- file.path(td, "skyline.xlsx")
  samp <- paste0("S", 1:3)
  sheet <- data.frame(Molecule = c("PGE2", "PGD2"),
                      matrix(c(10, 0, 12,
                               0,  0, 0), nrow = 2, byrow = TRUE,
                             dimnames = list(NULL, samp)),
                      check.names = FALSE)
  openxlsx::write.xlsx(list(skyline_data = sheet), xlsx)
  fdata <- data.frame(Processing_name = sheet$Molecule,
                      Compound = sheet$Molecule, Report = "YES", Comment = "")
  meta  <- data.frame(Name = samp, Include = "YES", Sample_type = "Sample")

  out <- suppressWarnings(suppressMessages(
    processDataset(fdata, meta, xlsx_path = xlsx, data_source = "skyline",
                   signal_filter = FALSE, blank_filter = FALSE,
                   processing_batch = FALSE)))
  d <- as.data.frame(out$measured$data)
  expect_equal(sum(is.na(d$PGE2)), 1)
  expect_false(any(d == 0, na.rm = TRUE))
  ex <- out$excluded_features
  expect_equal(ex$Comment[ex$Compound == "PGD2"],
               "No values in the included samples")
})

test_that("processDataset keeps the measured dataset beside the imputed one", {
  xlsx <- system.file("extdata", "example_data.xlsx", package = "MRManalyzeR")
  skip_if(!nzchar(xlsx), "Bundled example_data.xlsx not installed")

  out <- tl_run(xlsx, blank_filter = FALSE, replace_MVs = 0.2)
  expect_equal(dim(out$dataset$data), dim(out$measured$data))
  expect_lt(sum(is.na(as.data.frame(out$dataset$data))),
            sum(is.na(as.data.frame(out$measured$data))))

  plain <- tl_run(xlsx, blank_filter = FALSE)
  expect_equal(as.data.frame(plain$measured$data),
               as.data.frame(plain$dataset$data))
})

test_that("a processing-batch heading renamed on import still works", {
  td   <- withr::local_tempdir()
  xlsx <- file.path(td, "skyline.xlsx")
  samp <- paste0("S", 1:4)
  sheet <- data.frame(Molecule = c("PGE2", "TXB2"),
                      rbind(c("10", "#N/A", "12", "14"),
                            c("5", "6", "#N/A", "8")),
                      check.names = FALSE)
  colnames(sheet) <- c("Molecule", samp)
  openxlsx::write.xlsx(list(skyline_data = sheet), xlsx)

  fdata <- data.frame(Processing_name = sheet$Molecule,
                      Compound = sheet$Molecule, Report = "YES", Comment = "")
  meta  <- data.frame(Name = samp, Include = "YES", Sample_type = "Sample",
                      `Proc-Batch` = c("b1", "b1", "b2", "b2"),
                      check.names = FALSE)

  out <- suppressWarnings(suppressMessages(
    processDataset(fdata, meta, xlsx_path = xlsx, data_source = "skyline",
                   signal_filter = FALSE, blank_filter = FALSE,
                   processing_batch = "Proc-Batch", replace_MVs = 0.5)))
  vm <- as.data.frame(out$dataset$variable_meta)
  expect_true(all(c("impute_max_b1", "impute_max_b2") %in% colnames(vm)))
})
