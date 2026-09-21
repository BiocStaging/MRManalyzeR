# correctBatch(): median-ratio batch correction.
#
# Two batches, each with five QCs and two samples. The values are chosen so
# every expectation below can be checked by hand:
#   A  b1 QCs c(NA, NA, 1, 2, 3) -> median 1 counting NA as 0 (2 if NAs were
#      dropped); b2 QCs all 2; grand median of all ten QCs is 2.
#      So b1 is divided by 1/2 (doubled) and b2 is unchanged.
#   B  no missing values: b1 QCs 4, b2 QCs 8, grand median 6.
#      Samples at 6 (b1) and 12 (b2) both correct to 9.
#   C  b2 QCs c(NA, NA, NA, 5, 5) -> median 0. b1 alone would be rescaled
#      (median 10 against a grand median of 7.5), but a compound that fails
#      in one batch is corrected in none.
#   D  every median negative (-2, -4, grand -3): the ratio would be positive,
#      but a negative median is no usable reference.

mk_bc <- function() {
  ids <- c(paste0("Q", 1:5), "S1", "S2", paste0("Q", 6:10), "S3", "S4")
  struct::DatasetExperiment(
    data = data.frame(
      A = c(NA, NA, 1, 2, 3, 10, 11, 2, 2, 2, 2, 2, 10, 11),
      B = c(4, 4, 4, 4, 4, 6, 6, 8, 8, 8, 8, 8, 12, 12),
      C = c(10, 10, 10, 10, 10, 7, 7, NA, NA, NA, 5, 5, 7, 7),
      D = c(-2, -2, -2, -2, -2, -3, -3, -4, -4, -4, -4, -4, -5, -5),
      row.names = ids),
    sample_meta = data.frame(
      Sample_type = rep(c(rep("QC", 5), "Sample", "Sample"), 2),
      Chrom_Batch = rep(c("b1", "b2"), each = 7),
      row.names = ids),
    variable_meta = data.frame(
      Compound = c("A", "B", "C", "D"),
      Correct  = c("YES", " no ", "", "NO"),
      row.names = c("A", "B", "C", "D")))
}

bc <- function(de, ...)
  correctBatch(de, qc_label = "QC", factor_name = "Sample_type",
               batch_head = "Chrom_Batch", ...)

test_that("a missing value counts as zero in the median", {
  d <- as.data.frame(suppressMessages(suppressWarnings(bc(mk_bc())))$data)
  expect_equal(d["S1", "A"], 20)    # b1 median 1 against grand 2: doubled
  expect_equal(d["S3", "A"], 10)    # b2 median 2 = grand: unchanged
  expect_true(is.na(d["Q1", "A"]))  # missing values are not filled in
})

test_that("compounds without missing values get the plain median ratio", {
  d <- as.data.frame(suppressMessages(suppressWarnings(bc(mk_bc())))$data)
  expect_equal(d["S1", "B"], 9)
  expect_equal(d["S3", "B"], 9)
})

test_that("a compound that fails in one batch is corrected in none", {
  expect_warning(out <- suppressMessages(bc(mk_bc())), "not batch-corrected")
  d  <- as.data.frame(out$data)
  vm <- as.data.frame(out$variable_meta)
  # Left as measured in every batch - not rescaled in b1, not NA in b2.
  expect_equal(d[c("S1", "S2", "S3", "S4"), "C"], c(7, 7, 7, 7))
  expect_equal(d["S1", "D"], -3)
  expect_equal(vm["C", "bc_failed_batches"], "b2")
  expect_equal(vm["D", "bc_failed_batches"], "grand median, b1, b2")
  expect_equal(vm["A", "bc_failed_batches"], "")
  expect_equal(vm$batch_corrected, c(TRUE, TRUE, FALSE, FALSE))
})

test_that("each failing compound and its batches are printed", {
  msgs <- testthat::capture_messages(suppressWarnings(bc(mk_bc())))
  expect_match(msgs, "  C: b2", all = FALSE)
  expect_match(msgs, "  D: grand median, b1, b2", all = FALSE)
})

test_that("feature_col corrects only the YES compounds", {
  expect_warning(out <- bc(mk_bc(), feature_col = "Correct"), "no YES/NO")
  d  <- as.data.frame(out$data)
  vm <- as.data.frame(out$variable_meta)
  expect_equal(d["S1", "A"], 20)   # YES
  expect_equal(d["S1", "B"], 6)    # " no ": left as measured
  expect_equal(d["S3", "C"], 7)    # blank: left as measured
  expect_equal(vm$batch_corrected, c(TRUE, FALSE, FALSE, FALSE))
})

test_that("feature_col rejects a missing column or an unreadable flag", {
  expect_error(bc(mk_bc(), feature_col = "Nope"), "not a variable_meta column")

  de <- mk_bc()
  vm <- as.data.frame(de$variable_meta)
  vm$Correct <- c("YES", "maybe", "NO", "NO")
  de$variable_meta <- vm
  expect_error(bc(de, feature_col = "Correct"), "maybe")
})

# A two-batch Skyline workbook: two QCs and two samples per batch. "#N/A" is
# how Skyline writes a non-detect; the reader turns it into NA.
sky_fixture <- function(dir) {
  xlsx <- file.path(dir, "skyline.xlsx")
  samp <- c("Q1", "Q2", "S1", "S2", "Q3", "Q4", "S3", "S4")
  sheet <- data.frame(
    Molecule = c("PGE2", "TXB2"),
    rbind(c("#N/A", "4", "10", "#N/A", "8", "8", "20", "22"),
          c("5", "5", "6", "7", "9", "9", "11", "12")),
    check.names = FALSE)
  colnames(sheet) <- c("Molecule", samp)
  openxlsx::write.xlsx(list(skyline_data = sheet), xlsx)
  list(
    xlsx  = xlsx,
    fdata = data.frame(Processing_name = sheet$Molecule,
                       Compound = sheet$Molecule,
                       Report = "YES", Comment = ""),
    meta  = data.frame(Name = samp, Include = "YES",
                       Sample_type = rep(c("QC", "QC", "Sample", "Sample"), 2),
                       Chrom_Batch = rep(c("b1", "b2"), each = 4)))
}

sky_run <- function(fx, ...) suppressWarnings(suppressMessages(
  processDataset(fx$fdata, fx$meta, xlsx_path = fx$xlsx,
                 data_source = "skyline", signal_filter = FALSE,
                 blank_filter = FALSE, processing_batch = FALSE, ...)))[[1]]

test_that("processDataset corrects batches before it imputes", {
  fx <- sky_fixture(withr::local_tempdir())
  imputed   <- sky_run(fx, batch_correction = TRUE, replace_MVs = 0.5)
  corrected <- sky_run(fx, batch_correction = TRUE, replace_MVs = FALSE)
  # The pipeline's imputed result is the corrected data, imputed afterwards.
  expect_equal(
    as.data.frame(imputed$data),
    as.data.frame(suppressWarnings(imputeMissing(corrected, 0.5))$data))
  # Correcting without imputing keeps the non-detects missing.
  expect_true(anyNA(as.data.frame(corrected$data)))
})

test_that("batch_correction accepts a header as typed or as renamed", {
  fx <- sky_fixture(withr::local_tempdir())
  fx$fdata[["BC-flag"]] <- c("YES", "NO")
  # assembleDataset() stores the heading as BC.flag and says to use that
  # name; the spelling in the workbook has to keep working too.
  for (nm in c("BC-flag", "BC.flag")) {
    vm <- as.data.frame(sky_run(fx, batch_correction = nm)$variable_meta)
    expect_true(vm["PGE2", "batch_corrected"])
    expect_false(vm["TXB2", "batch_corrected"])
  }
})

test_that("processDataset rejects an absent batch_correction column", {
  xlsx <- system.file("extdata", "example_data.xlsx", package = "MRManalyzeR")
  skip_if(!nzchar(xlsx), "Bundled example_data.xlsx not installed")
  fdata <- openxlsx::read.xlsx(xlsx, sheet = "feature_metadata")
  meta  <- openxlsx::read.xlsx(xlsx, sheet = "sample_metadata")
  expect_error(
    suppressMessages(processDataset(fdata, meta, xlsx_path = xlsx,
                                    batch_correction = "Nope")),
    "feature_metadata column")
})

test_that("an injection with no batch label stops the correction", {
  de <- mk_bc()
  sm <- as.data.frame(de$sample_meta)
  sm["S4", "Chrom_Batch"] <- NA
  de$sample_meta <- sm
  expect_error(bc(de), "no 'Chrom_Batch' value: S4")
})

test_that("processDataset batch-corrects before it normalises", {
  fx <- sky_fixture(withr::local_tempdir())
  # Pooled QCs have no protein amount. Normalising first would turn them into
  # NA and leave the correction without a reference.
  fx$meta$protein <- c(NA, NA, 2, 2, NA, NA, 2, 2)
  vm <- as.data.frame(
    sky_run(fx, batch_correction = TRUE, normalize = "protein")$variable_meta)
  # Under the old order every QC is NA by now, both compounds fail, and this
  # is c(FALSE, FALSE).
  expect_equal(vm$batch_corrected, c(TRUE, TRUE))
})
