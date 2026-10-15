# loadDataset() reads back what a run wrote, from either the RDS or the xlsx.
# Most of it is alignment: the matrix, the sample metadata and the feature
# metadata have to end up describing each other, and where they cannot the
# message has to say which file and which names are at fault rather than
# leaving a silently mismatched dataset to be analysed.

.de_load <- function() {
  ids <- c("S1", "S2", "S3")
  struct::DatasetExperiment(
    name = "fixture",
    data = data.frame(A = c(1, 2, 3), B = c(4, 5, 6), C = c(7, 8, 9),
                      row.names = ids),
    sample_meta = data.frame(Sample_ID = ids, Treatment = c("x", "y", "y"),
                             row.names = ids),
    variable_meta = data.frame(Compound = c("A", "B", "C"),
                               row.names = c("A", "B", "C")))
}

# The three sheets runMRManalyzeR() writes. The matrix sheet carries its
# sample ids in the first column, which is what rowNames = TRUE reads back.
.xlsx_load <- function(dir, file = "run.xlsx", ids = c("S1", "S2", "S3"),
                       mat_ids = ids, compounds = c("A", "B", "C")) {
  path <- file.path(dir, file)
  mat  <- data.frame(Sample = mat_ids,
                     A = c(1, 2, 3)[seq_along(mat_ids)],
                     B = c(4, 5, 6)[seq_along(mat_ids)],
                     C = c(7, 8, 9)[seq_along(mat_ids)],
                     check.names = FALSE)
  openxlsx::write.xlsx(list(
    feature_metadata = data.frame(Compound = compounds, Report = "YES"),
    sample_metadata  = data.frame(Sample_ID = ids, Treatment = "x"),
    matrix           = mat), path)
  path
}

test_that("an RDS comes back as it went in", {
  td <- withr::local_tempdir()
  f  <- file.path(td, "de.RDS")
  saveRDS(.de_load(), f)

  de <- loadDataset(f)
  expect_s4_class(de, "DatasetExperiment")
  expect_equal(as.data.frame(de$data), as.data.frame(.de_load()$data))
  expect_equal(colnames(as.data.frame(de$data)), c("A", "B", "C"))
})

# The realignment path is for RDS files written by older versions of struct,
# where variable_meta and the matrix were separate data frames and could fall
# out of step. Current struct keeps the two in lockstep - assigning different
# row names renames the matrix columns instead - so the fixture here is a
# plain list, which is how loadDataset() reaches those slots either way.
.rds_legacy <- function(dir, file, vm) {
  f <- file.path(dir, file)
  saveRDS(list(
    name = "legacy",
    data = data.frame(A = c(1, 2, 3), B = c(4, 5, 6), C = c(7, 8, 9),
                      row.names = c("S1", "S2", "S3")),
    sample_meta = data.frame(Sample_ID = c("S1", "S2", "S3"),
                             row.names = c("S1", "S2", "S3")),
    variable_meta = vm), f)
  f
}

test_that("variable_meta is put back in step with the matrix columns", {
  td <- withr::local_tempdir()

  # The right compounds, in the wrong order, under row names of their own.
  vm  <- data.frame(Compound = c("C", "A", "B"), LOD = c(3, 1, 2),
                    row.names = c("row_one", "row_two", "row_three"))
  out <- loadDataset(.rds_legacy(td, "drift.RDS", vm))
  ovm <- as.data.frame(out$variable_meta)

  expect_equal(rownames(ovm), c("A", "B", "C"))
  expect_equal(as.character(ovm$Compound), c("A", "B", "C"))
  expect_equal(ovm$LOD, c(1, 2, 3))     # each row brought its own values
  expect_equal(as.data.frame(out$data)[["B"]], c(4, 5, 6))
})

test_that("an RDS that cannot be realigned says which names are at fault", {
  td <- withr::local_tempdir()

  # A matrix column that no Compound row claims.
  expect_error(
    loadDataset(.rds_legacy(td, "lost.RDS",
                            data.frame(Compound = c("A", "B", "Z"),
                                       row.names = c("A", "B", "Z")))),
    "no matching Compound row")

  # Nothing to align on at all.
  expect_error(
    loadDataset(.rds_legacy(td, "nokey.RDS",
                            data.frame(Label = c("A", "B", "C"),
                                       row.names = c("x", "y", "z")))),
    "no Compound column")
})

test_that("an xlsx is rebuilt from its three sheets", {
  td <- withr::local_tempdir()
  de <- loadDataset(.xlsx_load(td))

  expect_s4_class(de, "DatasetExperiment")
  expect_equal(de$name, "run")
  expect_equal(rownames(as.data.frame(de$data)), c("S1", "S2", "S3"))
  expect_equal(colnames(as.data.frame(de$data)), c("A", "B", "C"))
  expect_equal(as.data.frame(de$data)[["B"]], c(4, 5, 6))
  expect_equal(rownames(as.data.frame(de$sample_meta)), c("S1", "S2", "S3"))
  expect_equal(rownames(as.data.frame(de$variable_meta)), c("A", "B", "C"))
})

test_that("the matrix is keyed by sample id, and falls back to row order", {
  td <- withr::local_tempdir()

  # Named, but written in a different order from sample_metadata: the ids win.
  f  <- .xlsx_load(td, "shuffled.xlsx", mat_ids = c("S3", "S1", "S2"))
  d  <- as.data.frame(loadDataset(f)$data)
  expect_equal(rownames(d), c("S1", "S2", "S3"))
  expect_equal(d[["A"]], c(2, 3, 1))     # S1 was row 2 of the sheet

  # No usable ids in the matrix sheet, but one row per sample: take them in
  # the order given rather than refuse to load.
  f <- .xlsx_load(td, "unnamed.xlsx", mat_ids = c("1", "2", "3"))
  expect_equal(rownames(as.data.frame(loadDataset(f)$data)),
               c("S1", "S2", "S3"))
})

test_that("a mismatch names the file and what did not line up", {
  td <- withr::local_tempdir()

  # Fewer matrix rows than samples, and no ids to align on.
  f <- .xlsx_load(td, "short.xlsx", mat_ids = c("1", "2"))
  expect_error(loadDataset(f), "cannot align matrix rows")

  # A compound in the matrix that feature_metadata does not list.
  f <- .xlsx_load(td, "unlisted.xlsx", compounds = c("A", "B", "Z"))
  expect_error(loadDataset(f), "no matching row in feature_metadata")

  # sample_metadata without the id column asked for.
  f <- .xlsx_load(td, "nosid.xlsx")
  expect_error(loadDataset(f, sample_id_col = "Vial"),
               "sample_metadata lacks column 'Vial'")
})

test_that("a workbook that is not one of ours says which sheet is missing", {
  td <- withr::local_tempdir()
  f  <- file.path(td, "other.xlsx")
  openxlsx::write.xlsx(list(Sheet1 = data.frame(a = 1)), f)

  expect_error(loadDataset(f), "missing required sheet")
  expect_error(loadDataset(f), "feature_metadata")
})

test_that("anything else is refused by name", {
  td <- withr::local_tempdir()
  f  <- file.path(td, "peaks.csv")
  writeLines("a,b", f)

  expect_error(loadDataset(f), "unsupported extension")
  expect_error(loadDataset(file.path(td, "nope.RDS")), "file not found")
})
