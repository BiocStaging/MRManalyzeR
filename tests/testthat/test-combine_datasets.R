# Tests for combineDatasets() — merging multiple panels by Sample_ID.
# Builds two synthetic struct::DatasetExperiment RDS files in tempdir(),
# runs the merge, and checks the row/column outcomes.

.make_de <- function(samples, features, name = "x", values = NULL){
  M <- if(is.null(values)) {
    matrix(seq_len(length(samples) * length(features)),
           nrow = length(samples), ncol = length(features),
           dimnames = list(samples, features))
  } else values
  smeta <- data.frame(Sample_ID  = samples,
                      Sample_type = ifelse(grepl("^QC", samples), "QC", "Sample"),
                      Group      = rep(c("A", "B"), length.out = length(samples)),
                      row.names  = samples,
                      stringsAsFactors = FALSE)
  vmeta <- data.frame(Compound = features,
                      Class    = "lipid",
                      row.names = features,
                      stringsAsFactors = FALSE)
  struct::DatasetExperiment(
    name          = name,
    data          = as.data.frame(M),
    sample_meta   = smeta,
    variable_meta = vmeta
  )
}

.write_de_rds <- function(de, path){
  saveRDS(de, path)
  path
}

test_that("intersects samples by Sample_ID and rbinds variable_meta", {
  td   <- withr::local_tempdir()
  p1 <- .write_de_rds(.make_de(c("S1", "S2", "S3"),         c("A", "B")),
                      file.path(td, "panel1.RDS"))
  p2 <- .write_de_rds(.make_de(c("S2", "S3", "S4"),         c("C", "D")),
                      file.path(td, "panel2.RDS"))
  out <- combineDatasets(c(panel1 = p1, panel2 = p2))

  # Common samples = {S2, S3}; combined features = {A, B, C, D}
  expect_equal(sort(rownames(out$data)), c("S2", "S3"))
  expect_equal(sort(colnames(out$data)), c("A", "B", "C", "D"))
  expect_equal(nrow(out$variable_meta), 4)
})

test_that("prefix_features tags column names with the dataset", {
  td <- withr::local_tempdir()
  p1 <- .write_de_rds(.make_de(c("S1", "S2"), c("X", "Y")),
                      file.path(td, "panel1.RDS"))
  p2 <- .write_de_rds(.make_de(c("S1", "S2"), c("X", "Y")),
                      file.path(td, "panel2.RDS"))
  out <- combineDatasets(c(GOM = p1, SPM = p2), prefix_features = TRUE)

  expect_setequal(colnames(out$data),
                  c("GOM__X", "GOM__Y", "SPM__X", "SPM__Y"))
})

test_that("qc_remap renames listed QCs and drops unmapped ones", {
  td <- withr::local_tempdir()
  p1 <- .write_de_rds(.make_de(c("S1", "QC1", "QC2"),       c("A")),
                      file.path(td, "panel1.RDS"))
  p2 <- .write_de_rds(.make_de(c("S1", "QC_a", "QC_b"),     c("B")),
                      file.path(td, "panel2.RDS"))

  out <- combineDatasets(
    c(p1 = p1, p2 = p2),
    qc_remap = list(p1 = c("QC1" = "QC_pool"),
                    p2 = c("QC_a" = "QC_pool"))
  )
  # QC2 and QC_b are unmapped → dropped. QC1 / QC_a → QC_pool, surviving
  # the intersect across panels along with S1.
  expect_setequal(rownames(out$data), c("S1", "QC_pool"))
})

test_that("drop_samples filter removes named IDs from the merge", {
  td <- withr::local_tempdir()
  p1 <- .write_de_rds(.make_de(c("S1", "S2", "S3"), c("A")),
                      file.path(td, "panel1.RDS"))
  p2 <- .write_de_rds(.make_de(c("S1", "S2", "S3"), c("B")),
                      file.path(td, "panel2.RDS"))
  out <- combineDatasets(c(p1 = p1, p2 = p2), drop_samples = "S2")
  expect_setequal(rownames(out$data), c("S1", "S3"))
})

test_that("feature_meta_cols defaults to intersection across inputs", {
  td <- withr::local_tempdir()
  d1 <- .make_de(c("S1", "S2"), c("A"))
  d1$variable_meta$Pathway <- "P1"            # extra col only in panel 1
  d2 <- .make_de(c("S1", "S2"), c("B"))

  p1 <- .write_de_rds(d1, file.path(td, "panel1.RDS"))
  p2 <- .write_de_rds(d2, file.path(td, "panel2.RDS"))
  out <- combineDatasets(c(p1 = p1, p2 = p2))

  expect_true(all(c("Compound", "Class") %in% colnames(out$variable_meta)))
  expect_false("Pathway" %in% colnames(out$variable_meta))
})

test_that("prefix_features='auto' only prefixes colliding names", {
  td <- withr::local_tempdir()
  # X collides; Y and Z are unique
  p1 <- .write_de_rds(.make_de(c("S1", "S2"), c("X", "Y")),
                      file.path(td, "panel1.RDS"))
  p2 <- .write_de_rds(.make_de(c("S1", "S2"), c("X", "Z")),
                      file.path(td, "panel2.RDS"))
  out <- combineDatasets(c(GOM = p1, SPM = p2), prefix_features = "auto")

  expect_setequal(colnames(out$data), c("GOM__X", "Y", "SPM__X", "Z"))
})

test_that("errors when no samples are common across inputs", {
  td <- withr::local_tempdir()
  p1 <- .write_de_rds(.make_de(c("S1", "S2"), c("A")),
                      file.path(td, "panel1.RDS"))
  p2 <- .write_de_rds(.make_de(c("S3", "S4"), c("B")),
                      file.path(td, "panel2.RDS"))
  expect_error(combineDatasets(c(p1 = p1, p2 = p2)),
               "No samples in common")
})

test_that("the combine run keeps panel names and writes a measured merge", {
  td  <- withr::local_tempdir()
  ids <- c("S1", "S2", "S3")
  imp  <- .make_de(ids, c("A", "B"))
  meas <- imp
  d <- as.data.frame(meas$data); d[1, 1] <- NA; meas$data <- d
  .write_de_rds(imp,  file.path(td, "pa.RDS"))
  .write_de_rds(meas, file.path(td, "pa_measured.RDS"))
  .write_de_rds(.make_de(ids, c("C", "D")), file.path(td, "pb.RDS"))
  cfg <- list(
    combine = list(
      datasets        = list(list(path = file.path(td, "pa.RDS")),
                             list(path = file.path(td, "pb.RDS"))),
      prefix_features = TRUE,
      sample_id_col   = "Sample_ID",
      output_stub     = file.path(td, "comb")),
    stats_report = list(execute = FALSE))
  yml <- file.path(td, "comb.yml")
  yaml::write_yaml(cfg, yml)
  res <- suppressMessages(runMRManalyzeRCombine(yml))

  # Untagged panels are named after the file they were exported to, whether
  # or not their measured copy was read.
  expect_setequal(colnames(res$datasetExperiment$data),
                  c("pa__A", "pa__B", "pb__C", "pb__D"))
  expect_identical(colnames(res$measured$data),
                   colnames(res$datasetExperiment$data))
  expect_false(anyNA(res$datasetExperiment$data))
  expect_true(anyNA(res$measured$data))
  expect_true(file.exists(file.path(td, "comb_measured.RDS")))
  expect_true("matrix_measured" %in%
                openxlsx::getSheetNames(file.path(td, "comb.xlsx")))
})
