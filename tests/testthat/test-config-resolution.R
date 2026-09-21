test_that(".resolve_datatype() reads the enabled source block", {
  pmp = list(
    skyline_data = list(enabled = TRUE,  signal_filter = "LOD"),
    tl_data      = list(enabled = FALSE, datatype = "ng/mL"),
    matrix_data  = list(enabled = FALSE))

  # The disabled tl_data block must not win just because it is listed first.
  expect_equal(.resolve_datatype(pmp, list()), "LOD")

  pmp$skyline_data$enabled = FALSE
  pmp$tl_data$enabled      = TRUE
  expect_equal(.resolve_datatype(pmp, list()), "ng/mL")

  pmp$tl_data$enabled     = FALSE
  pmp$matrix_data         = list(enabled = TRUE, datatype = "Area")
  expect_equal(.resolve_datatype(pmp, list()), "Area")
})

test_that("skyline names its runs by signal_filter, not a second key", {
  pmp = list(skyline_data = list(enabled = TRUE, signal_filter = "LOQ"))
  expect_equal(.resolve_datatype(pmp, list()), "LOQ")

  # signal_filter: False means "apply no floor". That is a valid choice but
  # not a label - naming the outputs "FALSE" would be worse than falling back.
  pmp = list(skyline_data = list(enabled = TRUE, signal_filter = FALSE),
             datatype     = "Response")
  expect_equal(.resolve_datatype(pmp, list()), "Response")
})

test_that(".resolve_datatype() honours the older layouts", {
  # PeakMatrixProcessing level, as configs written before the move have it.
  pmp = list(datatype = "Response",
             tl_data  = list(enabled = TRUE))
  expect_equal(.resolve_datatype(pmp, list()), "Response")

  # paths: level, older still.
  expect_equal(.resolve_datatype(list(tl_data = list(enabled = TRUE)),
                                 list(datatype = "Area")), "Area")
})

test_that(".resolve_datatype() warns rather than silently choosing", {
  # Two different answers in the config: the source block wins, but quietly
  # renaming every output file on the strength of that is not acceptable.
  pmp = list(datatype = "Area",
             tl_data  = list(enabled = TRUE, datatype = "ng/mL"))
  expect_warning(dt <- .resolve_datatype(pmp, list()), "set in both")
  expect_equal(dt, "ng/mL")

  # Same value in both places is not a conflict.
  pmp$datatype = "ng/mL"
  expect_silent(.resolve_datatype(pmp, list()))

  # Absent everywhere: defaulting to "Area" renames the outputs, so it warns.
  expect_warning(dt <- .resolve_datatype(list(tl_data = list(enabled = TRUE)),
                                         list()),
                 "defaulting")
  expect_equal(dt, "Area")
})

test_that(".resolve_datatype() keeps a vector for the datatype loop", {
  expect_equal(
    .resolve_datatype(list(tl_data = list(enabled = TRUE,
                                          datatype = c("Area", "Response"))),
                      list()),
    c("Area", "Response"))
  expect_equal(
    .resolve_datatype(list(skyline_data = list(enabled = TRUE,
                                               signal_filter = c("LOD", "LOQ"))),
                      list()),
    c("LOD", "LOQ"))
})

test_that(".narrow_datatype() pins the enabled block to one value", {
  # readSkyline() indexes feature_metadata by signal_filter, so a two-element
  # list surviving into the run would be used as a column name.
  pmp = list(skyline_data = list(enabled = TRUE,
                                 signal_filter = c("LOD", "LOQ")))
  got = .narrow_datatype(pmp, "LOQ")
  expect_equal(got$skyline_data$signal_filter, "LOQ")
  expect_equal(got$datatype, "LOQ")

  pmp = list(tl_data = list(enabled = TRUE,
                            datatype = c("Area", "Response")))
  expect_equal(.narrow_datatype(pmp, "Area")$tl_data$datatype, "Area")

  pmp = list(matrix_data = list(enabled = TRUE, datatype = "Area"))
  expect_equal(.narrow_datatype(pmp, "Area")$matrix_data$datatype, "Area")
})

test_that(".resolve_meta_col() finds columns mangled by make.names()", {
  # struct stores `S-group` as `S.group`, so a config naming the spreadsheet
  # heading must still resolve or the plot silently loses its colouring.
  df = data.frame(Compound = "A", S.group = "x", pubchem.KEGG_id = "1",
                  stringsAsFactors = FALSE)

  expect_equal(.resolve_meta_col("S-group", df),        "S.group")
  expect_equal(.resolve_meta_col("pubchem/KEGG_id", df), "pubchem.KEGG_id")
  expect_equal(.resolve_meta_col("Compound", df),        "Compound")

  expect_null(.resolve_meta_col("Nope", df))
  expect_null(.resolve_meta_col(NULL, df))
  expect_null(.resolve_meta_col("", df))
})

test_that(".resolve_meta_col() finds a heading from its stored spelling", {
  df <- data.frame(`Chrom-Batch` = 1, check.names = FALSE)
  expect_equal(.resolve_meta_col("Chrom.Batch", df), "Chrom-Batch")
  both <- data.frame(`S-group` = 1, `S group` = 2, check.names = FALSE)
  expect_null(.resolve_meta_col("S.group", both))   # ambiguous
})

test_that("a stored dataset with no imputation record follows replace_MVs", {
  de <- struct::DatasetExperiment(
    data          = data.frame(A = c(1, 2), row.names = c("s1", "s2")),
    sample_meta   = data.frame(x = 1:2, row.names = c("s1", "s2")),
    variable_meta = data.frame(Compound = "A", row.names = "A"))
  expect_false(.values_imputed(de, FALSE))
  expect_false(.values_imputed(de, NULL))
  expect_true(.values_imputed(de, 0.5))
  vm <- as.data.frame(de$variable_meta)
  vm$n_imputed <- 0L
  de$variable_meta <- vm
  expect_false(.values_imputed(de, 0.5))   # recorded: nothing was filled
})

test_that(".resolve_meta_col() prefers an exact match over the mangled one", {
  # A workbook holding both spellings must not be silently redirected.
  df = data.frame(check.names = FALSE, `S-group` = "raw", S.group = "clean",
                  stringsAsFactors = FALSE)
  expect_equal(.resolve_meta_col("S-group", df), "S-group")
})

test_that(".read_sheet() names the key and the sheets actually present", {
  skip_if_not_installed("openxlsx")

  tmp = tempfile(fileext = ".xlsx")
  on.exit(unlink(tmp), add = TRUE)
  openxlsx::write.xlsx(list(samples = data.frame(Name = "S1")), tmp)

  expect_equal(.read_sheet(tmp, "samples", "sample_meta_tab")$Name, "S1")

  # A misconfigured sheet name is the most likely failure now that it is
  # configurable, so the error has to name the key and list the alternatives.
  err = expect_error(.read_sheet(tmp, "sample_metadata", "sample_meta_tab"))
  expect_match(conditionMessage(err), "sample_meta_tab")
  expect_match(conditionMessage(err), "samples")
})

test_that(".narrow_datatype() leaves a disabled signal filter off", {
  # The datatype came from elsewhere; writing it into signal_filter would
  # switch on a floor the config deliberately declined.
  pmp = list(skyline_data = list(enabled = TRUE, signal_filter = FALSE))
  expect_false(.narrow_datatype(pmp, "Response")$skyline_data$signal_filter)
})

test_that("a feature grouping column the data does not have is reported", {
  ids <- paste0("S", 1:3)
  de <- struct::DatasetExperiment(
    data = data.frame(A = c(1, 2, 3), row.names = ids),
    sample_meta = data.frame(g = rep("a", 3), row.names = ids),
    variable_meta = data.frame(Compound = "A", Class = "x", row.names = "A"))
  st <- list(pca     = list(loadings = list(color_by = "Enzymatic_pathway")),
             heatmap = list(group_features_by = "Class"))
  expect_warning(.check_feature_cols(st, de), "Enzymatic_pathway")

  st$pca$loadings$color_by <- "Class"
  expect_silent(.check_feature_cols(st, de))
})

test_that("injection order is checked only when the data-quality report runs", {
  cfg <- list(project = list(data_quality_report = list(execute = TRUE)))
  expect_equal(.order_head(cfg), "Injection_order")

  cfg$project$data_quality_report$injection_order_head <- "Inj"
  expect_equal(.order_head(cfg), "Inj")

  cfg$project$data_quality_report$execute <- FALSE
  expect_null(.order_head(cfg))

  # A legacy UVA_report: block is read the way the report reads it.
  legacy <- list(project = list(UVA_report = list(
    execute = TRUE, injection_order_head = "Run")))
  expect_equal(.order_head(legacy), "Run")
})

test_that("combine reads a panel's measured values when they are stored", {
  td <- withr::local_tempdir()
  a <- file.path(td, "a.RDS")
  b <- file.path(td, "b.RDS")
  file.create(a, b, file.path(td, "a_measured.RDS"))
  p <- suppressMessages(.prefer_measured(c(A = a, B = b)))
  expect_equal(unname(p), c(file.path(td, "a_measured.RDS"), b))
  expect_named(p, c("A", "B"))
})

test_that("a correlation on fewer than three shared samples is left out", {
  X <- cbind(a = c(1, 2, NA, NA, 5),
             b = c(2, 4, NA, 7, NA),
             c = c(1, 3, 2, 5, 4))
  out <- .pairwise_cor_long(X, "pearson", "t", "all")
  ab  <- out$feature_a == "b" & out$feature_b == "a"
  ca  <- out$feature_a == "c" & out$feature_b == "a"
  expect_equal(out$n[ab], 2)
  expect_true(is.na(out$estimate[ab]))
  expect_false(is.na(out$estimate[ca]))
})
