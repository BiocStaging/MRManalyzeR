# newConfig() writes the file a first run starts from. What matters is that
# what it writes parses - including the !concat_path anchors - that the data
# block matches the source asked for, and that it will not quietly replace a
# configuration someone has already edited.

test_that("the minimal study config parses, anchors and all", {
  td  <- withr::local_tempdir()
  cfg <- file.path(td, "study.yml")

  expect_message(out <- newConfig(cfg), "runMRManalyzeR")
  expect_equal(out, cfg)
  expect_true(file.exists(cfg))

  p <- loadConfig(cfg)
  expect_equal(p$project$name, "my_study")

  # !concat_path resolved against the &project_dir anchor.
  expect_equal(p$project$paths$data_dir, "C:/path/to/study/data")
  expect_equal(p$project$paths$result_dir, "C:/path/to/study/results")

  pmp <- p$project$PeakMatrixProcessing
  expect_true(pmp$execute)
  expect_equal(pmp$name_col, "Name")
  expect_equal(pmp$feature_meta_tab, "feature_metadata")

  # Every processing step off until asked for: a first run reads the data and
  # does nothing else to it.
  expect_false(pmp$blank_filter)
  expect_false(pmp$filter_features$enabled)
  expect_false(pmp$batch_correction)
  expect_false(pmp$normalize)
  expect_false(pmp$replace_MVs)

  expect_true(p$project$stats_report$execute)
  expect_true(p$project$data_quality_report$execute)
})

test_that("each source writes its own data block, and only its own", {
  td <- withr::local_tempdir()

  for(src in c("targetlynx", "skyline", "matrix")){
    f <- file.path(td, paste0(src, ".yml"))
    suppressMessages(newConfig(f, source = src))
    cfg <- loadConfig(f)
    pmp <- cfg$project$PeakMatrixProcessing

    blocks <- c(targetlynx = "tl_data", skyline = "skyline_data",
                matrix = "matrix_data")
    expect_true(pmp[[blocks[[src]]]]$enabled)
    expect_false(any(setdiff(blocks, blocks[[src]]) %in% names(pmp)))

    # And the output naming resolves from it, which is what the block is for.
    expect_equal(.resolve_datatype(pmp, cfg$project$paths),
                 if(identical(src, "skyline")) "LOD" else "Area")
  }

  expect_error(newConfig(file.path(td, "x.yml"), source = "mzml"), "arg")
})

test_that("a directory gets project_config.yml inside it", {
  td  <- withr::local_tempdir()
  out <- suppressMessages(newConfig(td))
  expect_equal(out, file.path(td, "project_config.yml"))
  expect_true(file.exists(out))
})

test_that("an edited configuration is not replaced by accident", {
  td  <- withr::local_tempdir()
  cfg <- file.path(td, "mine.yml")
  suppressMessages(newConfig(cfg))
  writeLines("project: {name: edited}", cfg)

  expect_error(newConfig(cfg), "overwrite = TRUE")
  expect_equal(readLines(cfg), "project: {name: edited}")

  suppressMessages(newConfig(cfg, overwrite = TRUE))
  expect_equal(loadConfig(cfg)$project$name, "my_study")
})

test_that("the combine template writes the merge form", {
  td  <- withr::local_tempdir()
  cfg <- file.path(td, "combine.yml")
  suppressMessages(newConfig(cfg, template = "combine"))

  p <- loadConfig(cfg)
  expect_null(p$project)
  expect_length(p$combine$datasets, 2)
  expect_equal(vapply(p$combine$datasets, function(d) d$tag, character(1)),
               c("A", "B"))
  expect_equal(p$combine$sample_id_col, "Sample_ID")
  expect_true(p$stats_report$execute)
})

test_that("minimal = FALSE copies the bundled reference, comments and all", {
  td <- withr::local_tempdir()

  for(tmpl in c("study", "combine")){
    f <- file.path(td, paste0(tmpl, "_full.yml"))
    suppressMessages(newConfig(f, template = tmpl, minimal = FALSE))
    ref <- system.file("extdata",
                       if(identical(tmpl, "combine"))
                         "example_combine_config.yml" else "example_config.yml",
                       package = "MRManalyzeR")
    expect_equal(readLines(f), readLines(ref))
  }
})

test_that("path has to be one name", {
  expect_error(newConfig(c("a.yml", "b.yml")), "single file or directory")
  expect_error(newConfig(""), "single file or directory")
  expect_error(newConfig(42), "single file or directory")
})
