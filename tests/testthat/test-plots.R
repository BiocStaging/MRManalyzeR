# The plot layer: each function returns a ggplot built from the data it was
# given, refuses input it cannot use, and returns NULL rather than erroring
# where there is nothing to draw. The bundled synthetic dataset is used so
# the plots are built from a realistic shape.

.de_plots <- function() {
  de <- loadDataset(system.file("extdata", "example_synthetic.RDS",
                                package = "MRManalyzeR"))
  subsetDataset(de, conditions = list(Sample_type = "Sample"))
}

.stats_plots <- function() {
  runStats(.de_plots(), list(
    sig_threshold = 0.05,
    comparisons = list(enabled = TRUE, entries = list(
      list(name = "PBS_vs_HDM", method = "welch",
           compare = list(factor = "Treatment", levels = c("PBS", "HDM"))))),
    correlations = list(enabled = TRUE, entries = list(
      list(name = "pairs", methods = "pearson", subsets = list(list()))))))
}

test_that("plotBoxplot draws a compound by group, as box or bar", {
  de <- .de_plots()
  # The points carry a `text` aesthetic for the plotly tooltips in the
  # reports, which ggplot2 does not know and warns about on every build.
  expect_s3_class(suppressWarnings(
    plotBoxplot(de, "Analyte_02", "Treatment")), "ggplot")
  expect_s3_class(suppressWarnings(
    plotBoxplot(de, "Analyte_02", "Treatment", type = "bar",
                error_bar = "SE")), "ggplot")
  expect_error(plotBoxplot(data.frame(a = 1), "A", "B"), "DatasetExperiment")
})

test_that("plotQQ draws one compound, transformed or not", {
  de <- .de_plots()
  expect_s3_class(plotQQ(de, "Analyte_02"), "ggplot")
  expect_s3_class(plotQQ(de, "Analyte_02", transform = "log2"), "ggplot")
  expect_error(plotQQ(de, "Analyte_02", transform = "cuberoot"), "arg")
})

test_that("plotDrift needs a numeric injection order", {
  de <- .de_plots()
  expect_s3_class(suppressWarnings(plotDrift(de, "Analyte_02")), "ggplot")

  sm <- as.data.frame(de$sample_meta)
  sm$Injection_order <- paste0("run", seq_len(nrow(sm)))
  de$sample_meta <- sm
  expect_error(plotDrift(de, "Analyte_02"), "numeric")
})

test_that("the PCA plots take the fit and the dataset", {
  de  <- .de_plots()
  pca <- mrmRunPCA(de, transform = "log2")

  expect_s3_class(mrmPlotPCA(pca, de, colour_by = "Treatment"), "ggplot")
  expect_s3_class(
    mrmPlotPCA(pca, de, colour_by = "Treatment", ellipse = "group"), "ggplot")
  expect_s3_class(plotLoadings(pca, de), "ggplot")
  expect_s3_class(plotLoadings(pca, de, colour_by = "Enzymatic_pathway"),
                  "ggplot")

  # A fit that found nothing to decompose draws nothing, and says so by
  # returning NULL rather than failing.
  expect_null(plotLoadings(list(pr = NULL), de))
})

test_that("plotVolcano takes the stats table or the whole list", {
  st <- .stats_plots()
  expect_s3_class(plotVolcano(st, comparison = "PBS_vs_HDM"), "ggplot")
  expect_s3_class(plotVolcano(st$stats, comparison = "PBS_vs_HDM"), "ggplot")
})

test_that("plotCorrelation takes the correlations table, and says when it cannot", {
  st <- .stats_plots()
  p  <- plotCorrelation(st, name = "pairs")
  expect_s3_class(p, "ggplot")

  # The dataset is not what this one takes - the message has to say so.
  expect_error(plotCorrelation(.de_plots()), "correlations table")
  expect_error(plotCorrelation("pairs"), "correlations table")

  # Nothing to draw is NULL, not an error.
  expect_null(plotCorrelation(st$correlations[0, ]))
})
