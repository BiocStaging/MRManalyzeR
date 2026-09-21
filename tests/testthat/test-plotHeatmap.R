test_that("plotHeatmap keeps gaps as gaps and still clusters", {
  skip_if_not_installed("pheatmap")
  ids <- paste0("S", 1:6)
  m <- matrix(c(5, 7, 6, 9, 8, 7,
                3, 4, 5, 3, 6, 4,
                10, 12, 11, 15, 13, 14,
                2, 3, 2, 4, 3, 5), nrow = 6,
              dimnames = list(ids, paste0("F", 1:4)))
  m[2, 1] <- NA
  m[5, 3] <- NA
  de <- struct::DatasetExperiment(
    data          = as.data.frame(m),
    sample_meta   = data.frame(G = rep(c("a", "b"), 3), row.names = ids),
    variable_meta = data.frame(Compound = colnames(m),
                               row.names = colnames(m)))
  p <- plotHeatmap(de, color_samples_by = "G", cluster_rows = TRUE)
  expect_s3_class(p, "ggplot")
})
