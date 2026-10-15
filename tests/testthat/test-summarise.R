# summariseGroups() and summariseFeatures(): the tables behind the boxplots
# and the "compounds excluded" section of the report.

.de_sum <- function() {
  ids <- paste0("S", 1:6)
  struct::DatasetExperiment(
    data = data.frame(A = c(1, 2, 3, 11, 12, 13),
                      B = c(5, 5, NA, 5, 5, 5),
                      row.names = ids),
    sample_meta = data.frame(Treatment = rep(c("Ctrl", "Trt"), each = 3),
                             row.names = ids),
    variable_meta = data.frame(Compound = c("A", "B"),
                               row.names = c("A", "B")))
}

test_that("group summaries give n, mean, sd and se per compound per group", {
  gs <- summariseGroups(.de_sum(), "Treatment")

  expect_setequal(colnames(gs), c("Compound", "group", "n", "mean", "sd", "se"))
  expect_equal(nrow(gs), 4)                    # 2 compounds x 2 groups

  a_ctrl <- gs[gs$Compound == "A" & gs$group == "Ctrl", ]
  expect_equal(a_ctrl$n, 3L)
  expect_equal(a_ctrl$mean, 2)
  expect_equal(a_ctrl$sd, 1)
  expect_equal(a_ctrl$se, 1 / sqrt(3))

  # A missing value is dropped from n, so se uses the count that was summarised.
  b_ctrl <- gs[gs$Compound == "B" & gs$group == "Ctrl", ]
  expect_equal(b_ctrl$n, 2L)
  expect_equal(b_ctrl$sd, 0)
})

test_that("compounds can be restricted, and bad input is named", {
  gs <- summariseGroups(.de_sum(), "Treatment", compounds = "B")
  expect_equal(unique(gs$Compound), "B")

  expect_error(summariseGroups(.de_sum(), "Nope"), "not a sample_meta column")
  expect_error(summariseGroups(.de_sum(), "Treatment", compounds = "Z"), "Z")
  expect_error(summariseGroups("not a dataset", "Treatment"),
               "DatasetExperiment")
})

test_that("the feature summary marks what was kept and what was dropped", {
  removed <- data.frame(Compound = c("X", "Y"),
                        Comment  = c("All values below LOD threshold",
                                     "No values in the included samples"),
                        stringsAsFactors = FALSE)
  fs <- summariseFeatures(.de_sum(), removed)

  expect_equal(nrow(fs), 4)
  expect_equal(fs$Status[fs$Compound == "A"], "included")
  expect_equal(fs$Status[fs$Compound == "X"], "excluded")
  expect_match(fs$Comment[fs$Compound == "Y"], "No values")

  # Nothing removed is the normal case, and must not error.
  expect_equal(nrow(summariseFeatures(.de_sum())), 2)
  expect_equal(nrow(summariseFeatures(.de_sum(), removed[0, ])), 2)
})
