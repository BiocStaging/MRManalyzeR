# Tests for plotGroupHeatmap() - the group-mean heatmap whose opacity
# carries each sample's contribution to the mean it is drawn from.
#
# The encoding is the thing worth testing: a group whose samples all move
# together must come out uniformly opaque, and a group whose mean rests on one
# sample must come out with one opaque cell and the rest at the floor. Those
# two cases are constructed explicitly rather than asserted about real data.

.de_for_heatmap <- function(){
  de <- loadDataset(system.file("extdata", "example_synthetic.RDS",
                                 package = "MRManalyzeR"))
  subsetDataset(de, conditions = list(Sample_type = "Sample"))
}

test_that("blocked compounds put their names between the bar and the cells", {
  de <- .de_for_heatmap()
  p  <- plotGroupHeatmap(de, color_samples_by = "Treatment",
                         group_features_by = "Enzymatic_pathway")
  expect_s3_class(p, "ggplot")
  # The names are a text layer now, and the axis text is blank.
  expect_true(any(vapply(p$layers,
                         function(l) inherits(l$geom, "GeomText"), logical(1))))
  expect_s3_class(p$theme$axis.text.y, "element_blank")
  # The columns reserved for them sit between the annotation bar and the cells.
  lv <- levels(p$layers[[1]]$data$sample)
  expect_true(all(c(".ann", ".lab1", ".gapLeft") %in% lv))
  expect_lt(match(".ann", lv), match(".lab1", lv))
  expect_lt(max(grep("^\\.lab", lv)), match(".gapLeft", lv))

  # Without a class column there is nowhere to put them: back to the axis.
  q <- plotGroupHeatmap(de, color_samples_by = "Treatment")
  expect_false(inherits(q$theme$axis.text.y, "element_blank"))

  # hide still drops them entirely.
  r <- plotGroupHeatmap(de, color_samples_by = "Treatment",
                        group_features_by = "Enzymatic_pathway",
                        feature_labels = "hide")
  expect_false(any(vapply(r$layers,
                          function(l) inherits(l$geom, "GeomText"), logical(1))))
})

test_that("returns a ggplot and honours the blocking arguments", {
  de <- .de_for_heatmap()

  p <- plotGroupHeatmap(de, color_samples_by = "Treatment")
  expect_s3_class(p, "ggplot")

  p2 <- plotGroupHeatmap(de, color_samples_by = "Treatment",
                           group_features_by = "Enzymatic_pathway")
  expect_s3_class(p2, "ggplot")

  # one cell per sample x feature, not one per group x feature: the whole point
  # is that individual samples are still drawn. The cells are the first layer;
  # later layers are the two annotation bars.
  n_feat <- ncol(as.data.frame(de$data))
  expect_equal(nrow(p$layers[[1]]$data),
               nrow(as.data.frame(de$data)) * n_feat)
})

test_that("alpha reaches the floor and the ceiling for the right samples", {
  de <- .de_for_heatmap()
  p  <- plotGroupHeatmap(de, color_samples_by = "Treatment", alpha_floor = 0.25)
  cells <- p$layers[[1]]$data

  expect_gte(min(cells$alpha, na.rm = TRUE), 0.25)
  expect_lte(max(cells$alpha, na.rm = TRUE), 1)

  # the floor must actually be reachable, or it is not doing anything
  expect_true(any(abs(cells$alpha - 0.25) < 1e-8))
})

# Two-group toy panels. The scaling is across ALL samples, so which sample
# "contributes" depends on where its group's mean sits relative to the global
# mean - not on which sample has the largest raw value. The next two tests pin
# down both sides of that, because it is the part of the encoding easiest to
# reason about wrongly.
.toy_de <- function(a_values, b_values, seed = 42L){
  set.seed(seed)
  M <- matrix(0, nrow = 8, ncol = 3,
              dimnames = list(paste0("S", 1:8), paste0("F", 1:3)))
  M[1:4, ] <- a_values
  M[5:8, ] <- b_values
  M <- M + matrix(rnorm(24, sd = 0.01), nrow = 8)   # break exact ties
  sm <- data.frame(Treatment = rep(c("A", "B"), each = 4),
                   row.names = rownames(M))
  vm <- data.frame(Compound = colnames(M), row.names = colnames(M))
  struct::DatasetExperiment(data = as.data.frame(M), sample_meta = sm,
                            variable_meta = vm, name = "toy",
                            description = "toy")
}

test_that("the sample that sets a group's mean is opaque, the rest are not", {
  # A = one sample high, three just above the global centre. S1's deviation is
  # 7x the others', so it outvotes them and A's mean comes out positive - S1 is
  # what produced it.
  de <- .toy_de(c(10, 1, 1, 1), rep(2, 4))
  d  <- plotGroupHeatmap(de, color_samples_by = "Treatment",
                           alpha_floor = 0.2)$layers[[1]]$data

  a_carrier <- mean(d$alpha[d$sample == "S1"])
  a_rest    <- mean(d$alpha[d$sample %in% c("S2", "S3", "S4")])

  expect_gt(a_carrier, a_rest)
  expect_equal(a_rest, 0.2, tolerance = 1e-6)   # they pull the other way
  expect_gt(a_carrier, 0.9)                     # it alone carries the block
})

test_that("an extreme sample outvoted by its own group fades", {
  # The mirror case, and the one that caught me out writing these tests: A has
  # one sample far up and three far DOWN, so the three set the group mean and
  # the extreme sample is the one opposing it. Largest raw value does not mean
  # largest contribution.
  de <- .toy_de(c(10, 0, 0, 0), rep(4, 4))
  d  <- plotGroupHeatmap(de, color_samples_by = "Treatment",
                           alpha_floor = 0.2)$layers[[1]]$data

  a_extreme  <- mean(d$alpha[d$sample == "S1"])
  a_majority <- mean(d$alpha[d$sample %in% c("S2", "S3", "S4")])

  expect_gt(a_majority, a_extreme)
  expect_equal(a_extreme, 0.2, tolerance = 1e-6)
})

test_that("degenerate input returns NULL rather than erroring", {
  M  <- matrix(1, nrow = 4, ncol = 3,
               dimnames = list(paste0("S", 1:4), paste0("F", 1:3)))
  sm <- data.frame(Treatment = rep(c("A", "B"), each = 2),
                   row.names = rownames(M))
  vm <- data.frame(Compound = colnames(M), row.names = colnames(M))
  de <- struct::DatasetExperiment(data = as.data.frame(M), sample_meta = sm,
                                  variable_meta = vm, name = "flat",
                                  description = "flat")

  # every compound is zero-variance, so nothing is left to scale
  expect_null(plotGroupHeatmap(de, color_samples_by = "Treatment"))
})

test_that("an unknown grouping column fails with a useful message", {
  de <- .de_for_heatmap()
  expect_error(plotGroupHeatmap(de, color_samples_by = "NoSuchColumn"),
               "not a sample_meta column")
})
