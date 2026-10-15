# adjustConcentration(): the calibrant-to-sample internal-standard correction.
#
# The factor is (cal_vol / sample_vol) * (sample_IS / cal_IS), per injection.
# S1 carries 100/50 * 10/20 = 1 and S2 carries 100/200 * 40/20 = 1, so two
# different vial setups land on the same correction; S3 doubles it.
# starting_vol_col adds sample_vol / starting_vol on top.

.de_conc <- function() {
  ids <- c("S1", "S2", "S3")
  struct::DatasetExperiment(
    data = data.frame(A = c(10, 10, 10), B = c(4, 4, 4), row.names = ids),
    sample_meta = data.frame(
      sample_volume_uL = c(50, 200, 50),
      cal_vol_uL       = c(100, 100, 100),
      sample_IS_vol_uL = c(10, 40, 20),
      cal_IS_vol_uL    = c(20, 20, 20),
      starting_vol_uL  = c(100, 100, 100),
      row.names        = ids),
    variable_meta = data.frame(Compound = c("A", "B"),
                               row.names = c("A", "B")))
}

test_that("the correction is the ratio of internal-standard concentrations", {
  d <- as.data.frame(adjustConcentration(.de_conc())$data)
  expect_equal(d$A, c(10, 10, 20))
  expect_equal(d$B, c(4, 4, 8))
})

test_that("starting_vol_col takes the result back to the original sample", {
  # sample_vol / starting_vol = 0.5, 2 and 0.5 on top of the factors above.
  d <- as.data.frame(
    adjustConcentration(.de_conc(), starting_vol_col = "starting_vol_uL")$data)
  expect_equal(d$A, c(5, 20, 10))
})

test_that("a missing volume gives NA for that injection and says so", {
  de <- .de_conc()
  sm <- as.data.frame(de$sample_meta)
  sm$cal_IS_vol_uL[2] <- NA
  de$sample_meta <- sm

  expect_message(out <- adjustConcentration(de), "missing volume")
  d <- as.data.frame(out$data)
  expect_true(all(is.na(d[2, ])))
  expect_false(anyNA(d[c(1, 3), ]))
})

test_that("a missing column stops the run and names it", {
  de <- .de_conc()
  sm <- as.data.frame(de$sample_meta)
  sm$cal_vol_uL <- NULL
  de$sample_meta <- sm

  expect_error(adjustConcentration(de), "cal_vol_uL")
  expect_error(adjustConcentration(.de_conc(), starting_vol_col = "nope"),
               "nope")
  expect_error(adjustConcentration(list(a = 1)), "DatasetExperiment")
})
