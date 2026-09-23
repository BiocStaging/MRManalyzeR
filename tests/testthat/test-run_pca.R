# Tests for runPCA() — the PCA preprocessing pipeline used by
# both reports. Pure function on a numeric matrix, easy to unit test.

# Helper: deterministic synthetic matrix with planted structure.
.synth_matrix <- function(n_samples = 20, n_features = 12, seed = 1L){
  set.seed(seed)
  M <- matrix(rlnorm(n_samples * n_features, meanlog = 4, sdlog = 0.6),
              nrow = n_samples, ncol = n_features)
  rownames(M) <- paste0("S", seq_len(n_samples))
  colnames(M) <- paste0("F", seq_len(n_features))
  M
}

test_that("returns expected list structure", {
  M <- .synth_matrix()
  res <- runPCA(M)

  expect_named(res, c("pr", "X_clean", "n_samples", "n_features",
                      "dropped_all_na", "dropped_feat_filter",
                      "dropped_sample_filter", "dropped_zero_var", "steps"))
  expect_s3_class(res$pr, "prcomp")
  expect_equal(res$n_samples, 20)
  expect_equal(res$n_features, 12)
  expect_type(res$steps, "character")
  expect_gt(length(res$steps), 0)
})

test_that("drops all-NA features regardless of feat_na_max", {
  M <- .synth_matrix()
  M[, 3] <- NA_real_
  res <- runPCA(M, feat_na_max = 1)
  expect_equal(res$dropped_all_na, 1)
  expect_equal(res$n_features, 11)
})

test_that("feat_na_max drops features above threshold", {
  M <- .synth_matrix()
  # Make column 5 have 70% NAs
  M[1:14, 5] <- NA_real_
  res <- runPCA(M, feat_na_max = 0.5, impute = "min")
  expect_equal(res$dropped_feat_filter, 1)
  expect_false("F5" %in% colnames(res$X_clean))
})

test_that("sample_na_max drops samples above threshold", {
  M <- .synth_matrix()
  # Make row 1 have 90% NAs
  M[1, 1:11] <- NA_real_
  res <- runPCA(M, sample_na_max = 0.5, impute = "min")
  expect_equal(res$dropped_sample_filter, 1)
  expect_false("S1" %in% rownames(res$X_clean))
})

test_that("frac_min imputation uses per-feature minimum × frac", {
  M <- .synth_matrix()
  feat_min <- min(M[, 1])
  M[5, 1] <- NA_real_
  res <- runPCA(M, impute = "frac_min", impute_frac = 0.2,
                          transform = "none", center = FALSE, scale = FALSE)
  # Imputed value should equal min × 0.2
  expect_equal(res$X_clean["S5", "F1"], feat_min * 0.2)
})

test_that("transform options apply correctly", {
  M <- matrix(c(1, 2, 4, 8), nrow = 2, ncol = 2,
              dimnames = list(c("S1", "S2"), c("F1", "F2")))

  res_none <- runPCA(M, transform = "none", impute = "none",
                               center = FALSE, scale = FALSE)
  res_log2 <- runPCA(M, transform = "log2", impute = "none",
                               center = FALSE, scale = FALSE)
  res_sqrt <- runPCA(M, transform = "sqrt", impute = "none",
                               center = FALSE, scale = FALSE)

  expect_equal(unname(res_none$X_clean), unname(M))
  expect_equal(unname(res_log2$X_clean), unname(log2(M + 1)))
  expect_equal(unname(res_sqrt$X_clean), unname(sqrt(M)))
})

test_that("insufficient data yields NULL pr without error", {
  M <- matrix(1, nrow = 2, ncol = 2,
              dimnames = list(c("S1", "S2"), c("F1", "F2")))
  res <- runPCA(M)
  expect_null(res$pr)                        # too few samples
  expect_true(any(grepl("Insufficient", res$steps)))
})

test_that("zero-variance features are dropped after transform", {
  M <- .synth_matrix()
  M[, 4] <- 5                                # constant column
  res <- runPCA(M, transform = "none", center = FALSE, scale = FALSE)
  expect_gte(res$dropped_zero_var, 1)
  expect_false("F4" %in% colnames(res$X_clean))
})

test_that("gaussian imputation fills below the minimum, reproducibly", {
  M <- .synth_matrix()
  M[1:3, 2] <- NA_real_
  run <- function() runPCA(M, impute = "gaussian", impute_frac = 0.2,
                           transform = "none", center = FALSE, scale = FALSE)
  r1 <- run()
  r2 <- run()
  fills   <- r1$X_clean[1:3, 2]
  ceiling <- 0.2 * min(M[, 2], na.rm = TRUE)
  expect_true(all(fills > 0 & fills < ceiling))
  expect_identical(r1$X_clean, r2$X_clean)
})

test_that("gaussian imputation treats zeros as missing", {
  M <- .synth_matrix()
  M[1, 2] <- 0
  M[2, 2] <- NA_real_
  r <- runPCA(M, impute = "gaussian", impute_frac = 0.2, transform = "none",
              center = FALSE, scale = FALSE)
  fills <- r$X_clean[1:2, 2]
  expect_true(all(fills > 0))   # a zero minimum would have filled with 0
})
