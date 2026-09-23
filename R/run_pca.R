#' Pre-process a sample x feature matrix and run PCA
#'
#' Standard metabolomics-style preprocessing pipeline:
#' \enumerate{
#'   \item Drop all-NA features (always). Optionally apply a stricter
#'         feature missing-value filter (`feat_na_max`).
#'   \item Optionally drop samples with too many missing values (`sample_na_max`).
#'   \item Impute remaining NAs below each feature's lowest measured value
#'         (random draws by default - see `impute`).
#'   \item Optional transform (log2(x+1) or sqrt(x)).
#'   \item Drop zero-variance features.
#'   \item `prcomp(center, scale.)` - mean-centre and (optionally) autoscale.
#' }
#'
#' Returns the prcomp object plus a human-readable trace of the steps that
#' were actually applied - reports embed this so the user can see at a
#' glance what transformation went into the PCA.
#'
#' @param X A `struct::DatasetExperiment` (its `data` is used), or a numeric
#'   matrix / data frame of samples x features. May contain NAs.
#' @param impute `"gaussian"` (default), `"frac_min"`, `"min"`, `"half_min"`
#'   or `"none"`. `"gaussian"` fills each NA with its own random value between
#'   0 and `impute_frac` x the feature's minimum (a normal centred at half
#'   that ceiling, SD a sixth of it - see [imputeMissing()]); `"frac_min"`,
#'   `"min"` and `"half_min"` fill every NA with `impute_frac`, 1 or 0.5 x the
#'   minimum. Zeros count as missing for every option except `"none"`.
#' @param impute_frac Fraction of the per-feature minimum used by
#'   `"gaussian"` and `"frac_min"` (default 0.5).
#' @param transform one of `"none"`, `"log2"`, `"sqrt"` (default `"log2"`).
#' @param center logical, passed to [stats::prcomp()].
#' @param scale logical, passed to [stats::prcomp()] as `scale.`.
#' @param feat_na_max Maximum allowed fraction of NA values per feature (0-1).
#'   Features exceeding this threshold are dropped before imputation.
#'   Default `1` retains the original behaviour (only all-NA features are dropped).
#'   Set e.g. `0.5` to also drop features missing in more than 50 percent of samples.
#' @param sample_na_max Maximum allowed fraction of NA values per sample (0-1).
#'   Samples exceeding this threshold are dropped before imputation.
#'   Default `1` applies no sample filtering (original behaviour).
#'   Set e.g. `0.8` to drop samples missing more than 80 percent of features.
#' @param seed Random seed for `impute = "gaussian"`, so a rerun gives the
#'   same PCA; `NULL` uses the current random stream.
#' @return A list with elements
#'   \describe{
#'     \item{`pr`}{the `prcomp` object, or `NULL` if PCA could not be fit.}
#'     \item{`X_clean`}{the matrix actually fed to `prcomp`.}
#'     \item{`n_samples`, `n_features`}{dimensions after cleaning.}
#'     \item{`dropped_all_na`, `dropped_feat_filter`, `dropped_sample_filter`,
#'           `dropped_zero_var`}{counts of features/samples dropped at each step.}
#'     \item{`steps`}{character vector of human-readable steps applied.}
#'   }
#' @examples
#' X <- matrix(abs(rnorm(60, 100, 20)), nrow = 12,
#'             dimnames = list(paste0("S", 1:12), paste0("F", 1:5)))
#' res <- runPCA(X, transform = "log2")
#' res$steps
#' @family QC check
#' @export
runPCA = function(X,
                            impute        = c("gaussian", "frac_min", "min",
                                              "half_min", "none"),
                            impute_frac   = 0.5,
                            transform     = c("log2", "sqrt", "none"),
                            center        = TRUE,
                            scale         = TRUE,
                            feat_na_max   = 1,
                            sample_na_max = 1,
                            seed          = 42){

  impute    = match.arg(impute)
  transform = match.arg(transform)
  steps     = character(0)

  feat_na_max   = as.numeric(feat_na_max)
  sample_na_max = as.numeric(sample_na_max)

  if(methods::is(X, "DatasetExperiment")) X = X$data
  X = as.matrix(X)

  # 1. Drop all-NA features (always) + optional stricter feature filter
  na_frac_feat = colMeans(is.na(X))
  all_na       = na_frac_feat >= 1
  dropped_all_na = sum(all_na)
  X = X[, !all_na, drop = FALSE]
  if(dropped_all_na > 0)
    steps = c(steps, sprintf("Dropped %d all-NA features.", dropped_all_na))

  dropped_feat_filter = 0L
  if(feat_na_max < 1){
    na_frac_feat2   = colMeans(is.na(X))
    drop_f          = na_frac_feat2 > feat_na_max
    dropped_feat_filter = sum(drop_f)
    X = X[, !drop_f, drop = FALSE]
    steps = c(steps, sprintf(
      "Feature missing filter (max %g%% NA): dropped %d feature(s).",
      feat_na_max * 100, dropped_feat_filter))
  } else {
    steps = c(steps, "No feature missing-value filter (feat_na_max = 1).")
  }

  # 2. Optional sample filter
  dropped_sample_filter = 0L
  if(sample_na_max < 1){
    na_frac_samp = rowMeans(is.na(X))
    drop_s       = na_frac_samp > sample_na_max
    dropped_sample_filter = sum(drop_s)
    X = X[!drop_s, , drop = FALSE]
    steps = c(steps, sprintf(
      "Sample missing filter (max %g%% NA): dropped %d sample(s).",
      sample_na_max * 100, dropped_sample_filter))
  } else {
    steps = c(steps, "No sample missing-value filter (sample_na_max = 1).")
  }

  # 3. Impute remaining NAs below each feature's lowest measured value
  if(impute != "none"){
    fac = switch(impute,
                 min      = 1,
                 half_min = 0.5,
                 as.numeric(impute_frac))
    random = identical(impute, "gaussian")
    X = if(random && !is.null(seed))
      withr::with_seed(seed, .fill_below_min(X, fac, TRUE)) else
      .fill_below_min(X, fac, random)
    steps = c(steps, if(random)
      sprintf(paste0("Imputed NAs with random values between 0 and %g x the ",
                     "per-feature minimum (normal, mean at half that ceiling, ",
                     "SD a sixth of it; seed %s)."),
              fac, if(is.null(seed)) "none" else seed)
      else sprintf("Imputed NAs with per-feature %s.",
                   switch(impute, min = "minimum", half_min = "half-minimum",
                          sprintf("%g x minimum", fac))))
  } else {
    steps = c(steps, "No NA imputation.")
  }

  # 3. Transform
  if(transform == "log2"){
    X[X < 0 & !is.na(X)] = 0
    X = log2(X + 1)
    steps = c(steps, "log2(x + 1) transform.")
  } else if(transform == "sqrt"){
    X[X < 0 & !is.na(X)] = 0
    X = sqrt(X)
    steps = c(steps, "sqrt(x) transform.")
  } else {
    steps = c(steps, "No transform.")
  }

  # 4. Drop zero-variance / still-NA features
  vars = apply(X, 2, stats::var, na.rm = TRUE)
  keep = is.finite(vars) & vars > 0
  dropped_zero_var = sum(!keep)
  X = X[, keep, drop = FALSE]
  if(dropped_zero_var > 0)
    steps = c(steps, sprintf("Dropped %d zero-variance features.", dropped_zero_var))

  # 5. PCA
  if(center) steps = c(steps, "Mean-centred per feature.")
  if(scale)  steps = c(steps, "Autoscaled (divide by per-feature SD).")

  # The handler returns the condition rather than writing to `steps` from
  # inside it: prcomp() cannot itself return something inheriting from
  # "error", so testing the result is unambiguous, and the audit trail is
  # then appended in the frame that owns it.
  if(ncol(X) >= 2 && nrow(X) >= 3){
    fit = tryCatch(stats::prcomp(X, center = center, scale. = scale),
                   error = function(e) e)
    if(inherits(fit, "error")){
      steps = c(steps, sprintf("prcomp() failed: %s.",
                               conditionMessage(fit)))
      pr = NULL
    } else {
      pr = fit
    }
  } else {
    steps = c(steps, sprintf("Insufficient data: n_samples=%d, n_features=%d.",
                             nrow(X), ncol(X)))
    pr = NULL
  }

  list(pr = pr, X_clean = X,
       n_samples             = nrow(X),
       n_features            = ncol(X),
       dropped_all_na        = dropped_all_na,
       dropped_feat_filter   = dropped_feat_filter,
       dropped_sample_filter = dropped_sample_filter,
       dropped_zero_var      = dropped_zero_var,
       steps                 = steps)
}
