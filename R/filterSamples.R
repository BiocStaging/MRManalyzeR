#' Filter injections by missing-value fraction
#'
#' Drops study samples in which too large a proportion of features is `NA`.
#' Where [filterFeatures()] asks a question about the biology, this asks one
#' about the run: an injection that missed most of the panel is usually a
#' failure rather than a finding.
#'
#' **Run this after [filterFeatures()], not before.** Judging an injection
#' against a feature list that still contains everything the instrument barely
#' saw gives a meaningless denominator - every sample looks poor because the
#' feature set is padded with noise. Once the undetected features are gone, the
#' fraction answers the question worth asking: of the features that are real,
#' how many did this injection get?
#'
#' **Only study samples are eligible for removal.** A QC or blank that fails the
#' threshold is reported but kept, because dropping a QC can leave a batch
#' without one and [correctBatch()] would then have nothing to anchor to. If a
#' QC really did fail, exclude it deliberately in the sample metadata rather
#' than letting a threshold do it silently.
#'
#' The fraction computed for every row is recorded in `sample_meta` as
#' `na_frac`, including for the QCs and blanks that were exempt, so a failing
#' QC is visible afterwards.
#'
#' @param de A `struct::DatasetExperiment`.
#' @param max_na Numeric in `[0, 1]`. A study sample is dropped when the
#'   proportion of `NA` features exceeds this.
#' @param sample_col `sample_meta` column identifying the sample type.
#' @param sample_labels Values of `sample_col` marking study samples - the only
#'   rows this function will remove.
#' @return `de` with failing study samples removed from `data` and
#'   `sample_meta`, and an `na_frac` column added to `sample_meta`.
#' @family peak-matrix processing
#' @seealso [filterFeatures()], which should run first.
#' @examples
#' de <- struct::DatasetExperiment(
#'   data = data.frame(
#'     f1 = c(10, 11, NA, 12),
#'     f2 = c(20, 21, NA, 22),
#'     f3 = c(30, 31, NA, 33),
#'     row.names = c("S1", "S2", "S3_failed", "QC1")),
#'   sample_meta = data.frame(
#'     Sample_type = c("Sample", "Sample", "Sample", "QC"),
#'     row.names   = c("S1", "S2", "S3_failed", "QC1")),
#'   variable_meta = data.frame(
#'     Compound  = c("f1", "f2", "f3"),
#'     row.names = c("f1", "f2", "f3")))
#'
#' # S3_failed is dropped; a QC at the same missingness would only be reported
#' filterSamples(de, max_na = 0.5)
#' @export
filterSamples <- function(de, max_na = 0.8,
                         sample_col = "Sample_type",
                         sample_labels = "Sample"){
  .check_de(de, "filterSamples")

  if(!is.numeric(max_na) || length(max_na) != 1L ||
     is.na(max_na) || max_na < 0 || max_na > 1)
    stop("filterSamples: max_na must be a single number in [0, 1].")

  X    <- as.data.frame(de$data)
  smeta <- as.data.frame(de$sample_meta)
  if(!nrow(X) || !ncol(X)) return(de)

  na_frac <- rowMeans(is.na(X))

  eligible <- rep(TRUE, nrow(X))
  if(!is.null(sample_col)){
    if(!sample_col %in% colnames(smeta)){
      fmt <- paste0("filterSamples: '%s' not in sample_meta - every ",
                    "injection is eligible for removal, QCs and ",
                    "blanks included.")
      warning(sprintf(fmt, sample_col))
    } else {
      eligible <- as.character(smeta[[sample_col]]) %in% sample_labels
    }
  }

  failing <- na_frac > max_na
  to_drop <- failing & eligible
  exempt  <- failing & !eligible

  .warn_overwrite(smeta, "na_frac", "filterSamples", "sample_meta")
  smeta$na_frac <- round(na_frac, 4)
  de$sample_meta <- smeta

  if(any(exempt)){
    fmt <- paste0("filterSamples: %d non-study injection(s) exceed ",
                  "max_na = %g and were KEPT (%s). Exclude them in ",
                  "the sample metadata if they really failed.")
    which_ones <- paste(utils::head(rownames(X)[exempt], 5), collapse = ", ")
    warning(sprintf(fmt, sum(exempt), max_na, which_ones))
  }

  if(!any(to_drop)){
    message(sprintf("[filterSamples] all %d injections retained (max_na = %g).",
                    nrow(X), max_na))
    return(de)
  }

  if(all(to_drop[eligible]))
    stop("filterSamples: max_na removed every study sample.")

  message(sprintf("[filterSamples] dropped %d of %d study injections above max_na = %g",
                  sum(to_drop), sum(eligible), max_na))
  message("   ", paste(sprintf("%s (%.0f%% missing)",
                               rownames(X)[to_drop],
                               100 * na_frac[to_drop]), collapse = ", "))

  keep <- !to_drop
  struct::DatasetExperiment(
    data          = X[keep, , drop = FALSE],
    sample_meta   = de$sample_meta[keep, , drop = FALSE],
    variable_meta = as.data.frame(de$variable_meta)
  )
}
