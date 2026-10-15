#' Adjust a dataset to true sample concentrations
#'
#' A calibration curve returns the concentration in the vial, and it does so
#' based on the calibrant vial: it was built at the calibrant's
#' internal-standard concentration and assumes the sample carries the same.
#' This corrects that assumption, given the same volume of calibration and
#' study sample being injected into the instrument:
#'
#' \deqn{conc \times (cal\_vol / sample\_vol) \times (sample\_IS / cal\_IS)}
#'
#' Ultimately leading to
#' \deqn{(sample\_IS / sample\_vol) / (cal\_IS / cal\_vol)}
#' which is the ratio of internal-standard *concentrations*, sample vial over
#' calibrant vial. More vial volume per unit of IS means a more dilute IS, a
#' smaller IS peak, an inflated analyte/IS response ratio and therefore a curve
#' that reads high - so the factor falls below 1 and scales it back down.
#'
#' When `starting_vol_col` is supplied, a second step scales back to the
#' concentration in the original, undried sample
#' (\eqn{conc \times sample\_vol / starting\_vol}).
#'
#' @param de A `struct::DatasetExperiment`.
#' @param sample_vol_col,cal_vol_col,sample_IS_col,cal_IS_col `sample_meta`
#'   columns holding, respectively, the sample's reconstitution volume, the
#'   calibration-standard volume, the IS volume added to the sample, and the IS
#'   volume added to the calibration standard.
#' @param starting_vol_col `sample_meta` column holding the pre-dry-down
#'   starting volume, or `FALSE` to skip the starting-volume correction.
#' @return `de` with `data` scaled to concentration per sample.
#' @family peak-matrix processing
#' @examples
#' de <- struct::DatasetExperiment(
#'   data = data.frame(PGE2 = c(10, 20), PGD2 = c(30, 40),
#'                     row.names = c("S1", "S2")),
#'   sample_meta = data.frame(sample_volume_uL = c(100, 50),
#'                            cal_vol_uL       = c(100, 100),
#'                            sample_IS_vol_uL = c(10, 10),
#'                            cal_IS_vol_uL    = c(10, 10),
#'                            row.names = c("S1", "S2")),
#'   variable_meta = data.frame(Compound = c("PGE2", "PGD2"),
#'                              row.names = c("PGE2", "PGD2")))
#' adjustConcentration(de)$data
#' @export
adjustConcentration <- function(de,
                                sample_vol_col   = "sample_volume_uL",
                                cal_vol_col      = "cal_vol_uL",
                                sample_IS_col    = "sample_IS_vol_uL",
                                cal_IS_col       = "cal_IS_vol_uL",
                                starting_vol_col = FALSE){
  .check_de(de, "adjustConcentration")

  smeta <- as.data.frame(de$sample_meta)
  req  <- c(sample_vol_col, cal_vol_col, sample_IS_col, cal_IS_col)
  miss <- setdiff(req, colnames(smeta))
  if(length(miss))
    stop(sprintf(
      "[adjustConcentration] sample_meta is missing required column(s): %s",
      paste(miss, collapse = ", ")))

  fact <- (smeta[[cal_vol_col]]   / smeta[[sample_vol_col]]) *
         (smeta[[sample_IS_col]] / smeta[[cal_IS_col]])

  n_missing <- sum(is.na(fact))
  if(n_missing)
    message(sprintf(
      "[adjustConcentration] %d sample(s) have missing volume/IS metadata and will be NA: %s",
      n_missing, paste(head(rownames(smeta)[is.na(fact)], 5), collapse = ", ")))

  x <- as.data.frame(de$data) * fact

  if(!isFALSE(starting_vol_col)){
    if(!starting_vol_col %in% colnames(smeta))
      stop(sprintf("[adjustConcentration] sample_meta has no '%s' column.",
                   starting_vol_col))
    sv_fact <- smeta[[sample_vol_col]] / smeta[[starting_vol_col]]
    if(sum(is.na(sv_fact)))
      message(sprintf(
        "[adjustConcentration] %d sample(s) have missing %s/%s and will be NA.",
        sum(is.na(sv_fact)), sample_vol_col, starting_vol_col))
    x <- x * sv_fact
  }

  de$data <- x
  de
}
