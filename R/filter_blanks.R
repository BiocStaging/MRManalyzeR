#' Blank-filter a dataset (mask values at or below the blank level)
#'
#' For each feature, computes a summary of the blank injections times
#' `blank_filter` and sets any value below that threshold to `NA`. Blank
#' injections are identified from `blank_head` in the dataset's `sample_meta`,
#' so nothing has to be passed separately.
#'
#' `summary = "median"` is worth considering over the default mean: a single
#' contaminated blank drags a mean threshold up and can mask most of a
#' compound's real measurements, and with the two or three blanks a typical
#' run carries there is no way to see that from the result alone.
#'
#' The per-feature threshold and the number of values it masked are recorded in
#' `variable_meta` as `blank_threshold` and `n_masked_blank`, so the filter's
#' effect is inspectable afterwards rather than only visible as missing data.
#'
#' @param de A `struct::DatasetExperiment`.
#' @param blank_filter Numeric fold-change multiplier applied to the blank
#'   summary (e.g. `3` keeps only peaks at least three times the blank level).
#' @param blank_head `sample_meta` column identifying blank injections.
#' @param blank_name Value in `blank_head` marking a blank injection.
#' @param summary `"mean"` (default) or `"median"` of the blank injections.
#' @param min_blanks Minimum blank injections required. Below this the
#'   threshold rests on too little to be meaningful, and `missing_blanks`
#'   decides what happens.
#' @param missing_blanks What to do when fewer than `min_blanks` are found:
#'   `"warn"` (default) filters anyway and says so, `"error"` stops, `"skip"`
#'   returns the dataset unfiltered.
#' @return `de` with sub-threshold values in `data` set to `NA`, and
#'   `blank_threshold` / `n_masked_blank` added to `variable_meta`.
#' @family peak-matrix processing
#' @examples
#' de <- struct::DatasetExperiment(
#'   data = data.frame(PGE2 = c(100, 90, 5, 4), PGD2 = c(50, 60, 1, 2),
#'                     row.names = c("S1", "S2", "B1", "B2")),
#'   sample_meta = data.frame(
#'     Sample_type = c("Sample", "Sample", "Blank", "Blank"),
#'     row.names = c("S1", "S2", "B1", "B2")),
#'   variable_meta = data.frame(Compound = c("PGE2", "PGD2"),
#'                              row.names = c("PGE2", "PGD2")))
#' out <- filterBlanks(de, blank_filter = 3)
#' out$data
#' as.data.frame(out$variable_meta)[, c("blank_threshold", "n_masked_blank")]
#' @export
filterBlanks = function(de, blank_filter, blank_head = "Sample_type",
                         blank_name = "Blank",
                         summary = c("mean", "median"),
                         min_blanks = 2,
                         missing_blanks = c("warn", "error", "skip")){

  summary        = match.arg(summary)
  missing_blanks = match.arg(missing_blanks)

  df    = as.data.frame(de$data)
  smeta = as.data.frame(de$sample_meta)
  vm    = as.data.frame(de$variable_meta)
  blank_head = .resolve_meta_col(blank_head, smeta) %||% blank_head

  if(!blank_head %in% colnames(smeta))
    stop(sprintf("[filterBlanks] '%s' is not a sample_meta column.",
                 blank_head))

  blank_samples = rownames(smeta)[smeta[[blank_head]] %in% blank_name]
  df_blank      = df[rownames(df) %in% blank_samples, , drop = FALSE]
  n_blank       = nrow(df_blank)

  if(n_blank < min_blanks){
    msg = sprintf(
      "[filterBlanks] found %d blank injection(s) (%s in {%s}); %d required.",
      n_blank, blank_head, paste(blank_name, collapse = ", "), min_blanks)
    if(missing_blanks == "error") stop(msg)
    if(missing_blanks == "skip"){
      message(msg, " Skipping the blank filter.")
      return(de)
    }
    warning(msg, " Filtering anyway - thresholds rest on very few values.")
  }

  fun = if(summary == "median")
    function(x) stats::median(x, na.rm = TRUE) else
    function(x) mean(x, na.rm = TRUE)

  blank_thresh = if(n_blank == 0)
    stats::setNames(rep(NA_real_, ncol(df)), colnames(df)) else
    vapply(df_blank, fun, numeric(1)) * blank_filter

  # Features with no usable blank measurement give NA/NaN and are left alone:
  # no threshold is better than a threshold of zero, which would mask nothing,
  # or of NA, which would mask everything.
  masked = stats::setNames(integer(ncol(df)), colnames(df))
  for(j in which(is.finite(blank_thresh))){
    col       = df[[j]]
    hit       = !is.na(col) & col < blank_thresh[j]
    masked[j] = sum(hit)
    col[hit]  = NA
    df[[j]]   = col
  }

  key = rownames(vm)
  vm$blank_threshold = unname(blank_thresh[key])
  vm$n_masked_blank  = unname(masked[key])

  de$data          = df
  de$variable_meta = vm
  de
}
