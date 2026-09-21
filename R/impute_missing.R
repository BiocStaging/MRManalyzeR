#' Impute missing values below each compound's lowest measurement
#'
#' Zeros are treated as missing. For each feature, the `NA`s among the
#' non-blank injections are filled with values below the lowest one actually
#' measured, because a non-detect sits under the detection floor rather than
#' at zero. The ceiling is `scalar` times that minimum.
#'
#' With `method = "gaussian"` (the default) every missing cell gets its own
#' random draw from a normal distribution centred at half the ceiling, with a
#' standard deviation of a sixth of it, truncated to lie between 0 and the
#' ceiling. A single constant would give every non-detect the same value,
#' which shrinks the variance of a compound and ties its ranks - both of which
#' flatter a statistical test. `method = "constant"` fills every cell with the
#' ceiling itself, as earlier versions did. The draws use `seed`, so a rerun
#' gives the same values; `seed = NULL` uses the session's random stream.
#'
#' `by` names a `sample_meta` column of processing batches. Each batch then
#' gets its own ceiling from its own minimum, so a batch measured at lower
#' response is not filled with another batch's values. Injections with an
#' empty `by` value are treated as one batch of their own.
#'
#' Imputation is the one step that invents numbers, so it records what it did
#' in `variable_meta`: `impute_max`, the ceiling (one column per batch,
#' `impute_max_<batch>`, when `by` gives more than one); `n_imputed`, how many
#' cells were filled; and `frac_imputed`, the share of the compound's
#' non-blank injections that were filled. A compound whose `frac_imputed` is
#' high is mostly invented, and a difference found in it is an artefact of the
#' fill. Check that column before interpreting anything, and consider
#' excluding compounds above `warn_frac`.
#'
#' Blank injections (identified from `blank_head` in `sample_meta`) are left
#' untouched and do not count towards the minimum.
#'
#' @param de A `struct::DatasetExperiment`.
#' @param scalar Fraction of the per-feature minimum used as the ceiling,
#'   e.g. `0.2`.
#' @param blank_head `sample_meta` column identifying blank injections.
#' @param blank_name Value in `blank_head` marking a blank injection.
#' @param warn_frac Warn about features imputed in more than this fraction of
#'   the non-blank injections. `NULL` disables the warning.
#' @param method `"gaussian"` (default) for random draws below the ceiling, or
#'   `"constant"` to fill every cell with the ceiling.
#' @param seed Random seed for `method = "gaussian"`; `NULL` uses the current
#'   random stream.
#' @param by Optional `sample_meta` column of processing batches, each given
#'   its own ceiling. `NULL` (default) treats all injections together.
#' @return `de` with `data` imputed in the original row order, and
#'   `impute_max` (or `impute_max_<batch>`), `n_imputed` and `frac_imputed`
#'   added to `variable_meta`.
#' @family peak-matrix processing
#' @examples
#' de <- struct::DatasetExperiment(
#'   data = data.frame(PGE2 = c(10, 5, NA, 8), PGD2 = c(2, NA, 4, NA),
#'                     row.names = c("S1", "S2", "S3", "S4")),
#'   sample_meta = data.frame(Sample_type = rep("Sample", 4),
#'                            row.names = c("S1", "S2", "S3", "S4")),
#'   variable_meta = data.frame(Compound = c("PGE2", "PGD2"),
#'                              row.names = c("PGE2", "PGD2")))
#' imputeMissing(de, scalar = 0.5)$data
#' imputeMissing(de, scalar = 0.5, method = "constant")$data
#' @export
imputeMissing = function(de, scalar, blank_head = "Sample_type",
                          blank_name = "Blank", warn_frac = 0.5,
                          method = c("gaussian", "constant"),
                          seed = 42, by = NULL){

  method = match.arg(method)
  if(isTRUE(scalar))
    stop("`scalar` must be numeric (e.g. 0.2), not TRUE.")

  df    = as.data.frame(de$data)
  smeta = as.data.frame(de$sample_meta)
  vm    = as.data.frame(de$variable_meta)
  # A heading typed as in the workbook (Sample-type) is stored as Sample.type.
  blank_head = .resolve_meta_col(blank_head, smeta) %||% blank_head
  if(!is.null(by)) by = .resolve_meta_col(by, smeta) %||% by

  blank_samples = if(blank_head %in% colnames(smeta))
    rownames(smeta)[smeta[[blank_head]] %in% blank_name] else character(0)

  df[df == 0] = NA
  is_blank = rownames(df) %in% blank_samples

  # Processing batches: each gets its own ceiling from its own minimum.
  grp = rep("_all_", nrow(df))
  if(!is.null(by)){
    if(!by %in% colnames(smeta))
      stop(sprintf("[imputeMissing] by = '%s' is not a sample_meta column.",
                   by), call. = FALSE)
    grp = as.character(smeta[[by]])
    grp[is.na(grp) | !nzchar(trimws(grp))] = "_unbatched_"
  }
  groups = unique(grp)

  # Ceiling per feature (rows) and batch (columns): scalar x the lowest
  # non-blank value measured in that batch.
  ceil = vapply(groups, function(g)
    suppressWarnings(vapply(df[grp == g & !is_blank, , drop = FALSE], min,
                            numeric(1), na.rm = TRUE)) * scalar,
    numeric(ncol(df)))
  ceil = matrix(ceil, nrow = ncol(df), dimnames = list(colnames(df), groups))

  fill = function(){
    n_imp = stats::setNames(integer(ncol(df)), colnames(df))
    for(g in groups){
      rows = grp == g & !is_blank
      for(j in seq_len(ncol(df))){
        mx  = ceil[j, g]
        hit = rows & is.na(df[[j]])
        k   = sum(hit)
        if(!k || !is.finite(mx)) next
        # A ceiling at or below zero (signed data) leaves no range to draw
        # from, so those cells get the ceiling itself.
        df[[j]][hit] = if(identical(method, "gaussian") && mx > 0)
          .draw_below(k, mx) else mx
        n_imp[j] = n_imp[j] + k
      }
    }
    list(df = df, n_imp = n_imp)
  }
  res = if(identical(method, "gaussian") && !is.null(seed))
    withr::with_seed(seed, fill()) else fill()
  df    = res$df
  n_imp = res$n_imp

  n_sample = sum(!is_blank)
  frac_imp = if(n_sample > 0) n_imp / n_sample else n_imp * NA_real_

  key = rownames(vm)
  at  = match(key, rownames(ceil))
  ceil_col = function(g){ v = unname(ceil[at, g]); v[!is.finite(v)] = NA; v }
  if(length(groups) == 1L){
    vm$impute_max = ceil_col(groups)
  } else {
    for(g in groups) vm[[paste0("impute_max_", g)]] = ceil_col(g)
  }
  vm$n_imputed    = unname(n_imp[key])
  vm$frac_imputed = round(unname(frac_imp[key]), 3)

  if(!is.null(warn_frac)){
    bad = names(frac_imp)[is.finite(frac_imp) & frac_imp > warn_frac]
    if(length(bad))
      warning(sprintf(
        "[imputeMissing] %d feature(s) imputed in more than %.0f%% of samples and are mostly invented: %s",
        length(bad), warn_frac * 100,
        paste(utils::head(bad, 5), collapse = ", ")))
  }

  de$data          = df
  de$variable_meta = vm
  de
}

#' Random values between 0 and a ceiling, clustered mid-range
#'
#' A normal distribution centred at `mx / 2` with standard deviation `mx / 6`,
#' truncated to (0, mx) by inverse-CDF sampling - the same as redrawing any
#' value outside the range, without a loop. At that spread the bounds sit at
#' +/- 3 SD, so the truncation trims only the outer 0.3 percent.
#' @param n Number of values.
#' @param mx The ceiling; must be positive.
#' @return A numeric vector of length `n`, every value in (0, mx).
#' @keywords internal
#' @noRd
.draw_below = function(n, mx){
  mu = mx / 2
  s  = mx / 6
  lo = stats::pnorm(0,  mu, s)
  hi = stats::pnorm(mx, mu, s)
  stats::qnorm(stats::runif(n, lo, hi), mu, s)
}

#' Fill each column's NAs below its minimum, for runPCA()
#'
#' Zeros count as missing, as in [imputeMissing()]: a zero is a non-detect,
#' and taking it as the minimum would put the ceiling at 0.
#' @param X Numeric matrix.
#' @param fac Fraction of the per-column minimum used as the ceiling.
#' @param random `TRUE` for draws from [.draw_below()], `FALSE` for the
#'   ceiling itself.
#' @return `X` with its NAs filled wherever the column has a finite minimum.
#' @keywords internal
#' @noRd
.fill_below_min = function(X, fac, random){
  for(j in seq_len(ncol(X))){
    x = X[, j]
    x[!is.na(x) & x == 0] = NA
    miss = is.na(x)
    if(!any(miss)) next
    mx = suppressWarnings(min(x, na.rm = TRUE)) * fac
    if(!is.finite(mx)) next
    x[miss] = if(random && mx > 0) .draw_below(sum(miss), mx) else mx
    X[, j] = x
  }
  X
}

#' How many values in a dataset were filled in by imputation
#'
#' Read from the `n_imputed` record [imputeMissing()] leaves in
#' `variable_meta`, so the reports can say what they were actually given -
#' whatever the configuration asked for.
#' @keywords internal
#' @noRd
.n_imputed_total = function(de){
  vm = as.data.frame(de$variable_meta)
  if(!"n_imputed" %in% colnames(vm)) return(0L)
  sum(vm$n_imputed, na.rm = TRUE)
}

#' Whether a stored dataset holds imputed values
#'
#' Read from the `n_imputed` record. A dataset stored before 0.99.4 by a run
#' with more than one processing batch has no record, and is then taken as
#' imputed when the configuration sets `replace_MVs`.
#' @param de A `DatasetExperiment`.
#' @param replace_MVs The configured `replace_MVs`.
#' @keywords internal
#' @noRd
.values_imputed = function(de, replace_MVs = FALSE){
  vm = colnames(as.data.frame(de$variable_meta))
  if("n_imputed" %in% vm) return(.n_imputed_total(de) > 0)
  !isFALSE(replace_MVs %||% FALSE)
}

#' Whether any combine input holds imputed values, from its `n_imputed` record
#' @keywords internal
#' @noRd
.inputs_imputed = function(paths, sample_id_col = "Sample_ID"){
  any(vapply(paths, function(p){
    de = suppressWarnings(suppressMessages(
      loadDataset(p, sample_id_col = sample_id_col)))
    .n_imputed_total(de) > 0
  }, logical(1)))
}
