#' Batch-correct a DatasetExperiment (median-ratio)
#'
#' Scales each batch so its per-feature median matches the grand median across
#' the reference (`qc_label`) samples, and returns the corrected
#' `DatasetExperiment`.
#'
#' The reference has to be material whose true value does not differ between
#' batches, because any difference that remains is treated as technical and
#' removed. Pooled QC injections satisfy that by construction - the same
#' extract, injected repeatedly - which is why `"QC"` is the default.
#'
#' Using the biological samples as the reference is a legitimate alternative,
#' and often a more stable one because there are more of them, but it rests on
#' the study being balanced across batches: it assumes the biology is the same
#' on average in each batch. Where that holds it is fine. Where batch is
#' confounded with the study factor - all treated animals run on day 2, say -
#' the between-batch difference is partly the effect you are looking for, and
#' correcting removes it silently. Pass `check_factor` to have that assumption
#' tested rather than assumed.
#'
#' Medians count a missing value as zero, so the median of
#' `c(NA, NA, 1, 2, 3)` is 1, not 2. A compound missing from most of a batch's
#' reference injections is at or below detection there, and a median of only
#' the injections where it was found would describe that batch by its highest
#' values. When more than half are missing the median is zero and the ratio is
#' undefined. A compound whose median is zero or negative in any batch, or
#' across all the reference injections, is corrected in no batch at all:
#' scaling only some of its batches would leave them on different footings,
#' which is the problem correction exists to remove. Its values are left as
#' measured - never turned into `NA` - the failing batches are recorded in
#' `bc_failed_batches`, and each such compound is listed on screen. Negative
#' medians only arise from signed data such as blank-subtracted values. Missing
#' values themselves are not filled in: impute afterwards, if at all, with
#' [imputeMissing()].
#'
#' Every injection needs a `batch_head` value. One with an empty value - a
#' blank run outside the sequence, say - cannot be scaled against anything,
#' so the function stops and names it: give it a batch, or exclude it from
#' processing.
#'
#' `feature_col` restricts the correction to part of the panel. It names a
#' `variable_meta` column holding `YES` or `NO` per compound; case and
#' surrounding spaces are ignored, and `TRUE` / `FALSE` work too. A `NO`
#' compound is left as measured. A blank cell is treated as `NO`, with a
#' warning, because leaving a value alone is the choice that can be undone.
#'
#' @param de A `struct::DatasetExperiment`.
#' @param qc_label Label identifying reference samples in `factor_name`.
#' @param factor_name `sample_meta` column holding `qc_label`.
#' @param batch_head `sample_meta` column identifying the batch. Every
#'   injection needs a value; an empty one stops with an error naming it.
#' @param check_factor Optional `sample_meta` column holding the study factor.
#'   When given, the reference samples are cross-tabulated against it and a
#'   warning is issued if any batch carries only one level - the signature of
#'   a design where correcting on biological samples would remove real
#'   differences. `NULL` skips the check.
#' @param min_ref Minimum reference samples required in each batch. Below this
#'   the batch median is not estimable and the correction would be noise.
#' @param feature_col `NULL` (default) to correct every compound, or the name
#'   of a `variable_meta` column of `YES` / `NO` selecting the compounds to
#'   correct.
#' @return The batch-corrected `struct::DatasetExperiment`, with two columns
#'   added to `variable_meta`: `batch_corrected`, `TRUE` where the compound was
#'   corrected, and `bc_failed_batches`, which for a selected compound that
#'   could not be corrected names the batches whose reference median was zero
#'   or negative - `"grand median"` when the median across all reference
#'   injections was - and is empty otherwise.
#' @examples
#' m  <- data.frame(A = c(10, 12, 20, 24), B = c(5, 6, 10, 12),
#'                  row.names = c("S1", "S2", "S3", "S4"))
#' fm <- data.frame(Compound = c("A", "B"), Correct = c("YES", "NO"),
#'                  row.names = c("A", "B"))
#' sm <- data.frame(Sample_type = "Sample",
#'                  Chrom_Batch = c("b1", "b1", "b2", "b2"),
#'                  row.names = c("S1", "S2", "S3", "S4"))
#' de <- struct::DatasetExperiment(data = m, sample_meta = sm, variable_meta = fm)
#' correctBatch(de, qc_label = "Sample", factor_name = "Sample_type",
#'               batch_head = "Chrom_Batch")
#'
#' # Only the compounds marked YES
#' out <- correctBatch(de, qc_label = "Sample", factor_name = "Sample_type",
#'                     batch_head = "Chrom_Batch", feature_col = "Correct")
#' out$data
#' @family peak-matrix processing
#' @export
correctBatch = function(de, qc_label = "QC", factor_name = "Sample_type",
                         batch_head = "Chrom_Batch",
                         check_factor = NULL, min_ref = 2,
                         feature_col = NULL){

  smeta = as.data.frame(de$sample_meta)
  # Headings are stored in make.names() form (Chrom-Batch -> Chrom.Batch).
  factor_name = .resolve_meta_col(factor_name, smeta) %||% factor_name
  batch_head  = .resolve_meta_col(batch_head,  smeta) %||% batch_head
  if(!is.null(check_factor))
    check_factor = .resolve_meta_col(check_factor, smeta) %||% check_factor
  for(h in c(factor_name, batch_head)){
    if(!h %in% colnames(smeta))
      stop(sprintf("[correctBatch] '%s' is not a sample_meta column.", h))
  }

  # Every injection being corrected has to belong to a batch. One with an
  # empty label cannot be scaled against anything, and leaving it as measured
  # beside corrected neighbours would put it on a different footing - so the
  # run stops and names it.
  batch    = as.character(smeta[[batch_head]])
  no_batch = is.na(batch) | !nzchar(trimws(batch))
  if(any(no_batch))
    stop(sprintf(
      "[correctBatch] %d injection(s) have no '%s' value: %s%s. Give each a batch, or set Include to NO for injections that should not be processed.",
      sum(no_batch), batch_head,
      paste(utils::head(rownames(smeta)[no_batch], 5), collapse = ", "),
      if(sum(no_batch) > 5) sprintf(" and %d more", sum(no_batch) - 5) else ""),
      call. = FALSE)
  is_ref = smeta[[factor_name]] %in% qc_label

  if(!any(is_ref))
    stop(sprintf(
      "[correctBatch] no reference samples: no row has %s in {%s}.",
      factor_name, paste(qc_label, collapse = ", ")))

  # Every batch needs enough reference samples for its median to mean
  # anything; a batch with one is scaled by a single observation.
  per_batch = table(batch[is_ref])
  thin = setdiff(unique(batch), names(per_batch)[per_batch >= min_ref])
  if(length(thin))
    stop(sprintf(
      "[correctBatch] batch(es) %s have fewer than %d reference samples (%s in {%s}). Add reference samples, lower min_ref, or use a reference present in every batch.",
      paste(sprintf("'%s'", thin), collapse = ", "), min_ref,
      factor_name, paste(qc_label, collapse = ", ")))

  # Is the "biology is the same on average in each batch" assumption safe?
  if(!is.null(check_factor)){
    if(!check_factor %in% colnames(smeta)){
      warning(sprintf("[correctBatch] check_factor '%s' is not a sample_meta column; confounding check skipped.",
                      check_factor))
    } else {
      lv = table(batch[is_ref],
                 as.character(smeta[[check_factor]][is_ref]))
      n_lv = rowSums(lv > 0)
      if(ncol(lv) > 1 && any(n_lv < 2))
        warning(sprintf(
          "[correctBatch] batch(es) %s contain only one level of '%s' among the reference samples, so batch is confounded with it. Median-ratio correction will remove that difference along with the technical shift. Use pooled QCs as the reference if you have them.",
          paste(sprintf("'%s'", names(n_lv)[n_lv < 2]), collapse = ", "),
          check_factor))
    }
  }

  X  = as.data.frame(de$data)
  vm = as.data.frame(de$variable_meta)

  # --- which compounds ------------------------------------------------------
  selected = rep(TRUE, ncol(X))
  if(!is.null(feature_col)){
    col = .resolve_meta_col(feature_col, vm)
    if(is.null(col))
      stop(sprintf(
        "[correctBatch] feature_col '%s' is not a variable_meta column.",
        feature_col), call. = FALSE)
    flag = .yes_no(vm[[col]][match(colnames(X), rownames(vm))], feature_col)
    if(any(is.na(flag)))
      warning(sprintf(
        "[correctBatch] %d compound(s) have no YES/NO in '%s' and were left uncorrected: %s",
        sum(is.na(flag)), feature_col,
        paste(utils::head(colnames(X)[is.na(flag)], 5), collapse = ", ")),
        call. = FALSE)
    selected = flag %in% TRUE
  }

  # --- medians --------------------------------------------------------------
  # A missing value counts as zero in every median - see Details.
  med0 = function(rows)
    vapply(X[rows, , drop = FALSE], function(v){
      v[is.na(v)] = 0
      stats::median(v)
    }, numeric(1))

  batches = unique(batch)
  grand   = med0(is_ref)
  # compounds x batches
  bmed = matrix(vapply(batches, function(b) med0(batch == b & is_ref),
                       numeric(ncol(X))),
                nrow = ncol(X), dimnames = list(colnames(X), batches))

  # A zero or negative median gives no usable ratio. Correcting a compound in
  # some batches but not others would leave its batches on different
  # footings - the problem this step exists to remove - so one failing batch
  # means the compound is corrected in none.
  bad_grand = !(is.finite(grand) & grand > 0)
  bad_batch = !(is.finite(bmed) & bmed > 0)
  failed    = selected & (bad_grand | rowSums(bad_batch) > 0)
  correct   = selected & !failed

  # --- correction -----------------------------------------------------------
  for(b in batches){
    in_b = which(batch == b)
    for(j in which(correct))
      X[[j]][in_b] = X[[j]][in_b] / (bmed[j, b] / grand[j])
  }

  cause = vapply(seq_len(ncol(X)), function(j){
    if(!failed[j]) return("")
    paste(c(if(bad_grand[j]) "grand median", batches[bad_batch[j, ]]),
          collapse = ", ")
  }, character(1))

  if(any(failed)){
    # On screen as it happens, one line per compound; the warning repeats the
    # count at the end of the run, where R collects warnings.
    message("[correctBatch] ", sum(failed), " compound(s) not corrected in ",
            "any batch - reference median zero or negative (missing values ",
            "count as zero) in:\n",
            paste(sprintf("  %s: %s", colnames(X)[failed], cause[failed]),
                  collapse = "\n"))
    warning(sprintf(
      "[correctBatch] %d compound(s) not batch-corrected because a reference median is zero or negative; see bc_failed_batches: %s",
      sum(failed),
      paste(utils::head(colnames(X)[failed], 5), collapse = ", ")),
      call. = FALSE)
  }

  key = match(rownames(vm), colnames(X))
  vm$batch_corrected   = correct[key]
  vm$bc_failed_batches = cause[key]

  de$data          = X
  de$variable_meta = vm
  de
}

#' Read a YES / NO flag column
#'
#' @param x Values from the column.
#' @param col Column name, for the error message.
#' @return Logical: `TRUE` for yes, `FALSE` for no, `NA` for a blank cell.
#'   Any other value stops, naming it.
#' @keywords internal
#' @noRd
.yes_no = function(x, col){
  if(is.logical(x)) return(x)
  v   = tolower(trimws(as.character(x)))
  out = rep(NA, length(v))
  out[v %in% c("yes", "true")] = TRUE
  out[v %in% c("no", "false")] = FALSE
  bad = !is.na(v) & nzchar(v) & is.na(out)
  if(any(bad))
    stop(sprintf("[batch correction] '%s' must hold YES or NO; found: %s.",
                 col, paste(sprintf("'%s'", unique(x[bad])), collapse = ", ")),
         call. = FALSE)
  out
}
