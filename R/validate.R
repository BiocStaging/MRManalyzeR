#' Collect validation findings
#'
#' A tiny accumulator so each check reports itself the same way and the caller
#' gets one consolidated answer instead of the first failure.
#'
#' @return An environment with `errors`, `warnings` and `rows`.
#' @keywords internal
#' @noRd
.chk_new = function(){
  e = new.env(parent = emptyenv())
  e$errors = character(); e$warnings = character(); e$rows = list()
  e
}

#' @param env Accumulator from [.chk_new()].
#' @param check Short check name.
#' @param status `"ok"`, `"warning"` or `"error"`.
#' @param message Explanation; `NA` for a passing check.
#' @keywords internal
#' @noRd
.chk_add = function(env, check, status, message = NA_character_){
  env$rows[[length(env$rows) + 1L]] = data.frame(
    check = check, status = status, message = message,
    stringsAsFactors = FALSE)
  if(identical(status, "error"))   env$errors   = c(env$errors, message)
  if(identical(status, "warning")) env$warnings = c(env$warnings, message)
  invisible(env)
}

#' @keywords internal
#' @noRd
.chk_out = function(env){
  structure(list(errors   = env$errors,
                 warnings = env$warnings,
                 summary  = do.call(rbind, env$rows)),
            class = "mrm_validation")
}

#' Print a validation result
#'
#' @param x An `mrm_validation` object.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.mrm_validation = function(x, ...){
  n_e = length(x$errors); n_w = length(x$warnings)
  cat(sprintf("MRManalyzeR validation: %d error(s), %d warning(s), %d check(s)\n",
              n_e, n_w, nrow(x$summary)))
  if(n_e) cat("\nErrors:\n",   paste0("  - ", x$errors,   collapse = "\n"), "\n")
  if(n_w) cat("\nWarnings:\n", paste0("  - ", x$warnings, collapse = "\n"), "\n")
  if(!n_e && !n_w) cat("All checks passed.\n")
  invisible(x)
}

#' Validate a workbook before processing it
#'
#' Checks the shape of the inputs - the things that make a run fail halfway
#' through, or worse, succeed on the wrong data. Every check runs, so one call
#' reports every problem rather than stopping at the first.
#'
#' @param feature_meta Feature metadata sheet.
#' @param sample_meta Sample metadata sheet.
#' @param data Optional wide matrix (samples x compounds), e.g. from
#'   [readTargetLynx()]. Cross-checks between data and metadata are skipped
#'   when absent.
#' @param name_col,compound_col,processing_name_col,report_col,include_col
#'   Column names, matching the arguments of [processDataset()].
#' @param include_value Value of `include_col` marking injections to process.
#' @param injection_order_head `sample_meta` column holding acquisition order,
#'   or `NULL` (default) to skip the check. When given, every included
#'   injection must have a numeric value there: drift plots and the QC-PCA
#'   colour scale read it as a number, and a text value (a formula returning
#'   `"_48"`, say) would otherwise leave them empty without an error. An absent
#'   column is only a warning, since the plots then fall back to row order.
#' @return An `mrm_validation` object with `errors`, `warnings` and a `summary`
#'   data frame of every check.
#' @examples
#' xlsx <- system.file("extdata", "example_data.xlsx",
#'                     package = "MRManalyzeR")
#' fdata <- openxlsx::read.xlsx(xlsx, sheet = "feature_metadata")
#' meta  <- openxlsx::read.xlsx(xlsx, sheet = "sample_metadata")
#' validateInput(fdata, meta)
#' @family data parse
#' @export
validateInput = function(feature_meta, sample_meta, data = NULL,
                          name_col            = "Name",
                          compound_col        = "Compound",
                          processing_name_col = "Processing_name",
                          report_col          = "Report",
                          include_col         = "Include",
                          include_value       = "YES",
                          injection_order_head = NULL){

  chk = .chk_new()
  fm  = as.data.frame(feature_meta)
  sm  = as.data.frame(sample_meta)

  need = function(df, cols, what, label){
    miss = setdiff(cols, colnames(df))
    if(length(miss))
      .chk_add(chk, label, "error", sprintf(
        "%s is missing required column(s): %s. Present: %s.",
        what, paste(miss, collapse = ", "),
        paste(colnames(df), collapse = ", ")))
    else .chk_add(chk, label, "ok")
    length(miss) == 0
  }

  ok_f = need(fm, c(compound_col, processing_name_col),
              "feature_metadata", "feature_metadata columns")
  ok_s = need(sm, c(name_col, include_col), "sample_metadata",
              "sample_metadata columns")

  # The report filter is optional: an externally curated feature table carries
  # no exclusion column, and every feature is then reported. Worth a warning
  # rather than silence, because the other way this happens is a typo in
  # report_col, and that would quietly reinstate features meant to be dropped.
  has_report = report_col %in% colnames(fm)
  if(!has_report)
    .chk_add(chk, "feature report filter", "warning", sprintf(
      "feature_metadata has no '%s' column, so all %d features will be reported. Set report_col to an existing column if some should be excluded.",
      report_col, nrow(fm)))
  else .chk_add(chk, "feature report filter", "ok")

  # --- identifiers must be unique ------------------------------------------
  if(ok_s){
    ids = as.character(sm[[name_col]])
    dup = unique(ids[duplicated(ids)])
    if(length(dup))
      .chk_add(chk, "unique sample names", "error", sprintf(
        "sample_metadata$%s has %d duplicate value(s): %s. Every injection needs its own identifier.",
        name_col, length(dup), paste(utils::head(dup, 5), collapse = ", ")))
    else .chk_add(chk, "unique sample names", "ok")

    blank_id = sum(is.na(ids) | !nzchar(trimws(ids)))
    if(blank_id)
      .chk_add(chk, "sample names present", "error", sprintf(
        "sample_metadata$%s has %d empty value(s).", name_col, blank_id))
    else .chk_add(chk, "sample names present", "ok")
  }

  if(ok_f){
    reported = if(has_report) as.character(fm[[report_col]]) == "YES"
               else           rep(TRUE, nrow(fm))
    for(cl in unique(c(compound_col, processing_name_col))){
      v   = as.character(fm[[cl]])[reported]
      dup = unique(v[duplicated(v)])
      if(length(dup))
        .chk_add(chk, sprintf("unique %s", cl), "error", sprintf(
          "feature_metadata$%s has %d duplicate value(s) among reported compounds: %s.",
          cl, length(dup), paste(utils::head(dup, 5), collapse = ", ")))
      else .chk_add(chk, sprintf("unique %s", cl), "ok")
    }
  }

  # --- injection order ------------------------------------------------------
  if(ok_s && !is.null(injection_order_head)){
    if(!injection_order_head %in% colnames(sm)){
      .chk_add(chk, "injection order numeric", "warning", sprintf(
        "sample_metadata has no '%s' column; drift plots will use row order.",
        injection_order_head))
    } else {
      incl = which(as.character(sm[[include_col]]) == include_value)
      raw  = sm[[injection_order_head]][incl]
      txt  = trimws(as.character(raw))
      bad  = is.na(raw) | !nzchar(txt) |
             is.na(suppressWarnings(as.numeric(txt)))
      if(any(bad))
        .chk_add(chk, "injection order numeric", "error", sprintf(
          "%d included injection(s) have a missing or non-numeric %s, e.g. %s. Injection order must be a number.",
          sum(bad), injection_order_head,
          paste(sprintf("'%s'", utils::head(unique(txt[bad]), 5)),
                collapse = ", ")))
      else .chk_add(chk, "injection order numeric", "ok")
    }
  }

  # --- data / metadata agreement -------------------------------------------
  if(!is.null(data) && ok_s){
    dm  = as.data.frame(data)
    ids = as.character(sm[[name_col]])

    orphan = setdiff(rownames(dm), ids)
    if(length(orphan))
      .chk_add(chk, "samples in data have metadata", "error", sprintf(
        "%d sample(s) in the data have no sample_metadata row: %s. Usually a trailing space or a renamed injection.",
        length(orphan), paste(utils::head(orphan, 5), collapse = ", ")))
    else .chk_add(chk, "samples in data have metadata", "ok")

    incl = ids[as.character(sm[[include_col]]) == include_value]
    absent = setdiff(incl, rownames(dm))
    if(length(absent))
      .chk_add(chk, "included samples are in the data", "warning", sprintf(
        "%d sample(s) marked %s = %s have no rows in the data: %s.",
        length(absent), include_col, include_value,
        paste(utils::head(absent, 5), collapse = ", ")))
    else .chk_add(chk, "included samples are in the data", "ok")

    non_num = names(dm)[!vapply(dm, is.numeric, logical(1))]
    if(length(non_num))
      .chk_add(chk, "measurements are numeric", "error", sprintf(
        "%d measurement column(s) are not numeric: %s. A stray text entry (e.g. 'n.d.') in one cell turns the whole compound into text.",
        length(non_num), paste(utils::head(non_num, 5), collapse = ", ")))
    else .chk_add(chk, "measurements are numeric", "ok")

    empty = names(dm)[colSums(!is.na(dm)) == 0]
    if(length(empty))
      .chk_add(chk, "no all-missing compounds", "warning", sprintf(
        "%d compound(s) have no values at all and will be dropped: %s.",
        length(empty), paste(utils::head(empty, 5), collapse = ", ")))
    else .chk_add(chk, "no all-missing compounds", "ok")
  }

  .chk_out(chk)
}

#' Validate a YAML configuration
#'
#' Structural checks on the config alone - blocks present, exactly one data
#' source enabled, enumerated values spelled correctly. This catches the
#' mistakes that would otherwise surface as an obscure error deep in report
#' rendering.
#'
#' @param config A config list, or a path to a YAML file.
#' @return An `mrm_validation` object.
#' @examples
#' validateConfig(system.file("extdata", "example_config.yml",
#'                             package = "MRManalyzeR"))
#' @family data parse
#' @export
validateConfig = function(config){

  if(is.character(config) && length(config) == 1L) config = loadConfig(config)
  chk = .chk_new()
  prj = config$project %||% config

  pmp = prj$PeakMatrixProcessing
  if(is.null(pmp))
    .chk_add(chk, "PeakMatrixProcessing block", "error",
             "config has no 'PeakMatrixProcessing:' block.")
  else .chk_add(chk, "PeakMatrixProcessing block", "ok")

  if(!is.null(pmp)){
    # Exactly one source, or the run silently reads the wrong sheet.
    on = c(skyline_data = isTRUE(pmp$skyline_data$enabled),
           tl_data       = isTRUE(pmp$tl_data$enabled),
           matrix_data   = isTRUE(pmp$matrix_data$enabled))
    if(sum(on) != 1)
      .chk_add(chk, "one data source enabled", "error", sprintf(
        "%d data sources enabled (%s); exactly one must be.",
        sum(on), if(any(on)) paste(names(on)[on], collapse = ", ")
                 else "none"))
    else .chk_add(chk, "one data source enabled", "ok")

    if(isTRUE(pmp$adjust_conc$enabled)){
      vols = c("sample_vol_col", "cal_vol_col", "sample_IS_col", "cal_IS_col")
      miss = vols[vapply(vols, function(k) is.null(pmp$adjust_conc[[k]]),
                         logical(1))]
      if(length(miss))
        .chk_add(chk, "concentration adjustment columns", "error", sprintf(
          "adjust_conc is enabled but %s not set.",
          paste(miss, collapse = ", ")))
      else .chk_add(chk, "concentration adjustment columns", "ok")
    }
  }

  st = prj$stats_report %||% prj$MVA_report
  # Only a report that will run, which means execute: True - the same
  # condition runMRManalyzeR() renders it under. An invalid method in a report
  # that is switched off must not stop the processing run.
  if(!is.null(st) && isTRUE(st$execute)){
    for(e in .section_entries(st$comparisons, "comparisons")){
      nm = e$name %||% "<unnamed>"
      lv = length(e$compare$levels %||% character(0))
      m  = tolower(as.character(e$method %||% "")[1])
      allowed = if(lv == 2) c("", "welch", "student", "wilcoxon", "t_test")
                else        c("", "anova", "kruskal")
      if(!m %in% allowed)
        .chk_add(chk, sprintf("comparison '%s' method", nm), "error", sprintf(
          "method '%s' is not valid for %d level(s). Use one of: %s.",
          m, lv, paste(setdiff(allowed, ""), collapse = ", ")))
      else if(!nzchar(m))
        .chk_add(chk, sprintf("comparison '%s' method", nm), "warning",
                 sprintf("comparison '%s' has no method:; the default will be used. Pre-specify the test.",
                         nm))
      else .chk_add(chk, sprintf("comparison '%s' method", nm), "ok")

      if(isTRUE(e$paired) && is.null(e$subject))
        .chk_add(chk, sprintf("comparison '%s' pairing", nm), "error",
                 sprintf("comparison '%s' is paired: true but has no subject:.",
                         nm))
    }
  }

  .chk_out(chk)
}

#' Validate a study design against the data
#'
#' The checks that need the data and the configuration together: are the
#' reference samples actually present, does every batch have them, is the
#' design balanced enough for the correction and the tests requested. These are
#' the failures that produce a plausible-looking result rather than an error,
#' which is what makes them worth checking up front.
#'
#' @param de A `struct::DatasetExperiment`.
#' @param config Optional config list or path; when supplied, the requested
#'   comparisons and normalisation are checked against the data.
#' @param sample_type_head,qc_label,blank_name Sample-type column and the
#'   values marking QC and blank injections.
#' @param batch_head Batch column.
#' @param min_group Minimum observations per group before a comparison is
#'   worth running.
#' @return An `mrm_validation` object.
#' @examples
#' de <- loadDataset(system.file("extdata", "example_synthetic.RDS",
#'                                package = "MRManalyzeR"))
#' validateDesign(de)
#' @family data parse
#' @export
validateDesign = function(de, config = NULL,
                           sample_type_head = "Sample_type",
                           qc_label   = "QC",
                           blank_name = "Blank",
                           batch_head = "Chrom_Batch",
                           min_group  = 3){

  chk   = .chk_new()
  smeta = as.data.frame(de$sample_meta)
  dm    = as.data.frame(de$data)
  # The config may name the batch heading as written in the workbook
  # (Chrom-Batch); the dataset stores it as Chrom.Batch.
  batch_head = .resolve_meta_col(batch_head, smeta) %||% batch_head

  if(!sample_type_head %in% colnames(smeta)){
    .chk_add(chk, "sample-type column", "error", sprintf(
      "'%s' is not a sample_meta column.", sample_type_head))
    return(.chk_out(chk))
  }
  .chk_add(chk, "sample-type column", "ok")

  st = as.character(smeta[[sample_type_head]])

  for(lab in c(qc_label, blank_name)){
    n = sum(st %in% lab)
    if(n == 0)
      .chk_add(chk, sprintf("%s injections present", lab), "warning", sprintf(
        "no injections with %s = '%s'; the checks that depend on them cannot run.",
        sample_type_head, lab))
    else .chk_add(chk, sprintf("%s injections present", lab), "ok")
  }

  # --- batch structure ------------------------------------------------------
  if(batch_head %in% colnames(smeta)){
    batch = as.character(smeta[[batch_head]])
    .chk_add(chk, "batch column", "ok")

    if(any(st %in% qc_label)){
      per = table(batch[st %in% qc_label])
      # An injection with no batch label sits batch correction out, so it is
      # not a batch that needs QCs.
      labelled = batch[!is.na(batch) & nzchar(trimws(batch))]
      thin = setdiff(unique(labelled), names(per)[per >= 2])
      if(length(thin))
        .chk_add(chk, "QCs in every batch", "warning", sprintf(
          "batch(es) %s have fewer than 2 QC injections, so QC-referenced batch correction is not available for them.",
          paste(sprintf("'%s'", thin), collapse = ", ")))
      else .chk_add(chk, "QCs in every batch", "ok")
    }
  } else {
    .chk_add(chk, "batch column", "warning", sprintf(
      "'%s' is not a sample_meta column; batch checks skipped.", batch_head))
    batch = NULL
  }

  # --- confounding, comparisons, normalisation ------------------------------
  if(!is.null(config)){
    if(is.character(config) && length(config) == 1L) config = loadConfig(config)
    prj = config$project %||% config
    pmp = prj$PeakMatrixProcessing
    sp  = prj$stats_report %||% prj$MVA_report

    norm = pmp$normalize
    if(!isFALSE(norm) && !is.null(norm)){
      if(!norm %in% colnames(smeta)){
        .chk_add(chk, "normalisation column", "error", sprintf(
          "normalize = '%s' is not a sample_meta column.", norm))
      } else {
        v = suppressWarnings(as.numeric(smeta[[norm]]))
        bad = sum(is.na(v) | v == 0)
        if(bad)
          .chk_add(chk, "normalisation divisors", "error", sprintf(
            "%d sample(s) have a missing or zero '%s'; dividing by it yields NA or Inf.",
            bad, norm))
        else .chk_add(chk, "normalisation divisors", "ok")
      }
    }

    # Only the comparisons of a stats report that will run: anything else
    # never tests them, so warning about their levels is noise.
    comps = list()
    if(isTRUE(sp$execute))
      comps = .section_entries(sp$comparisons, "comparisons")
    for(e in comps){
      nm  = e$name %||% "<unnamed>"
      fct = e$compare$factor
      lv  = e$compare$levels

      if(is.null(fct) || !fct %in% colnames(smeta)){
        .chk_add(chk, sprintf("comparison '%s' factor", nm), "error", sprintf(
          "factor '%s' is not a sample_meta column.", fct %||% "<NULL>"))
        next
      }

      have = as.character(smeta[[fct]])
      miss = setdiff(lv, unique(have))
      if(length(miss)){
        .chk_add(chk, sprintf("comparison '%s' levels", nm), "error", sprintf(
          "level(s) %s not present in '%s'. Present: %s.",
          paste(sprintf("'%s'", miss), collapse = ", "), fct,
          paste(utils::head(unique(have), 8), collapse = ", ")))
        next
      }

      n_per = table(have[have %in% lv])
      if(any(n_per < min_group))
        .chk_add(chk, sprintf("comparison '%s' group sizes", nm), "warning",
                 sprintf("group size(s) %s below %d; the test will run but is not worth much.",
                         paste(sprintf("%s=%d", names(n_per), n_per),
                               collapse = ", "), min_group))
      else .chk_add(chk, sprintf("comparison '%s' group sizes", nm), "ok")

      # Complete confounding: batch determines group, so nothing can separate
      # a batch effect from the effect being tested.
      if(!is.null(batch)){
        keep = have %in% lv
        tab  = table(batch[keep], have[keep])
        if(nrow(tab) > 1 && all(rowSums(tab > 0) == 1))
          .chk_add(chk, sprintf("comparison '%s' vs batch", nm), "error",
                   sprintf("'%s' is completely confounded with '%s' - each batch holds only one level, so a batch effect and the treatment effect cannot be told apart.",
                           fct, batch_head))
        else .chk_add(chk, sprintf("comparison '%s' vs batch", nm), "ok")
      }
    }
  }

  # --- features -------------------------------------------------------------
  empty = names(dm)[colSums(!is.na(dm)) == 0]
  if(length(empty))
    .chk_add(chk, "no all-missing compounds", "warning", sprintf(
      "%d compound(s) have no values: %s.", length(empty),
      paste(utils::head(empty, 5), collapse = ", ")))
  else .chk_add(chk, "no all-missing compounds", "ok")

  .chk_out(chk)
}
