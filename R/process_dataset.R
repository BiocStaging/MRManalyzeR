#' Build a sample x compound matrix from TargetLynx or Skyline input
#'
#' Applies (optionally) signal filtering (SNR for TargetLynx; LOD/LOQ for
#' Skyline), blank and detection filtering, batch correction, normalisation,
#' concentration adjustment and missing-value imputation to produce a
#' `struct::DatasetExperiment`.
#'
#' The peak-matrix workflow behind [runMRManalyzeR()], composed of the public
#' step functions, in this order: it reads the workbook ([readTargetLynx()] /
#' [readSkyline()]) into a wide sample x compound matrix and assembles it into
#' a `struct::DatasetExperiment` ([assembleDataset()]); applies
#' [filterBlanks()] per acquisition batch, then [filterFeatures()] and
#' [filterSamples()] to the whole dataset; optionally batch-corrects
#' ([correctBatch()]); applies [normaliseMatrix()] and
#' [adjustConcentration()]; and imputes last ([imputeMissing()], per batch).
#' Batch correction works on the instrument signal, before any per-sample
#' rescaling, and imputation comes last so every earlier step sees only
#' measured values.
#'
#' A value of 0 is read as a non-detect and set to `NA` straight after
#' reading, so no step counts it as a measurement: the detection filters count
#' it as missing, the blank filter leaves it out of the blank mean, and a
#' matrix that is not imputed shows a gap.
#'
#' @param fdata Feature metadata. Must contain the columns named by
#'   `compound_col`, `processing_name_col`, `report_col` (defaults
#'   `Compound` / `Processing_name` / `Report`). For `signal_filter = "LOD"`
#'   or `"LOQ"`, must also contain the matching `LOD`/`LOQ` column.
#' @param metadata Sample metadata. Must contain the columns named by
#'   `name_col`, `include_col` and the headers referenced by `bc_header`,
#'   `bc_factor_name`, `blank_head`.
#' @param xlsx_path Path to the source workbook (read by [readTargetLynx()]
#'   or [readSkyline()]).
#' @param data_source `"targetlynx"` (default), `"skyline"` or `"matrix"` -
#'   selects the reader. `"matrix"` is the vendor-neutral route via
#'   [readPeakMatrix()]: any software that can export a rectangular
#'   sample x analyte table can be used without a dedicated parser.
#' @param data_tab_names Sheet name(s) to read from `xlsx_path`. `NULL`
#'   (default) uses each reader's usual sheets: every sheet whose name
#'   contains `"lcms_data"` for `"targetlynx"`, `"skyline_data"` for
#'   `"skyline"` and `"matrix_data"` for `"matrix"`. Multiple Skyline sheets
#'   are treated as separate acquisition batches and row-bound together - see
#'   `.read_skyline_matrix()`. Sheet membership does not itself imply
#'   batch-correction grouping; that remains driven by `sample_metadata`.
#' @param datatype TargetLynx column to report (e.g. "Area", "Response",
#'   "ng/mL"). Only used when `data_source = "targetlynx"`.
#' @param tl_headers TargetLynx column headers to extract (passed to
#'   [readTargetLynx()]).
#' @param signal_filter `"SNR"`, `"LOD"`, `"LOQ"`, or `FALSE`. `"SNR"` is
#'   only valid for `data_source = "targetlynx"`; `"LOD"`/`"LOQ"` are only
#'   valid for `data_source = "skyline"` or `"matrix"`. Default `NULL` infers
#'   `"SNR"`
#'   unless `snr = FALSE`, for backwards compatibility with configs that
#'   only set `snr:`.
#' @param snr S/N threshold. Only consulted when `signal_filter = "SNR"`.
#' @param blank_filter Blank-filter factor; `FALSE` or `NULL` to skip.
#' @param replace_MVs Scalar passed to [`imputeMissing()`]; `FALSE` or `NULL`
#'   to skip. Applied last, after batch correction, normalisation and
#'   concentration adjustment, per `processing_batch`.
#' @param impute_method,impute_seed Passed to [imputeMissing()] as `method`
#'   and `seed`: `"gaussian"` (default) draws each fill at random below the
#'   per-feature ceiling, `"constant"` uses the ceiling itself.
#' @param batch_correction `FALSE` (default) to skip, `TRUE` to correct every
#'   feature, or the name of an `fdata` column holding `YES` / `NO` per
#'   feature, so that only the `YES` features are corrected. Runs after
#'   the blank and detection filters, on the instrument signal, and before
#'   normalisation, concentration adjustment and imputation - see
#'   [correctBatch()].
#' @param bc_qc_label,bc_factor_name,bc_header Batch-correction parameters:
#'   the reference-sample label, the column holding it, and the batch column.
#'   Defaults to pooled QCs, which are the only reference guaranteed not to
#'   differ between batches for biological reasons - see [correctBatch()].
#' @param bc_check_factor Optional study factor; when given, `correctBatch()`
#'   warns if batch is confounded with it.
#' @param matrix_id_col,matrix_orientation Passed to [readPeakMatrix()] when
#'   `data_source = "matrix"`: the sample-identifier column (`NULL` = first)
#'   and whether samples are rows or columns.
#' @param processing_batch `sample_meta` column used to partition the *per-batch
#'   processing* - blank filtering and imputation are each applied within a
#'   partition, so blanks and per-feature minima are taken from the same
#'   acquisition as the samples they apply to.
#'   `NULL` (default) reuses `bc_header`, which is what earlier versions did
#'   implicitly; `FALSE` processes every sample together. This is a separate
#'   question from which column defines the batch-correction reference, even
#'   when the same column answers both.
#' @param blank_head,blank_name `sample_meta` column and value identifying
#'   blank injections. The defaults are the canonical names this package
#'   documents (`Sample_type` / `Blank`), not a particular laboratory's sheet -
#'   set them to whatever your workbook uses.
#' @param filter_features `FALSE` (default) to skip, or a named list of
#'   arguments for [filterFeatures()] - `min_frac`, `method`, `group_col`,
#'   `sample_col`, `sample_labels`, `qc_label`. Applied once to the whole
#'   dataset after blank masking and before imputation, never per batch.
#' @param filter_samples `FALSE` (default) to skip, or a named list of
#'   arguments for [filterSamples()] - `max_na`, `sample_col`,
#'   `sample_labels`. Runs after `filter_features`, on the reduced feature
#'   set, so the missing-value fraction means something.
#' @param normalize Metadata column used as divisor, or `FALSE`/`NULL` for
#'   none.
#' @param adjust_conc Logical master toggle for concentration adjustment via
#'   [adjustConcentration()], using `sample_vol_col`, `cal_vol_col`,
#'   `sample_IS_col`, `cal_IS_col` and `starting_vol_col`.
#' @param starting_vol_col `FALSE` to skip the starting-volume correction, or
#'   the `metadata` column holding each sample's pre-dry-down starting volume.
#' @param sample_vol_col,cal_vol_col,sample_IS_col,cal_IS_col `metadata`
#'   columns holding, respectively: the vial reconstitution volume, the
#'   calibration-standard vial volume, the IS volume added to the sample,
#'   and the IS volume added to the calibration standard. Only consulted
#'   when `adjust_conc = TRUE`; different studies name these columns
#'   differently, hence configurable rather than hardcoded.
#' @param name_col `sample_metadata` column whose values match the data-matrix
#'   row identifiers (LC-MS injection / sample names). Default `"Name"`.
#' @param include_col,include_value `sample_metadata` inclusion filter: rows
#'   where `metadata[[include_col]] == include_value` are processed. Defaults
#'   `"Include"` / `"YES"`.
#' @param compound_col,processing_name_col `feature_metadata` columns holding
#'   the canonical feature name (becomes the data-matrix column names) and the
#'   raw name that joins to the source data. Defaults `"Compound"` /
#'   `"Processing_name"`.
#' @param report_col,report_value `feature_metadata` inclusion filter: features
#'   where `fdata[[report_col]] == report_value` are reported. Defaults
#'   `"Report"` / `"YES"`. Optional - when `report_col` is not a column of
#'   `fdata` every feature is reported, which is the normal case for a feature
#'   table curated outside a TargetLynx workflow.
#' @param comment_col `feature_metadata` column used to record why a feature
#'   was excluded. Default `"Comment"`.
#' @return A named list with `dataset` (the `DatasetExperiment`, imputed when
#'   `replace_MVs` is set), `excluded_features` (a data frame of features
#'   dropped along the way, with the reason recorded) and `measured` (the same
#'   dataset before imputation, gaps included; identical to `dataset` when
#'   nothing was imputed). The first two stay in that order, so existing code
#'   indexing `[[1]]` and `[[2]]` keeps working.
#' @examples
#' xlsx  <- system.file("extdata", "example_data.xlsx", package = "MRManalyzeR")
#' fdata <- openxlsx::read.xlsx(xlsx, sheet = "feature_metadata")
#' meta  <- openxlsx::read.xlsx(xlsx, sheet = "sample_metadata")
#' out <- processDataset(fdata, meta, xlsx_path = xlsx,
#'                        data_source = "targetlynx", datatype = "Area",
#'                        data_tab_names = "lcms_data", snr = 3,
#'                        blank_filter = FALSE, bc_header = "Chrom_Batch",
#'                        blank_head = "Sample_type")
#' out[[1]]
#' @family peak-matrix processing
#' @seealso [runMRManalyzeR()] to drive the whole workflow from a YAML config.
#' @export
processDataset = function(fdata,
                              metadata,
                              xlsx_path = NULL,
                              data_source = "targetlynx",
                              data_tab_names = NULL,
                              datatype = "Area",
                              tl_headers = c("ID", "Name", "Area", "ng/mL", "Response", "S/N"),
                              signal_filter = NULL,
                              snr = 5,
                              blank_filter = 5,
                              replace_MVs = FALSE,
                              impute_method = "gaussian",
                              impute_seed = 42,
                              batch_correction = FALSE,
                              bc_qc_label = "QC",
                              bc_factor_name = "Sample_type",
                              bc_header = "Chrom_Batch",
                              bc_check_factor = NULL,
                              processing_batch = NULL,
                              matrix_id_col = NULL,
                              matrix_orientation = "samples_rows",
                              blank_head = "Sample_type",
                              filter_features = FALSE,
                              filter_samples = FALSE,
                              blank_name = "Blank",
                              normalize = FALSE,
                              adjust_conc = FALSE,
                              starting_vol_col = FALSE,
                              sample_vol_col = "sample_volume_uL",
                              cal_vol_col = "cal_vol_uL",
                              sample_IS_col = "sample_IS_vol_uL",
                              cal_IS_col = "cal_IS_vol_uL",
                              name_col = "Name",
                              include_col = "Include",
                              include_value = "YES",
                              compound_col = "Compound",
                              processing_name_col = "Processing_name",
                              report_col = "Report",
                              report_value = "YES",
                              comment_col = "Comment"){

  for(.c in c(name_col, include_col)){
    if(!.c %in% colnames(metadata))
      stop(sprintf("[processDataset] sample_metadata has no '%s' column.", .c))
  }
  for(.c in c(compound_col, processing_name_col)){
    if(!.c %in% colnames(fdata))
      stop(sprintf("[processDataset] feature_metadata has no '%s' column.", .c))
  }

  # Canonicalise feature_metadata to the internal standard columns/values
  # (Compound / Processing_name / Report=YES|NO / Comment) so the rest of the
  # package - engine, reports, subsetDataset, combine - keeps using fixed names
  # regardless of what the input sheet calls them.
  fdata$Processing_name = fdata[[processing_name_col]]
  fdata$Compound        = fdata[[compound_col]]

  # Report is an inclusion filter, not measured data. A feature table curated
  # outside a TargetLynx workflow - the usual case for data_source = "matrix" -
  # carries no such column, and the only reading of "no exclusion information"
  # is that every feature is reported. Said out loud rather than assumed.
  fdata$Report = if(report_col %in% colnames(fdata)){
    ifelse(as.character(fdata[[report_col]]) == as.character(report_value),
           "YES", "NO")
  } else {
    message(sprintf(
      "[processDataset] feature_metadata has no '%s' column; reporting all %d features.",
      report_col, nrow(fdata)))
    rep("YES", nrow(fdata))
  }

  if(comment_col %in% colnames(fdata)) fdata$Comment = fdata[[comment_col]]

  if(!data_source %in% c("targetlynx", "skyline", "matrix"))
    stop("[processDataset] data_source must be 'targetlynx', 'skyline' or 'matrix'.")
  # Checked now rather than when imputation is reached, at the end of the run.
  impute_method = match.arg(impute_method, c("gaussian", "constant"))

  # Backwards compatibility: YAMLs that only set `snr:` (no `signal_filter:`)
  # keep working unchanged - infer the filter mode from `snr`.
  if(is.null(signal_filter))
    signal_filter = if(isFALSE(snr)) FALSE else "SNR"

  if(!identical(signal_filter, FALSE) && !signal_filter %in% c("SNR", "LOD", "LOQ"))
    stop("[processDataset] signal_filter must be one of: 'SNR', 'LOD', 'LOQ', FALSE.")

  if(identical(signal_filter, "SNR") && identical(data_source, "skyline"))
    stop("[processDataset] signal_filter='SNR' requires per-injection S/N values, which are not available for data_source='skyline'. Use signal_filter='LOD' or 'LOQ' instead.")

  if(signal_filter %in% c("LOD", "LOQ") &&
     !data_source %in% c("skyline", "matrix"))
    stop(sprintf(
      "[processDataset] signal_filter='%s' needs a per-compound threshold and is supported for data_source='skyline' or 'matrix'. For 'targetlynx', use signal_filter='SNR'.",
      signal_filter))

  if(signal_filter %in% c("LOD", "LOQ") && !signal_filter %in% colnames(fdata))
    stop(sprintf(
      "[processDataset] signal_filter='%s' but feature_metadata has no '%s' column.",
      signal_filter, signal_filter))

  # --- Batch-correction selection -----------------------------------------
  # TRUE corrects every feature; a string names an fdata column of YES / NO.
  # Checked before the data sheets are read, so a typo costs seconds rather
  # than a whole processing run. The header is matched as typed or in the
  # make.names() form assembleDataset() stores it under - correctBatch()
  # accepts both, so this check has to as well.
  bc_on = isTRUE(batch_correction) || is.character(batch_correction)
  if(is.character(batch_correction)){
    bc_col = NA_character_
    if(length(batch_correction) == 1L)
      bc_col = if(batch_correction %in% colnames(fdata)) batch_correction else
        colnames(fdata)[match(make.names(batch_correction),
                              make.names(colnames(fdata)))]
    if(is.na(bc_col))
      stop(sprintf(
        "[processDataset] batch_correction = '%s' is not True, False or a feature_metadata column (got: %s).",
        paste(batch_correction, collapse = ", "),
        paste(colnames(fdata), collapse = ", ")), call. = FALSE)
    invisible(.yes_no(fdata[[bc_col]][fdata$Report == "YES"],
                      batch_correction))
  }

  # --- Validation ---------------------------------------------------------
  # Only complain about a column the run will actually consult: warning that
  # blank_head is missing when blank_filter is off is noise, and noise is what
  # makes the real warnings easy to skip past.
  needed_heads = c(if(bc_on) c(bc_factor_name, bc_header),
                   if(!isFALSE(blank_filter) &&
                      !is.null(blank_filter)) blank_head)
  for(h in unique(needed_heads)){
    if(is.null(.resolve_meta_col(h, metadata)))
      warning(sprintf("'%s' not a metadata column.", h))
  }

  # --- Whitelists ----------------------------------------------------------
  fnames = fdata$Processing_name[fdata$Report == "YES"]
  snames = metadata[[name_col]][which(metadata[[include_col]] == include_value)]

  # --- Sheets -------------------------------------------------------------
  # NULL means each reader's usual sheets. TargetLynx keeps NULL, which
  # extractTable() reads as every sheet whose name contains "lcms_data".
  if(is.null(data_tab_names) && !identical(data_source, "targetlynx"))
    data_tab_names = if(identical(data_source, "skyline"))
      "skyline_data" else "matrix_data"

  # --- Read raw matrix (format-specific) + apply signal filter ------------
  if(identical(data_source, "targetlynx")){

    # readTargetLynx() reads + pivots + (optionally) SNR-masks the full sheet;
    # restrict to reportable compounds (fnames) and included samples (snames)
    # here, exactly as the old pre-pivot filter did.
    snr_use = if(identical(signal_filter, "SNR")) snr else FALSE
    tl = .read_targetlynx(xlsx_path, datatype = datatype,
                          tl_headers = tl_headers, snr = snr_use,
                          data_tab_names = data_tab_names)
    keep_r    = rownames(tl$masked) %in% snames
    keep_c    = colnames(tl$masked) %in% fnames
    out_table = tl$masked[keep_r, keep_c, drop = FALSE]
    raw_table = tl$raw[keep_r, keep_c, drop = FALSE]

  } else if(identical(data_source, "skyline")){

    if(is.null(xlsx_path))
      stop("[processDataset] xlsx_path is required when data_source='skyline'.")

    # readSkyline()'s two steps, split so the unmasked matrix is kept too.
    raw_table = .read_skyline_matrix(xlsx_path, data_tab_names, fdata)
    out_table = raw_table
    if(signal_filter %in% c("LOD", "LOQ"))
      out_table = .apply_lod_loq_mask(raw_table, fdata,
                                      floor_col = signal_filter)
    keep_r    = rownames(out_table) %in% snames
    out_table = out_table[keep_r, , drop = FALSE]
    raw_table = raw_table[keep_r, , drop = FALSE]

  } else {   # data_source == "matrix"

    if(is.null(xlsx_path))
      stop("[processDataset] xlsx_path is required when data_source='matrix'.")

    raw_table = readPeakMatrix(xlsx_path, data_tab_names,
                                 id_col      = matrix_id_col,
                                 orientation = matrix_orientation)
    # A generic matrix carries no S/N, but it may still have a per-compound
    # LOD/LOQ in feature_metadata, so the same mask applies.
    out_table = raw_table
    if(signal_filter %in% c("LOD", "LOQ"))
      out_table = .apply_lod_loq_mask(out_table, fdata,
                                      floor_col = signal_filter)
    keep_r    = rownames(out_table) %in% snames
    out_table = out_table[keep_r, , drop = FALSE]
    raw_table = raw_table[keep_r, , drop = FALSE]
  }

  # --- Zeros are non-detects ---------------------------------------------
  # TargetLynx reports an area of 0 when it finds no peak: nothing was
  # measured, rather than zero. Set to NA before any step counts it, so the
  # detection filters see it as missing and an unimputed matrix shows a gap.
  # raw_table too, so an all-zero compound is recorded as never measured.
  out_table[!is.na(out_table) & out_table == 0] = NA
  raw_table[!is.na(raw_table) & raw_table == 0] = NA

  # --- Mismatch diagnostics ----------------------------------------------
  # Track features that exist in fdata$Processing_name but have no row in
  # the source data (compound never reported for any sample, regardless of
  # data_source). Tracked here so they can be surfaced in the QC report and
  # the console.
  missing_in_data = setdiff(fnames, colnames(out_table))
  if(length(missing_in_data)){
    message(sprintf(
      "[processDataset] %d feature(s) in fdata$Processing_name have no rows in the source data:\n  %s",
      length(missing_in_data),
      paste(missing_in_data, collapse = ", ")))
  }

  # Rename Processing_name -> Compound via O(1) lookup. Any source column
  # not present in fdata$Processing_name maps to NA - those columns carry
  # data but have no feature_meta row, so we must drop them before
  # assembling the DatasetExperiment (otherwise SE complains that assay
  # colnames disagree with rowData rownames). For data_source="skyline",
  # this is also where unmatched Molecule names (e.g. internal standards
  # with no feature_metadata row) get caught.
  orig_names  = colnames(out_table)
  name_lookup = stats::setNames(fdata$Compound, fdata$Processing_name)
  new_names   = name_lookup[orig_names]
  unmapped    = is.na(new_names)
  unmapped_ids = orig_names[unmapped]   # surfaced in QC report
  if(length(unmapped_ids)){
    message(sprintf(
      "[processDataset] %d source column(s) have no matching fdata$Processing_name and were dropped:\n  %s",
      length(unmapped_ids),
      paste(unmapped_ids, collapse = ", ")))
    out_table = out_table[, !unmapped, drop = FALSE]
    new_names = new_names[!unmapped]
  }
  colnames(out_table) = new_names

  # Drop all-NA features (vectorized via colSums)
  keep_col = colSums(!is.na(out_table)) > 0
  remove_feats = names(out_table)[!keep_col]
  # Was each feature measured at all in the included samples, before any
  # signal filter masked it? Separates "never measured" from "filtered out"
  # in the reason recorded below. raw_table has out_table's rows and columns.
  had_value = stats::setNames(
    colSums(!is.na(raw_table[, orig_names[!unmapped], drop = FALSE])) > 0,
    colnames(out_table))
  out_table    = out_table[, keep_col, drop = FALSE]

  # --- Per-batch processing ----------------------------------------------
  meta_yes = metadata[which(metadata[[include_col]] == include_value), , drop = FALSE]

  # `bc_header` names the batch-correction column, but historically it also
  # silently partitioned every processing step. Those are different questions,
  # so `processing_batch` now asks the second one explicitly - defaulting to
  # bc_header, which preserves the previous behaviour.
  asked_for  = !is.null(processing_batch)
  proc_batch = if(asked_for) processing_batch else bc_header

  if(!isFALSE(proc_batch) && !is.null(proc_batch) &&
     is.null(.resolve_meta_col(proc_batch, meta_yes))){
    if(asked_for){
      # Named outright, so a missing column is a mistake worth stopping for.
      stop(sprintf(
        "processing_batch='%s' is not a column in sample_metadata. Set it to an existing column, or to False to process every sample together (got: %s).",
        proc_batch, paste(colnames(meta_yes), collapse = ", ")))
    }
    # Inherited from bc_header rather than requested. A study with no batch
    # column at all is an ordinary single-batch run, and failing it over a
    # value nobody typed would be wrong - so say what is happening and carry
    # on with every sample in one batch.
    warning(sprintf(
      "processing_batch defaulted to bc_header='%s', which is not a sample_metadata column; processing every sample together. Set processing_batch explicitly to silence this.",
      proc_batch), call. = FALSE)
    proc_batch = FALSE
  }
  if(is.null(proc_batch)) proc_batch = FALSE

  # The correction column is only needed if correction is actually requested.
  if(bc_on &&
     (is.null(bc_header) || isFALSE(bc_header) ||
      is.null(.resolve_meta_col(bc_header, meta_yes)))){
    stop(sprintf(
      "batch_correction is on but bc_header='%s' is not a column in sample_metadata (got: %s).",
      bc_header %||% "<NULL>",
      paste(colnames(meta_yes), collapse = ", ")))
  }

  # --- Assemble the DatasetExperiment, then process it --------------------
  # The DE is built first so every processing step downstream operates on the
  # dataset object (data + sample_meta together) rather than a bare matrix.
  if(nrow(out_table) == 0 || ncol(out_table) == 0){
    stop(sprintf(
      "[processDataset] Empty matrix after reading (%d samples x %d features). Check: processing_batch column, blank/signal filters, and that the source data actually contains data for the included samples.",
      nrow(out_table), ncol(out_table)))
  }
  lcms_experiment = assembleDataset(out_table, fdata, metadata, name_col = name_col)

  # assembleDataset() stores headings in make.names() form (Chrom-Batch ->
  # Chrom.Batch), so look the processing-batch column up under that name.
  if(!isFALSE(proc_batch)){
    smeta_st   = as.data.frame(lcms_experiment$sample_meta)
    proc_batch = .resolve_meta_col(proc_batch, smeta_st) %||% proc_batch
  }

  # Treat NA / "" bc_header values (e.g. blanks, QCs that weren't
  # cryosectioned) as a single "_unbatched_" group so they survive the
  # per-batch split instead of being silently dropped via `NA == "x"`.
  # Computed locally: the dataset's own sample_meta keeps the original values.
  bc_vec = if(isFALSE(proc_batch))
    rep("_all_", nrow(as.data.frame(lcms_experiment$data))) else
    as.character(as.data.frame(lcms_experiment$sample_meta)[[proc_batch]])
  bc_vec[is.na(bc_vec) | !nzchar(trimws(bc_vec))] = "_unbatched_"
  batches = unique(bc_vec)

  # --- pass 1: blank filter, per batch ----------------------------------
  # Blank thresholds are a property of a batch, so this stays per-batch. It
  # only masks values, never removes a row or column, so the rebind below is
  # safe.
  lcms_experiment = .run_per_batch(
    lcms_experiment, bc_vec, batches,
    function(d) if(isFALSE(blank_filter) || is.null(blank_filter)) d else
      filterBlanks(d, blank_filter = blank_filter,
                   blank_head = blank_head, blank_name = blank_name))

  # --- detection filters, GLOBAL ----------------------------------------
  # Deliberately not per batch: dropping a feature in one batch and keeping it
  # in another gives ragged column sets that cannot be rebound, and a feature's
  # detection is a property of the study rather than of a batch. Runs here
  # because it needs the NA pattern intact - imputation below destroys it.
  # Features first, then samples: judging an injection against a feature list
  # still padded with undetected features gives a meaningless denominator.
  if(!isFALSE(filter_features)){
    ff = as.list(filter_features)
    lcms_experiment = filterFeatures(
      lcms_experiment,
      min_frac      = ff$min_frac      %||% 0.5,
      method        = ff$method        %||% "within",
      group_col     = ff$group_col,
      sample_col    = ff$sample_col    %||% "Sample_type",
      sample_labels = ff$sample_labels %||% "Sample",
      qc_label      = ff$qc_label      %||% "QC")
  }
  if(!isFALSE(filter_samples)){
    fs = as.list(filter_samples)
    lcms_experiment = filterSamples(
      lcms_experiment,
      max_na        = fs$max_na        %||% 0.8,
      sample_col    = fs$sample_col    %||% "Sample_type",
      sample_labels = fs$sample_labels %||% "Sample")
  }

  # --- batch correction, GLOBAL ------------------------------------------
  # Straight after the filters, on the instrument signal: a batch effect is a
  # property of the measurement, and normalising first would rescale the QC
  # reference by per-sample amounts pooled QCs usually do not have. Before
  # imputation too, so a batch median describes what was measured rather than
  # values imputation filled in.
  if(bc_on)
    lcms_experiment = correctBatch(
      lcms_experiment,
      qc_label     = bc_qc_label,
      check_factor = bc_check_factor,
      factor_name  = bc_factor_name,
      batch_head   = bc_header,
      feature_col  = if(is.character(batch_correction)) batch_correction)

  # --- normalise / concentration adjustment --------------------------------
  # Per-sample divisors and factors, so no batch split is needed.
  lcms_experiment = .scale_samples(
    lcms_experiment,
    normalize        = normalize,
    adjust_conc      = adjust_conc,
    starting_vol_col = starting_vol_col,
    sample_vol_col   = sample_vol_col,
    cal_vol_col      = cal_vol_col,
    sample_IS_col    = sample_IS_col,
    cal_IS_col       = cal_IS_col)

  # --- imputation, per processing batch, LAST ------------------------------
  # Last, because it destroys the missingness every step above relies on. The
  # dataset as measured is kept as well: the reports use it by default.
  measured = lcms_experiment
  if(!isFALSE(replace_MVs) && !is.null(replace_MVs))
    lcms_experiment = imputeMissing(
      lcms_experiment, scalar = replace_MVs,
      blank_head = blank_head, blank_name = blank_name,
      method = impute_method, seed = impute_seed,
      by = if(isFALSE(proc_batch)) NULL else proc_batch)

  # --- Record removed features -------------------------------------------
  # Two different reasons: nothing was measured in the included samples, or
  # something was and the signal filter removed all of it.
  never_seen = remove_feats[!had_value[remove_feats] %in% TRUE]
  filtered   = setdiff(remove_feats, never_seen)
  fdata$Report [fdata$Compound %in% remove_feats] = "NO"
  fdata$Comment[fdata$Compound %in% never_seen] =
    "No values in the included samples"
  if(length(filtered))
    fdata$Comment[fdata$Compound %in% filtered] =
      sprintf("All values below %s threshold", signal_filter)

  # Flag features whose Processing_name appears in fdata but has no rows
  # in the source data (compound never reported for any sample).
  if(length(missing_in_data)){
    miss_row = fdata$Processing_name %in% missing_in_data
    fdata$Report [miss_row] = "NO"
    fdata$Comment[miss_row] = "No matching data in source file"
  }
  removed_features = dplyr::filter(fdata, Report == "NO")

  # Append unmapped source columns (no fdata entry at all) so the QC report
  # also surfaces those - fdata has no row for them, so we add stub rows.
  if(length(unmapped_ids)){
    extra = data.frame(
      Compound = unmapped_ids,
      Report   = "NO",
      Comment  = "No matching Processing_name in feature_metadata",
      stringsAsFactors = FALSE
    )
    # pad to the same columns as removed_features
    miss_cols = setdiff(colnames(removed_features), colnames(extra))
    for(m in miss_cols) extra[[m]] = NA
    extra = extra[, colnames(removed_features), drop = FALSE]
    removed_features = rbind(removed_features, extra)
  }

  # Console summary banner so mismatches are impossible to miss
  if(length(unmapped_ids) || length(missing_in_data)){
    message("\n--- Feature-mismatch summary ---")
    message(sprintf("  fdata Processing_name with no source data : %d", length(missing_in_data)))
    message(sprintf("  source columns with no fdata entry        : %d", length(unmapped_ids)))
    message("  (both groups listed in 'Compounds excluded' in the QC report)\n")
  }

  # Named, but the order is unchanged so out[[1]] / out[[2]] still work.
  list(dataset = lcms_experiment, excluded_features = removed_features,
       measured = measured)
}


# ---- Internal helpers ------------------------------------------------------

#' Pivot the long TargetLynx table to a wide sample x compound matrix
#' @keywords internal
#' @noRd
.build_wide_matrix = function(base, datatype){
  base %>%
    dplyr::select(ID, Name, dplyr::all_of(datatype)) %>%
    dplyr::distinct() %>%
    dplyr::mutate(dplyr::across(dplyr::all_of(datatype), as.numeric)) %>%
    # values_fill = NA, not 0. A sample x compound pair with no row in the
    # export was not measured as zero - it was not reported at all, whether
    # because the transition was not acquired, the peak was not integrated, or
    # the row was dropped from the export. Filling with 0 makes those absences
    # indistinguishable from a genuine zero and lets them behave as real
    # measurements in blank filtering, normalisation and concentration
    # adjustment, which all run before any imputation step.
    tidyr::pivot_wider(names_from = "ID",
                       values_from = dplyr::all_of(datatype),
                       id_cols = "Name", values_fill = NA_real_) %>%
    tibble::column_to_rownames("Name")
}

#' Apply an S/N threshold mask to a wide matrix
#' @keywords internal
#' @noRd
.apply_snr_mask = function(out_table, base, snr){
  snr_table = base %>%
    dplyr::select(ID, Name, `S/N`) %>%
    dplyr::distinct() %>%
    dplyr::mutate(`S/N` = suppressWarnings(as.numeric(`S/N`))) %>%
    tidyr::pivot_wider(names_from = "ID", values_from = "S/N",
                       id_cols = "Name") %>%
    tibble::column_to_rownames("Name")

  snr_table = snr_table[rownames(out_table), colnames(out_table), drop = FALSE]
  out_table[snr_table < snr] = NA
  out_table
}

#' Normalise and concentration-adjust every sample
#'
#' Both steps are per-sample divisors or factors, so they need no batch
#' split. They run after batch correction and before imputation - see
#' processDataset() for the order and why.
#'
#' Operates on a `DatasetExperiment` so each step reads what it needs
#' (divisors, volumes) straight from `sample_meta`.
#' @keywords internal
#' @noRd
.scale_samples = function(de, normalize, adjust_conc, starting_vol_col,
                          sample_vol_col, cal_vol_col,
                          sample_IS_col, cal_IS_col){

  if(!isFALSE(normalize) && !is.null(normalize))
    de = normaliseMatrix(de, column = normalize)

  # Vial-level concentration correction (+ optional starting-volume
  # correction) -- see adjustConcentration.R.
  if(isTRUE(adjust_conc))
    de = adjustConcentration(de,
                              sample_vol_col   = sample_vol_col,
                              cal_vol_col      = cal_vol_col,
                              sample_IS_col    = sample_IS_col,
                              cal_IS_col       = cal_IS_col,
                              starting_vol_col = starting_vol_col)

  de
}

#' Apply a per-batch transform and rebind, preserving row order
#'
#' `fn` must return a DatasetExperiment with the same rows and columns it was
#' given - it may change values but not the shape, because the batches are
#' rebound on a common column set. Single-batch datasets skip the split.
#' @keywords internal
#' @noRd
.run_per_batch = function(de, bc_vec, batches, fn){
  if(length(batches) <= 1L) return(fn(de))
  row_order = rownames(as.data.frame(de$data))
  frames = lapply(batches, function(b)
    as.data.frame(fn(.de_rows(de, which(bc_vec == b)))$data))
  merged = do.call(rbind, frames)
  de$data = merged[row_order, , drop = FALSE]
  de
}

#' Subset a DatasetExperiment to the given sample rows, keeping all features
#' @keywords internal
#' @noRd
.de_rows = function(de, idx){
  struct::DatasetExperiment(
    data          = as.data.frame(de$data)[idx, , drop = FALSE],
    sample_meta   = as.data.frame(de$sample_meta)[idx, , drop = FALSE],
    variable_meta = as.data.frame(de$variable_meta)
  )
}
