#' Run the post-acquisition MRManalyzeR workflow from a YAML config
#'
#' Single entry point that:
#' \itemize{
#'   \item loads a project YAML,
#'   \item (optionally) reads the TargetLynx xlsx workbook and builds the
#'     processed peak matrix - gated by `PeakMatrixProcessing.execute`,
#'   \item writes the processed matrix to xlsx + RDS and persists the YAML
#'     parameters alongside - imputed when `replace_MVs` is set, with the
#'     measured values added as a `matrix_measured` sheet and
#'     `<name>_measured.RDS`,
#'   \item runs the statistical analyses defined under `stats_report:`
#'     (comparisons / correlations / linear_models) and writes them as
#'     extra tabs in the results xlsx,
#'   \item (optionally) renders the bundled `data_quality_report` and
#'     `stats_report` HTML vignettes.
#' }
#'
#' Per-compound quality metrics are written into the returned dataset's
#' `variable_meta`, matched to the feature they describe: `CV_QC` (coefficient
#' of variation across the pooled QC injections - technical variability, and
#' the usual quantification filter), `CV_sample` (the same across the
#' biological samples, so biological spread plus that technical noise), and
#' `CV_sample_vs_QC`. A ratio near 1 flags a compound varying no more between
#' samples than between replicate injections of identical material. `CV_QC` is
#' `NA` when a run contains no QC injections. They are computed on the measured
#' values, never on imputed ones.
#'
#' Both reports and the statistics use the measured values unless their block
#' sets `use_imputed: True`; the PCAs impute for themselves.
#'
#' Output filenames are derived from `paths.fn` + `datatype` +
#' `paths.suffix`, so when `PeakMatrixProcessing.execute: False` the code can
#' still locate the existing RDS without re-reading the source xlsx. Those
#' three must therefore match the stored file; if they do not, the error
#' lists the `.RDS` files that are in `paths.result_dir`.
#'
#' `datatype` is read from the enabled data-source block: `datatype:` under
#' `tl_data:` / `matrix_data:`, and `signal_filter:` under `skyline_data:`,
#' where the LOD / LOQ floor is what distinguishes one run from another. Any
#' of them may be a list, in which case the whole workflow runs once per entry.
#'
#' Backwards compatibility:
#' \itemize{
#'   \item `datatype` falls back to `PeakMatrixProcessing:` and then `paths:`,
#'     where earlier configs put it. Setting it in two places warns.
#'   \item Legacy `UVA_report:` / `MVA_report:` blocks are recognised when
#'     `data_quality_report:` / `stats_report:` are missing.
#' }
#'
#' @param path_yaml Full path to a project YAML.
#' @return Invisibly, a list with the `DatasetExperiment` (imputed when
#'   `replace_MVs` is set), the `measured` one, `removed_features`,
#'   `stats_tables`, and the resolved output paths. When `datatype` is a
#'   vector (e.g. `["Area", "Response", "Conc"]`), the function loops over each
#'   datatype and returns a list of per-datatype results.
#' @examples
#' # runExample() writes a config for the bundled dataset and drives
#' # runMRManalyzeR(); render = FALSE keeps it to the peak-matrix build.
#' runExample(render = FALSE, open = FALSE)
#' @family entry points
#' @export
runMRManalyzeR = function(path_yaml){

  stopifnot(file.exists(path_yaml))

  # --- Load & unpack YAML -------------------------------------------------
  project_params = loadConfig(path_yaml)

  # Fail early, and once. A structural mistake in the config would otherwise
  # surface as an obscure error somewhere in processing or report rendering,
  # after minutes of work, with only the first problem visible.
  cfg_check = validateConfig(project_params)
  if(length(cfg_check$errors))
    stop("[runMRManalyzeR] configuration is not valid:
",
         paste0("  - ", cfg_check$errors, collapse = "
"),
         call. = FALSE)
  for(w in cfg_check$warnings) message("[validateConfig] ", w)

  project_paths = project_params$project$paths
  pmp_params    = project_params$project$PeakMatrixProcessing

  # Prefer new combined blocks; fall back to legacy UVA/MVA blocks
  dq_params = project_params$project$data_quality_report %||%
                project_params$project$UVA_report
  st_params = project_params$project$stats_report %||%
                project_params$project$MVA_report


  # --- Datatype loop ------------------------------------------------------
  # Accept either a scalar ("Area") or a vector (["Area", "Response", "Conc"]).
  datatypes = .resolve_datatype(pmp_params, project_paths)
  if(length(datatypes) > 1){
    results = lapply(datatypes, function(dt){
      message(sprintf("\n========== Datatype: %s ==========", dt))
      pp_one  = project_paths;            pp_one$datatype = dt
      pmp_one = .narrow_datatype(pmp_params, dt)
      pp_proj_one = project_params
      pp_proj_one$project$paths = pp_one
      pp_proj_one$project$PeakMatrixProcessing = pmp_one
      .run_MRManalyzeR_one(project_params = pp_proj_one,
                           project_paths  = pp_one,
                           pmp_params     = pmp_one,
                           dq_params      = dq_params,
                           st_params      = st_params,
                           datatype       = dt)
    })
    names(results) = datatypes
    return(invisible(results))
  }
  datatype = datatypes[[1]]

  invisible(.run_MRManalyzeR_one(
    project_params = project_params,
    project_paths  = project_paths,
    pmp_params     = pmp_params,
    dq_params      = dq_params,
    st_params      = st_params,
    datatype       = datatype
  ))
}


#' Resolve `datatype` - the reported value - from the config
#'
#' Read from the enabled data-source block first, because the thing that
#' distinguishes one run from another differs by vendor:
#' \itemize{
#'   \item `tl_data:` and `matrix_data:` use `datatype:` - a TargetLynx units
#'     tag ("Area", "Response", "ng/mL"), or a label for a generic matrix.
#'   \item `skyline_data:` uses `signal_filter:`. A Skyline export carries no
#'     units tag, so the LOD / LOQ floor being applied is the only thing that
#'     separates one run from another, and a second key naming the same
#'     distinction would only be a second place to get it wrong.
#' }
#' `PeakMatrixProcessing:` and then `paths:` are still honoured, because that
#' is where earlier configs put it.
#'
#' @param pmp_params The `PeakMatrixProcessing:` block.
#' @param project_paths The `paths:` block.
#' @return Character vector of one or more datatypes.
#' @keywords internal
#' @noRd
.resolve_datatype = function(pmp_params, project_paths){

  sources = c("tl_data", "skyline_data", "matrix_data")
  enabled = sources[vapply(sources, function(s)
    isTRUE(pmp_params[[s]]$enabled), logical(1))]

  # Exactly one source is expected to be enabled; processDataset() is what
  # enforces that, so here an ambiguous config simply falls through to the
  # block-level setting rather than guessing which reader was meant.
  from_source = if(length(enabled) != 1L) NULL
                else .source_datatype(enabled, pmp_params[[enabled]])
  from_block  = pmp_params$datatype %||% project_paths$datatype

  if(!is.null(from_source) && !is.null(from_block) &&
     !identical(from_source, from_block))
    warning(sprintf(
      "datatype is set in both %s (%s) and above it (%s). Using the %s value.",
      enabled, paste(from_source, collapse = ", "),
      paste(from_block, collapse = ", "), enabled), call. = FALSE)

  dt = from_source %||% from_block
  if(is.null(dt)){
    # Silently defaulting here renames every output file, which surfaces much
    # later as "expected RDS not found" rather than as the missing key it is.
    warning(sprintf(
      "No datatype found in the config; defaulting to \"Area\", which names the output files. %s",
      if(identical(enabled, "skyline_data"))
        "Set signal_filter: in the skyline_data block, or datatype: above it."
      else if(length(enabled) == 1L)
        sprintf("Add `datatype:` to the %s block.", enabled)
      else "Enable exactly one data source and give it a `datatype:`."),
      call. = FALSE)
    dt = "Area"
  }
  as.character(dt)
}

#' The field a given source block uses to name its reported value
#' @keywords internal
#' @noRd
.source_datatype = function(source_name, blk){
  if(identical(source_name, "skyline_data")){
    sf = blk$signal_filter
    # FALSE is "apply no floor", which is a valid choice but not a label -
    # fall through so the run is named by datatype: instead of "FALSE".
    if(is.null(sf) || isFALSE(sf)) NULL else sf
  } else {
    blk$datatype
  }
}

#' Pin a multi-valued config down to the single datatype being run
#'
#' The loop runs the whole workflow once per entry, so the source block itself
#' has to be narrowed and not just the top-level copy: [readSkyline()] indexes
#' `feature_metadata` by `signal_filter`, and a two-element list is not a
#' column name.
#'
#' @param pmp_params The `PeakMatrixProcessing:` block.
#' @param dt One datatype.
#' @return `pmp_params`, narrowed.
#' @keywords internal
#' @noRd
.narrow_datatype = function(pmp_params, dt){
  pmp_params$datatype = dt
  if(isTRUE(pmp_params$skyline_data$enabled)){
    # Only when signal_filter is what supplied the datatype: a config that
    # deliberately applies no floor must keep applying none.
    sf = pmp_params$skyline_data$signal_filter
    if(!is.null(sf) && !isFALSE(sf))
      pmp_params$skyline_data$signal_filter = dt
  } else if(isTRUE(pmp_params$tl_data$enabled)){
    pmp_params$tl_data$datatype = dt
  } else if(isTRUE(pmp_params$matrix_data$enabled)){
    pmp_params$matrix_data$datatype = dt
  }
  pmp_params
}

#' Read one named sheet, naming the config key when it is not there
#'
#' `openxlsx::read.xlsx()` fails with "Cannot find sheet named ..." and leaves
#' the reader to work out which setting produced the name. Since the sheet is
#' now configurable, the error has to say which key to change and what the
#' workbook actually holds.
#'
#' @param xlsx_path Workbook path.
#' @param sheet Sheet name.
#' @param key The YAML key that supplied `sheet`, for the error message.
#' @return A data frame.
#' @keywords internal
#' @noRd
.read_sheet = function(xlsx_path, sheet, key){
  present = openxlsx::getSheetNames(xlsx_path)
  if(!sheet %in% present)
    stop(sprintf(
      "[runMRManalyzeR] no sheet '%s' in %s.\nSheets present: %s.\nSet PeakMatrixProcessing.%s to one of these.",
      sheet, basename(xlsx_path), paste(present, collapse = ", "), key),
      call. = FALSE)
  openxlsx::read.xlsx(xlsx_path, sheet = sheet)
}

#' Single-datatype runner extracted from [`runMRManalyzeR()`].
#' @keywords internal
#' @noRd
.run_MRManalyzeR_one = function(project_params, project_paths, pmp_params,
                                dq_params, st_params, datatype){

  # The resolved datatype is written back so the report templates see it
  # wherever the config happened to declare it - they fall back to it for the
  # axis / matrix value label when `units:` is unset - and so the source block
  # holds a single value rather than the list that was looped over.
  pmp_params = .narrow_datatype(pmp_params, datatype)

  units_tag = gsub("/", "_", datatype)
  out_stub  = sprintf("%s/%s_%s%s",
                      project_paths$result_dir,
                      project_paths$fn,
                      units_tag,
                      project_paths$suffix %||% "")

  out_xlsx       = paste0(out_stub, ".xlsx")
  out_stats_xlsx = paste0(out_stub, "_stats.xlsx")
  out_RDS        = paste0(out_stub, ".RDS")
  out_pars       = paste0(out_stub, ".txt")

  dir.create(project_paths$result_dir, recursive = TRUE, showWarnings = FALSE)

  # --- Peak-matrix processing (toggle) -----------------------------------
  if(isTRUE(pmp_params$execute)){
    out = .run_pmp(project_paths, pmp_params, datatype,
                   out_xlsx = out_xlsx, out_RDS = out_RDS, out_pars = out_pars,
                   project_params = project_params)
    combined_datamatrices = out$combined_datamatrices
    measured_datamatrices = out$measured
    removed_features      = out$removed_features
    imputed_flags = c(imputed  = .n_imputed_total(combined_datamatrices) > 0,
                      measured = FALSE)
  } else {
    if(!file.exists(out_RDS)){
      # The name is built from fn / datatype / suffix, so a mismatch is a
      # config error rather than a missing file. Listing what is actually in
      # result_dir turns "work out the name" into "read it off this list".
      avail = grep("_measured\\.RDS$",
                   list.files(project_paths$result_dir, pattern = "\\.RDS$"),
                   value = TRUE, invert = TRUE)
      stop(sprintf(
        "PeakMatrixProcessing.execute is FALSE but no stored dataset was found:\n  %s\n\n%s\n\n%s",
        out_RDS,
        if(length(avail))
          paste0("Available in ", project_paths$result_dir, ":\n  ",
                 paste(avail, collapse = "\n  "))
        else
          paste0("No .RDS files in ", project_paths$result_dir, "."),
        "The name is <fn>_<datatype><suffix>.RDS - set paths.fn, datatype and paths.suffix to match one of these, or re-run with execute: True."
      ))
    }
    message("Skipping PeakMatrixProcessing (execute: False). Loading: ", out_RDS)
    combined_datamatrices = readRDS(out_RDS)
    # Written beside it when the run imputed, so the reports can still use
    # the measured values.
    out_meas = .measured_rds(out_RDS)
    has_meas = file.exists(out_meas)
    measured_datamatrices = if(has_meas) readRDS(out_meas) else
      combined_datamatrices
    stored_imp    = .values_imputed(combined_datamatrices,
                                    pmp_params$replace_MVs)
    imputed_flags = c(imputed = stored_imp, measured = stored_imp && !has_meas)
    if(stored_imp && !has_meas)
      warning(sprintf(
        "[runMRManalyzeR] %s holds imputed values and no %s was found (it is written from version 0.99.4 on), so the reports and statistics will use the imputed values. Re-run with execute: True to get the measured ones.",
        basename(out_RDS), basename(out_meas)), call. = FALSE)
    removed_features      = data.frame()    # not available without re-processing
  }

  # --- Statistics ---------------------------------------------------------
  # Each stats section is a block with `enabled:` + `entries:`;
  # .section_entries() honours the toggle and returns the runnable entries.
  # Tests and summaries use the measured values unless the stats report
  # asks for the imputed matrix.
  st_data = if(isTRUE(st_params$use_imputed)) combined_datamatrices else
    measured_datamatrices
  stats_tables  = NULL
  group_summary = NULL
  if(isTRUE(st_params$execute)){
    has_stats =
      length(.section_entries(st_params$comparisons,   "comparisons"))   > 0 ||
      length(.section_entries(st_params$correlations,  "correlations"))  > 0 ||
      length(.section_entries(st_params$linear_models, "linear_models")) > 0 ||
      length(.ion_entries(st_params$ion_ratios))                          > 0
    if(has_stats){
      message("Running statistics...")
      stats_tables = runStats(st_data, st_params)

      # Group summaries are what the boxplots and barplots are drawn from, so
      # they belong in the workbook next to the tests.
      gsf = st_params$global_summary_factor
      if(!is.null(gsf) &&
         gsf %in% colnames(as.data.frame(st_data$sample_meta)))
        group_summary = tryCatch(
          summariseGroups(st_data, gsf),
          error = function(e){
            warning("[runMRManalyzeR] group summary skipped: ",
                    conditionMessage(e)); NULL
          })

      append_stats_xlsx(out_stats_xlsx, stats_tables, group_summary)
      message("Wrote stats to: ", out_stats_xlsx)
    }
  }

  # --- Optional reports ---------------------------------------------------
  .render_reports(project_params   = project_params,
                  project_paths    = project_paths,
                  pmp_params       = pmp_params,
                  dq_params        = dq_params,
                  st_params        = st_params,
                  out_RDS          = out_RDS,
                  out_stub         = out_stub,
                  combined_datamatrices = combined_datamatrices,
                  measured_datamatrices = measured_datamatrices,
                  removed_features      = removed_features,
                  stats_tables          = stats_tables,
                  imputed               = imputed_flags)

  invisible(list(
    datasetExperiment = combined_datamatrices,
    measured          = measured_datamatrices,
    removed_features  = removed_features,
    stats_tables      = stats_tables,
    # A named element rather than an attribute: this is the first thing
    # anyone wants after a run, and `res$output_directory` is discoverable
    # from str() and tab-completion in a way attr() is not.
    output_directory  = project_paths$result_dir,
    out_xlsx       = out_xlsx,
    out_stats_xlsx = out_stats_xlsx,
    out_RDS        = out_RDS,
    out_pars       = out_pars,
    datatype       = datatype
  ))
}


# --- Internal helpers -------------------------------------------------------

`%||%` = function(a, b) if(is.null(a)) b else a

#' Read a detection-filter block from the YAML
#'
#' Returns FALSE (the off switch processDataset() understands) when the block
#' is absent, or present with enabled: False. Otherwise hands back the block
#' itself, minus the `enabled` flag, so its keys map onto the filter's
#' arguments by name.
#' @keywords internal
#' @noRd
.filter_block = function(x){
  if(is.null(x) || isFALSE(x)) return(FALSE)
  x = as.list(x)
  if(!isTRUE(x$enabled)) return(FALSE)
  x$enabled = NULL
  x
}

#' The injection-order column to validate, or NULL when nothing reads it
#'
#' Only the data-quality report uses injection order (drift plots, QC-PCA
#' colour scale), so a run that does not draw it is not stopped over the
#' column. Falls back to a legacy `UVA_report:` block, as the report does.
#' @keywords internal
#' @noRd
.order_head = function(project_params){
  dq = project_params$project$data_quality_report %||%
       project_params$project$UVA_report %||% list()
  if(!isTRUE(dq$execute)) return(NULL)
  dq$injection_order_head %||% "Injection_order"
}

#' @keywords internal
#' @noRd
.run_pmp = function(project_paths, pmp_params, datatype,
                    out_xlsx, out_RDS, out_pars, project_params){

  xlsx_path = sprintf("%s/%s.xlsx", project_paths$data_dir, project_paths$fn)
  message("Reading ", xlsx_path)

  # Sample-metadata contract columns (overridable for studies that name them
  # differently). name_col matches LC-MS injection / data-matrix row names.
  name_col       = pmp_params$name_col        %||% "Name"
  include_col    = pmp_params$include_col     %||% "Include"
  include_value  = pmp_params$include_value   %||% "YES"
  sample_meta_tab = pmp_params$sample_meta_tab %||% "sample_metadata"

  # feature_metadata contract columns (canonicalised internally to
  # Compound / Processing_name / Report / Comment).
  compound_col        = pmp_params$compound_col        %||% "Compound"
  processing_name_col = pmp_params$processing_name_col %||% "Processing_name"
  report_col          = pmp_params$report_col          %||% "Report"
  report_value        = pmp_params$report_value        %||% "YES"
  comment_col         = pmp_params$comment_col         %||% "Comment"
  feature_meta_tab    = pmp_params$feature_meta_tab    %||% "feature_metadata"

  # Path-based reads avoid an openxlsx::loadWorkbook() bug on workbooks
  # with certain styling/drawing XML, and are also faster.
  fdata    = .read_sheet(xlsx_path, feature_meta_tab, "feature_meta_tab")
  metadata = .read_sheet(xlsx_path, sample_meta_tab,  "sample_meta_tab")

  # Check the workbook before anything reads it. Every problem is reported at
  # once: a duplicated sample name and a text value in a measurement column
  # are both worth knowing about before a ten-minute run, not one after the
  # other across three attempts.
  in_check = validateInput(fdata, metadata,
                            name_col            = name_col,
                            compound_col        = compound_col,
                            processing_name_col = processing_name_col,
                            report_col          = report_col,
                            include_col         = include_col,
                            include_value       = include_value,
                            # Drift plots and the QC-PCA read it as a number,
                            # so a text value is stopped here, not drawn empty
                            # - checked only when that report will be drawn.
                            injection_order_head = .order_head(project_params))
  if(length(in_check$errors))
    stop("[runMRManalyzeR] input workbook is not valid:
",
         paste0("  - ", in_check$errors, collapse = "
"), call. = FALSE)
  for(w in in_check$warnings) message("[validateInput] ", w)

  metadata = metadata[which(metadata[[include_col]] == include_value), , drop = FALSE]

  # Data source is selected by two mutually-exclusive nested blocks, each
  # with its own `enabled:` toggle and source-specific params:
  #   skyline_data: { enabled, data_tab_names, signal_filter }   # LOD/LOQ,
  #                                                              # doubles as datatype
  #   tl_data:      { enabled, data_tab_names, tl_headers, snr }  # SNR
  #   matrix_data:  { enabled, data_tab_names, id_col, orientation,
  #                   signal_filter }                             # LOD/LOQ
  # Exactly one must be enabled.
  sky = pmp_params$skyline_data %||% list()
  tl  = pmp_params$tl_data      %||% list()
  mx  = pmp_params$matrix_data  %||% list()
  on  = c(skyline_data = isTRUE(sky$enabled),
          tl_data      = isTRUE(tl$enabled),
          matrix_data  = isTRUE(mx$enabled))
  if(sum(on) != 1)
    stop(sprintf(
      "[runMRManalyzeR] %d data sources enabled under PeakMatrixProcessing (%s) - enable exactly one.",
      sum(on), if(any(on)) paste(names(on)[on], collapse = ", ") else "none"))

  matrix_id_col      = NULL
  matrix_orientation = "samples_rows"

  if(on[["skyline_data"]]){
    data_source    = "skyline"
    data_tab_names = sky$data_tab_names %||% "skyline_data"
    signal_filter  = sky$signal_filter %||% FALSE   # LOD | LOQ | False
    snr            = FALSE
    tl_headers     = NULL
  } else if(on[["matrix_data"]]){
    # Vendor-neutral: any software that exports a rectangular
    # sample x analyte table, with no format-specific parsing.
    data_source        = "matrix"
    data_tab_names     = mx$data_tab_names %||% "matrix_data"
    signal_filter      = mx$signal_filter  %||% FALSE  # LOD | LOQ | False
    snr                = FALSE
    tl_headers         = NULL
    matrix_id_col      = mx$id_col
    matrix_orientation = mx$orientation %||% "samples_rows"
  } else {
    data_source    = "targetlynx"
    data_tab_names = tl$data_tab_names               # NULL => extractTable greps "lcms_data"
    tl_headers     = tl$tl_headers
    snr            = tl$snr %||% FALSE               # numeric threshold, or FALSE to skip
    signal_filter  = if(isFALSE(snr)) FALSE else "SNR"
  }

  # PeakMatrixProcessing.adjust_conc is a nested YAML block (enabled toggle
  # + 5 sample_metadata column names) rather than flat keys, so the study's
  # own column names never need to be hardcoded in R -- see adjustConcentration.R.
  ac_pars = pmp_params$adjust_conc %||% list()

  message("Generating data matrix...")
  combined_data = processDataset(
    fdata            = fdata,
    metadata         = metadata,
    xlsx_path        = xlsx_path,
    data_source      = data_source,
    # NULL is meaningful: for tl_data it means every 'lcms_data*' sheet. The
    # other sources set their own default above.
    data_tab_names   = data_tab_names,
    datatype         = datatype,
    tl_headers       = tl_headers %||% c("ID", "Name", "Area", "ng/mL", "Response", "S/N"),
    signal_filter    = signal_filter,
    snr              = snr,
    # An absent key means the step is off. Passing NULL through would reach
    # the step itself, which expects FALSE or a value.
    blank_filter     = pmp_params$blank_filter %||% FALSE,
    replace_MVs      = pmp_params$replace_MVs  %||% FALSE,
    impute_method    = pmp_params$impute_method %||% "gaussian",
    # `impute_seed: ~` asks for the session's random stream, so only an
    # absent key gets the default.
    impute_seed      = if("impute_seed" %in% names(pmp_params))
                         pmp_params$impute_seed else 42,
    batch_correction = pmp_params$batch_correction,
    bc_qc_label      = pmp_params$bc_qc_label,
    bc_factor_name   = pmp_params$bc_factor_name,
    bc_header        = pmp_params$bc_header,
    processing_batch = pmp_params$processing_batch,
    matrix_id_col      = matrix_id_col,
    matrix_orientation = matrix_orientation,
    # Lets correctBatch() test whether batch is confounded with the study
    # factor, which is the assumption a non-QC reference rests on.
    bc_check_factor  = pmp_params$bc_check_factor %||%
                         (project_params$project$stats_report %||%
                          list())$global_summary_factor,
    blank_head       = pmp_params$blank_head,
    blank_name       = pmp_params$blank_name,
    # Detection filters. Both default off, so an existing config reproduces
    # byte-for-byte; a block that is present but has enabled: False is also
    # off, which is what lets the example configs document the shape.
    filter_features  = .filter_block(pmp_params$filter_features),
    filter_samples   = .filter_block(pmp_params$filter_samples),
    normalize        = pmp_params$normalize %||% FALSE,
    adjust_conc      = isTRUE(ac_pars$enabled),
    starting_vol_col = ac_pars$starting_vol_col %||% FALSE,
    sample_vol_col   = ac_pars$sample_vol_col   %||% "sample_volume_uL",
    cal_vol_col      = ac_pars$cal_vol_col      %||% "cal_vol_uL",
    sample_IS_col    = ac_pars$sample_IS_col    %||% "sample_IS_vol_uL",
    cal_IS_col       = ac_pars$cal_IS_col       %||% "cal_IS_vol_uL",
    name_col            = name_col,
    include_col         = include_col,
    include_value       = include_value,
    compound_col        = compound_col,
    processing_name_col = processing_name_col,
    report_col          = report_col,
    report_value        = report_value,
    comment_col         = comment_col
  )

  combined_datamatrices = combined_data[[1]]
  removed_features      = combined_data[[2]]
  measured              = combined_data$measured %||% combined_datamatrices
  imputed               = !isFALSE(pmp_params$replace_MVs %||% FALSE)

  # Design checks need the data and the config together, so they run here
  # rather than at load time. These are warnings by default: a thin group or a
  # batch without QCs is a judgement call, not a malformed input - except
  # complete confounding, which validateDesign() reports as an error.
  dz_check = tryCatch(
    validateDesign(combined_datamatrices, project_params,
                    sample_type_head = (project_params$project$data_quality_report %||%
                                        list())$sample_type_head %||% "Sample_type",
                    qc_label   = (project_params$project$data_quality_report %||%
                                  list())$qc_label %||% "QC",
                    blank_name = pmp_params$blank_name %||% "Blank",
                    batch_head = pmp_params$bc_header  %||% "Chrom_Batch"),
    error = function(e){
      warning("[runMRManalyzeR] design validation skipped: ",
              conditionMessage(e)); NULL
    })
  if(!is.null(dz_check)){
    for(w in dz_check$warnings) message("[validateDesign] ", w)
    for(e in dz_check$errors)   warning("[validateDesign] ", e)
  }

  scale_fac = pmp_params$scale_fac %||% 1
  combined_datamatrices$data = combined_datamatrices$data * scale_fac
  measured$data              = measured$data * scale_fac

  # Per-feature CV metrics (QC, sample, sample/QC ratio) -> variable_meta,
  # matched to features by Compound. Uses the data_quality_report QC/sample
  # labels (defaults if that block is absent).
  dq_cv = project_params$project$data_quality_report %||%
            project_params$project$UVA_report %||% list()
  # Computed on the measured values: a filled-in value says nothing about how
  # reproducibly a compound was measured.
  measured = addCVMetrics(
    measured,
    sample_type_head = dq_cv$sample_type_head %||% "Sample_type",
    qc_label         = dq_cv$qc_label         %||% "QC",
    sample_labels    = dq_cv$sample_labels    %||% "Sample")
  combined_datamatrices = .copy_cv_cols(combined_datamatrices, measured)

  add_info = data.frame(sheet = "matrix", info = "units",
                        value = pmp_params$units %||% datatype)

  .write_output_xlsx(out_xlsx, combined_datamatrices, removed_features,
                     add_info, measured = if(imputed) measured)
  saveRDS(combined_datamatrices, file = out_RDS)
  # With imputation on, the measured values are stored beside it, so a later
  # execute: False run can still give the reports the real values.
  out_meas = .measured_rds(out_RDS)
  if(imputed){
    saveRDS(measured, file = out_meas)
    message("Wrote: ", out_meas)
  } else if(file.exists(out_meas)){
    invisible(file.remove(out_meas))   # left over from a run that imputed
  }

  con = file(out_pars, open = "wt"); on.exit(close(con), add = TRUE)
  utils::capture.output(project_params, file = con)

  message("Wrote: ", out_xlsx)
  message("Wrote: ", out_RDS)

  list(combined_datamatrices = combined_datamatrices,
       removed_features      = removed_features,
       measured              = measured)
}

#' Where the measured dataset is stored when the output matrix is imputed
#' @keywords internal
#' @noRd
.measured_rds = function(out_RDS) sub("\\.RDS$", "_measured.RDS", out_RDS)

#' Swap each combine input for its measured copy where one exists
#' @param paths Input paths, possibly named by tag.
#' @return `paths` with `<name>.RDS` replaced by `<name>_measured.RDS`
#'   wherever that file exists; names are kept.
#' @keywords internal
#' @noRd
.prefer_measured = function(paths){
  is_rds = grepl("\\.RDS$", paths, ignore.case = TRUE)
  meas   = sub("\\.RDS$", "_measured.RDS", paths, ignore.case = TRUE)
  use    = is_rds & file.exists(meas)
  if(any(use))
    message(sprintf("[combine] measured values stored for %d input(s): %s",
                    sum(use), paste(basename(meas[use]), collapse = ", ")))
  paths[use] = meas[use]
  paths
}

#' Copy the CV columns computed on the measured values onto another dataset
#' @keywords internal
#' @noRd
.copy_cv_cols = function(to, from){
  vt = as.data.frame(to$variable_meta)
  vf = as.data.frame(from$variable_meta)
  at = match(rownames(vt), rownames(vf))
  for(cc in intersect(c("CV_QC", "CV_sample", "CV_sample_vs_QC"),
                      colnames(vf)))
    vt[[cc]] = vf[[cc]][at]
  to$variable_meta = vt
  to
}

#' @keywords internal
#' @noRd
.write_output_xlsx = function(out_xlsx, combined_datamatrices,
                              removed_features, add_info, measured = NULL){
  wb = openxlsx::createWorkbook()

  openxlsx::addWorksheet(wb, "feature_metadata")
  openxlsx::writeData(wb, "feature_metadata",
                      combined_datamatrices$variable_meta,
                      colNames = TRUE, rowNames = FALSE, keepNA = FALSE)

  openxlsx::addWorksheet(wb, "sample_metadata")
  openxlsx::writeData(wb, "sample_metadata",
                      combined_datamatrices$sample_meta,
                      colNames = TRUE, rowNames = FALSE, keepNA = FALSE)

  mat = combined_datamatrices$data %>%
    dplyr::select(dplyr::where(function(x) any(!is.na(x))))
  openxlsx::addWorksheet(wb, "matrix")
  openxlsx::writeData(wb, "matrix", mat,
                      colNames = TRUE, rowNames = TRUE, keepNA = FALSE)

  # With replace_MVs on, 'matrix' holds imputed values. These are the values
  # as measured, gaps included - what the reports use - on the same columns.
  if(!is.null(measured)){
    openxlsx::addWorksheet(wb, "matrix_measured")
    openxlsx::writeData(wb, "matrix_measured",
                        as.data.frame(measured$data)[, colnames(mat),
                                                     drop = FALSE],
                        colNames = TRUE, rowNames = TRUE, keepNA = FALSE)
  }

  openxlsx::addWorksheet(wb, "removed_features")
  openxlsx::writeData(wb, "removed_features", removed_features,
                      colNames = TRUE, rowNames = FALSE, keepNA = FALSE)

  openxlsx::addWorksheet(wb, "key")
  openxlsx::writeData(wb, "key", add_info,
                      colNames = TRUE, rowNames = FALSE, keepNA = FALSE)

  openxlsx::saveWorkbook(wb, file = out_xlsx, overwrite = TRUE)
}


#' @keywords internal
#' @noRd
.render_reports = function(project_params, project_paths, pmp_params,
                           dq_params, st_params,
                           out_RDS, out_stub,
                           combined_datamatrices, removed_features,
                           stats_tables = NULL,
                           measured_datamatrices = combined_datamatrices,
                           imputed = NULL){

  render_one = function(rmd_file, out_html, params_block){

    rmd_path = system.file("rmd", rmd_file, package = "MRManalyzeR")
    if (rmd_path == "") {
      rmd_path = file.path("inst", "rmd", rmd_file)
    }
    if(rmd_path == ""){
      warning(sprintf("Report template '%s' not found in installed package.",
                      rmd_file))
      return(invisible(NULL))
    }

    # Parent on the package namespace so the Rmd can find package
    # functions (runPCA, subsetDataset, ...) even when the caller
    # invoked us via `MRManalyzeR::runMRManalyzeR()` without
    # `library(MRManalyzeR)`. Lookup chain: chunk env -> render env ->
    # MRManalyzeR namespace -> imports -> base.
    pkg_ns = tryCatch(asNamespace("MRManalyzeR"),
                      error = function(e) globalenv())
    env = new.env(parent = pkg_ns)
    env$project_params        = project_params
    env$project_paths         = project_paths
    env$pmp_params            = pmp_params
    env$report_pars           = params_block
    env$include_MVA           = isTRUE(params_block$include_MVA) ||
                                 isTRUE(params_block$pca$enabled)
    env$out_RDS               = out_RDS
    # Measured values unless this report asks for the imputed matrix.
    use_imp = isTRUE(params_block$use_imputed)
    env$combined_datamatrices = if(use_imp) combined_datamatrices else
      measured_datamatrices
    # Whether those values include imputed ones, for the report to say so.
    # `imputed` gives it per matrix when the n_imputed record cannot: a
    # dataset stored without one.
    env$values_imputed = if(is.null(imputed))
      .n_imputed_total(env$combined_datamatrices) > 0 else
      isTRUE(imputed[[if(use_imp) "imputed" else "measured"]])
    env$removed_features      = removed_features
    env$stats_tables          = stats_tables

    message("Rendering: ", out_html)
    rmarkdown::render(
      input       = rmd_path,
      output_file = out_html,
      envir       = env,
      quiet       = TRUE
    )
  }

  if(isTRUE(dq_params$execute)){
    render_one("data_quality_report.Rmd",
               paste0(out_stub, "_data_quality_report.html"),
               dq_params)
  }
  if(isTRUE(st_params$execute)){
    .check_feature_cols(st_params, combined_datamatrices)
    render_one("stats_report.Rmd",
               paste0(out_stub, "_stats_report.html"),
               st_params)
  }
}

#' Warn about stats_report feature columns the dataset does not have
#'
#' A grouping or colouring column that does not resolve is otherwise dropped
#' without a word: the loadings come out uncoloured and the heatmap
#' ungrouped. Only sections that will render are checked.
#' @param st_params The `stats_report:` block.
#' @param de The dataset the report is drawn from.
#' @return `NULL`, invisibly; warns once per missing column.
#' @keywords internal
#' @noRd
.check_feature_cols = function(st_params, de){
  vm  = as.data.frame(de$variable_meta)
  on  = function(x, default) isTRUE(x %||% default)
  pca = st_params$pca %||% list()
  hm  = st_params$heatmap %||% list()
  ld  = pca$loadings %||% list()

  # Correlation heatmaps group by heatmap.group_features_by too, whether or
  # not the heatmap section itself is on.
  corr_on = length(tryCatch(.section_entries(st_params$correlations,
                                             "correlations"),
                            error = function(e) list())) > 0

  want = list()
  if(on(pca$enabled, TRUE) && on(ld$enabled, TRUE))
    want[["pca: loadings: color_by"]] = ld$color_by
  if(on(hm$enabled, TRUE) || corr_on)
    want[["heatmap: group_features_by"]] = hm$group_features_by
  if(on(hm$enabled, TRUE)){
    if(on((hm$general %||% list())$enabled, TRUE))
      want[["heatmap: general: group_features_by"]] =
        hm$general$group_features_by
    if(on((hm$top_changing %||% list())$enabled, FALSE))
      want[["heatmap: top_changing: group_features_by"]] =
        hm$top_changing$group_features_by
    if(on((hm$split_by_class %||% list())$enabled, FALSE))
      want[["heatmap: split_by_class: subclass_col"]] =
        hm$split_by_class$subclass_col
  }

  for(k in names(want)){
    col = want[[k]]
    if(!is.null(col) && is.null(.resolve_meta_col(col, vm)))
      warning(sprintf(
        "[stats_report] %s '%s' is not a feature_metadata column, so that grouping is not applied. Columns: %s.",
        k, col, paste(colnames(vm), collapse = ", ")), call. = FALSE)
  }
  invisible(NULL)
}


#' Run the combine-mode workflow from a dedicated combine YAML
#'
#' Merges several previously-saved datasets (`.RDS` or `.xlsx` outputs of
#' [`runMRManalyzeR()`]) into one `DatasetExperiment`, then runs the
#' `stats_report` analyses on the merged data. Skips PeakMatrixProcessing
#' and the data_quality_report.
#'
#' Defaults: samples are matched by `Sample_ID` (intersect across inputs);
#' `feature_meta_cols` and `sample_meta_cols` default to the intersection
#' of column names across inputs (specify them only to *narrow* the set).
#'
#' Outputs:
#' \itemize{
#'   \item `<output_stub>.RDS`        - merged `DatasetExperiment`
#'   \item `<output_stub>.xlsx`       - feature_metadata / sample_metadata / matrix tabs
#'   \item `<output_stub>_measured.RDS` - the same merge of the panels'
#'     measured values, when any panel stored a `_measured.RDS`; the xlsx
#'     then gains a `matrix_measured` tab
#'   \item `<output_stub>_stats.xlsx` - key / stats / correlations / linear_models
#'   \item `<output_stub>_stats_report.html`
#' }
#'
#' A panel run with `replace_MVs:` exports an imputed matrix and stores its
#' measured values beside it as `<name>_measured.RDS`. The merged `.RDS` and
#' `matrix` tab hold what each panel exported; the statistics and the report
#' use the merge of the measured values unless `stats_report` sets
#' `use_imputed: True`.
#'
#' YAML schema (see `inst/extdata/example_combine_config.yml`):
#' \preformatted{
#' combine:
#'   datasets:
#'     - path: ".../GOM_..._.RDS"      # .RDS or .xlsx, mix is fine
#'       tag:  "GOM"
#'       qc_remap: { "QC1": "QC_pooled" }
#'     - path: ".../cysLT_..._.xlsx"
#'       tag:  "cysLT"
#'   feature_meta_cols: ~              # NULL = intersect across inputs
#'   sample_meta_cols:  ~              # NULL = intersect across inputs
#'   prefix_features:   True
#'   sample_id_col:     Sample_ID
#'   output_stub:       ".../Combined/combined_HDM"
#'   units:             "combined"
#' stats_report:
#'   execute: True
#'   use_imputed: False             # True = statistics on the imputed merge
#'   comparisons:   [...]
#'   correlations:  [...]
#'   linear_models: [...]
#' }
#'
#' @param path_yaml Full path to the combine YAML.
#' @return Invisibly, a list with the merged `DatasetExperiment`, its
#'   `measured` counterpart (the same object when no panel stored one), the
#'   `stats_tables`, and the output paths.
#' @examples
#' # Two tiny processed panels (normally .RDS/.xlsx outputs of
#' # runMRManalyzeR()) that share sample IDs, merged via a combine YAML.
#' # stats_report execute = FALSE keeps the example to the merge itself.
#' dir <- tempfile("combine_"); dir.create(dir)
#' mk <- function(feats) struct::DatasetExperiment(
#'   data          = as.data.frame(matrix(1, 3, length(feats),
#'                     dimnames = list(c("S1", "S2", "S3"), feats))),
#'   sample_meta   = data.frame(Sample_ID = c("S1", "S2", "S3"),
#'                     Sample_type = "Sample",
#'                     row.names = c("S1", "S2", "S3")),
#'   variable_meta = data.frame(Compound = feats, Class = "lipid",
#'                     row.names = feats))
#' p1 <- file.path(dir, "panelA.RDS"); saveRDS(mk(c("A", "B")), p1)
#' p2 <- file.path(dir, "panelB.RDS"); saveRDS(mk(c("C", "D")), p2)
#' cfg <- list(
#'   combine = list(
#'     datasets        = list(list(path = p1, tag = "A"),
#'                            list(path = p2, tag = "B")),
#'     prefix_features = TRUE,
#'     sample_id_col   = "Sample_ID",
#'     output_stub     = file.path(dir, "combined")),
#'   stats_report = list(execute = FALSE))
#' yml <- file.path(dir, "combine.yml"); yaml::write_yaml(cfg, yml)
#' runMRManalyzeRCombine(yml)
#' @family entry points
#' @export
runMRManalyzeRCombine = function(path_yaml){

  stopifnot(file.exists(path_yaml))
  cfg = loadConfig(path_yaml)

  # Tolerate either a flat top-level layout or a `project:` wrapper.
  root = cfg$project %||% cfg
  combine_params = root$combine %||%
                     stop("[combine] YAML lacks a top-level `combine:` block.")
  st_params      = root$stats_report %||% list(execute = FALSE)

  ds = combine_params$datasets %||% list()
  if(!length(ds)) stop("[combine] `combine.datasets:` is empty.")

  paths = vapply(ds, function(d) d$path, character(1))
  tags  = vapply(ds, function(d) d$tag %||% NA_character_, character(1))
  # Named by tag, or after the file when any tag is missing - as
  # combineDatasets() would name them - so a panel keeps its name, and its
  # feature prefix, when it is read from its _measured.RDS.
  names(paths) = if(all(!is.na(tags) & nzchar(tags))) tags else
    tools::file_path_sans_ext(basename(paths))
  # A panel run with replace_MVs stores its measured values beside its main
  # RDS. The merged matrix is built from what each panel exported, and a
  # measured merge from those _measured.RDS files.
  paths_meas = .prefer_measured(paths)
  has_meas   = any(paths_meas != paths)

  # Per-dataset config maps, keyed by tag and by both paths, so the lookup
  # in combineDatasets() finds them for either merge.
  build_map = function(field){
    out = list()
    for(i in seq_along(ds)){
      v = ds[[i]][[field]]
      if(is.null(v)) next
      keys = c(if(!is.na(tags[i]) && nzchar(tags[i])) tags[[i]],
               paths[[i]], paths_meas[[i]])
      for(k in unique(keys)) out[[k]] = v
    }
    if(!length(out)) NULL else out
  }
  qc_remap_list  = build_map("qc_remap")
  rename_list    = build_map("feature_meta_rename")

  output_stub = combine_params$output_stub %||%
                stop("[combine] `combine.output_stub:` is required.")
  dir.create(dirname(output_stub), recursive = TRUE, showWarnings = FALSE)

  sid_col = combine_params$sample_id_col %||% "Sample_ID"
  combine_from = function(p) combineDatasets(
    paths               = p,
    feature_meta_cols   = combine_params$feature_meta_cols,
    sample_meta_cols    = combine_params$sample_meta_cols,
    feature_meta_rename = rename_list,
    qc_remap            = qc_remap_list,
    sample_id_col       = sid_col,
    drop_samples        = combine_params$drop_samples,
    prefix_features     = combine_params$prefix_features %||% FALSE,
    combined_name       = combine_params$combined_name   %||% basename(output_stub),
    duplicate_samples   = combine_params$duplicate_samples %||% "error"
  )

  message(sprintf("[combine] Merging %d dataset(s) ...", length(ds)))
  combined = combine_from(paths)
  message(sprintf("[combine] Result: %d samples x %d features.",
                  nrow(combined$data), ncol(combined$data)))
  combined_meas = combined
  if(has_meas){
    message("[combine] Merging the measured values ...")
    combined_meas = combine_from(paths_meas)
  }

  out_xlsx       = paste0(output_stub, ".xlsx")
  out_stats_xlsx = paste0(output_stub, "_stats.xlsx")
  out_RDS        = paste0(output_stub, ".RDS")
  out_pars       = paste0(output_stub, ".txt")

  add_info = data.frame(
    sheet = c("matrix", "matrix"),
    info  = c("source", "n_input_datasets"),
    value = c("combineDatasets()", as.character(length(ds)))
  )
  .write_output_xlsx(out_xlsx, combined,
                     removed_features = data.frame(),
                     add_info          = add_info,
                     measured          = if(has_meas) combined_meas)
  saveRDS(combined, file = out_RDS)
  out_meas_RDS = .measured_rds(out_RDS)
  if(has_meas){
    saveRDS(combined_meas, file = out_meas_RDS)
  } else if(file.exists(out_meas_RDS)){
    invisible(file.remove(out_meas_RDS))   # left over from an earlier merge
  }
  con = file(out_pars, open = "wt"); on.exit(close(con), add = TRUE)
  utils::capture.output(cfg, file = con)
  message("Wrote: ", out_xlsx)
  message("Wrote: ", out_RDS)
  if(has_meas) message("Wrote: ", out_meas_RDS)

  # The statistics and report use the measured merge unless the report asks
  # for the imputed one.
  use_imp = isTRUE(st_params$use_imputed)
  st_data = if(use_imp) combined else combined_meas

  stats_tables = NULL
  if(isTRUE(st_params$execute)){
    has_stats =
      length(.section_entries(st_params$comparisons,   "comparisons"))   > 0 ||
      length(.section_entries(st_params$correlations,  "correlations"))  > 0 ||
      length(.section_entries(st_params$linear_models, "linear_models")) > 0 ||
      length(.ion_entries(st_params$ion_ratios))                          > 0
    if(has_stats){
      message("[combine] Running statistics on merged data ...")
      stats_tables = runStats(st_data, st_params)
      append_stats_xlsx(out_stats_xlsx, stats_tables)
      message("Wrote stats to: ", out_stats_xlsx)
    }
  }

  if(isTRUE(st_params$execute)){
    .check_feature_cols(st_params, st_data)
    rmd_path = system.file("rmd", "stats_report.Rmd", package = "MRManalyzeR")
    if(rmd_path == "")
      rmd_path = file.path("inst", "rmd", "stats_report.Rmd")
    pkg_ns = tryCatch(asNamespace("MRManalyzeR"),
                      error = function(e) globalenv())
    env = new.env(parent = pkg_ns)
    units_lbl = combine_params$units %||% "combined"
    env$project_params        = cfg
    env$project_paths         = list(datatype = units_lbl)
    env$pmp_params            = list(units = units_lbl, datatype = units_lbl)
    env$report_pars           = st_params
    env$out_RDS               = out_RDS
    env$combined_datamatrices = st_data
    env$values_imputed        = .inputs_imputed(
      if(use_imp) paths else paths_meas, sid_col)
    env$removed_features      = data.frame()
    env$stats_tables          = stats_tables

    out_html = paste0(output_stub, "_stats_report.html")
    message("Rendering: ", out_html)
    rmarkdown::render(input = rmd_path, output_file = out_html,
                      envir = env, quiet = TRUE)
  }

  invisible(list(
    datasetExperiment = combined,
    measured          = combined_meas,
    stats_tables      = stats_tables,
    out_xlsx          = out_xlsx,
    out_stats_xlsx    = out_stats_xlsx,
    out_RDS           = out_RDS,
    out_pars          = out_pars,
    n_input_datasets  = length(ds)
  ))
}
