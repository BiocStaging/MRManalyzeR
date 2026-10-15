#' Write a configuration file to start from
#'
#' A run is driven by one YAML file, so the first thing anyone needs is a
#' file to edit. This writes one: `minimal = TRUE` gives the keys a run
#' cannot do without, each on one line with its options beside it;
#' `minimal = FALSE` copies the fully commented reference configuration
#' shipped with the package, which documents every key the package reads.
#'
#' Nothing is guessed from your data. The paths, the sheet names and the
#' metadata columns in the file are placeholders, and the run stops with a
#' message naming whichever one does not match your workbook - so editing is
#' a matter of answering those messages, in order.
#'
#' @param path File to write. A directory is also accepted, and the file is
#'   then called `project_config.yml` inside it.
#' @param template `"study"` for a single acquisition, the usual case, or
#'   `"combine"` for merging several processed panels by sample.
#' @param source Data source for a `"study"` template: `"targetlynx"`,
#'   `"skyline"` or `"matrix"` (any sample x analyte table). Ignored for
#'   `"combine"`.
#' @param minimal `TRUE` (default) writes the short form. `FALSE` copies the
#'   bundled reference configuration, comments and all.
#' @param overwrite Overwrite `path` if it exists. `FALSE` (default) stops
#'   rather than replace a configuration you have already edited.
#' @return The path written, invisibly.
#' @family entry points
#' @seealso [runMRManalyzeR()], which runs the file this writes.
#' @examples
#' cfg <- newConfig(tempfile(fileext = ".yml"))
#' cat(readLines(cfg)[1:12], sep = "\n")
#' @export
newConfig <- function(path,
                     template  = c("study", "combine"),
                     source    = c("targetlynx", "skyline", "matrix"),
                     minimal   = TRUE,
                     overwrite = FALSE){

  template <- match.arg(template)
  source  <- match.arg(source)

  if(!is.character(path) || length(path) != 1L || !nzchar(path))
    stop("[newConfig] `path` must be a single file or directory name.",
         call. = FALSE)
  if(dir.exists(path))
    path <- file.path(path, "project_config.yml")
  if(file.exists(path) && !isTRUE(overwrite))
    stop(sprintf(
      "[newConfig] %s exists. Pass overwrite = TRUE to replace it.", path),
      call. = FALSE)

  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)

  if(!isTRUE(minimal)){
    ref <- system.file("extdata",
                      if(identical(template, "combine"))
                        "example_combine_config.yml" else "example_config.yml",
                      package = "MRManalyzeR")
    if(!nzchar(ref))
      stop("[newConfig] the bundled reference configuration is not installed.",
           call. = FALSE)
    file.copy(ref, path, overwrite = TRUE)
  } else {
    writeLines(if(identical(template, "combine")) .cfg_combine()
               else .cfg_study(source), path)
  }

  message("Wrote: ", path,
          "\nEdit paths:, then the column names under PeakMatrixProcessing:,",
          " then run:\n  runMRManalyzeR(\"", path, "\")")
  invisible(path)
}

#' The short single-acquisition configuration
#' @param source One of targetlynx, skyline, matrix.
#' @return Character vector of lines.
#' @keywords internal
#' @noRd
.cfg_study <- function(source){
  data_block <- switch(
    source,
    targetlynx = c(
      '    tl_data:',
      '      enabled: True',
      '      data_tab_names: [lcms_data]     # sheet(s) holding the export',
      '      tl_headers: ["ID", "Name", "Area", "ng/mL", "Response", "S/N"]',
      '      datatype: ["Area"]              # which of those to report',
      '      snr: False                      # number | False'),
    skyline = c(
      '    skyline_data:',
      '      enabled: True',
      '      data_tab_names: [skyline_data]  # molecule x sample sheet(s)',
      '      signal_filter: ["LOD"]          # LOD | LOQ | False'),
    matrix = c(
      '    matrix_data:',
      '      enabled: True',
      '      data_tab_names: [matrix_data]   # sheet holding the matrix',
      '      datatype: ["Area"]              # output label',
      '      orientation: samples_rows       # samples_rows | samples_cols',
      '      signal_filter: False            # LOD | LOQ | False'))

  c(
'# MRManalyzeR configuration. Every `enabled:` and `execute:` is True | False;',
'# a key with a list of values after `#` takes one of those; anything else is',
'# a column name in your metadata sheets, a number, or free text.',
'#',
'# The fully commented reference, documenting every key the package reads:',
'#   newConfig("full_config.yml", minimal = FALSE)',
'',
'project:',
'  name: my_study',
'  researcher: ""',
'',
'  # Where the workbook is, and where results go. fn is the workbook name',
'  # inside data_dir, without the .xlsx.',
'  paths:',
'    project_dir: &project_dir "C:/path/to/study"',
'    data_dir:    &data_dir    !concat_path [*project_dir, "data"]',
'    result_dir:  &result_dir  !concat_path [*project_dir, "results"]',
'    fn: &fn my_workbook',
'    suffix: _run1                 # appended to every output name',
'',
'  PeakMatrixProcessing:',
'    execute: True                 # False = reuse the stored RDS, re-render only',
'',
'    # 1. Data source',
    data_block,
'',
'    # 2. Sample metadata: which injections to process',
'    sample_meta_tab: sample_metadata',
'    name_col: Name                # injection id, matches the export',
'    include_col: Include',
'    include_value: "YES"',
'',
'    # 3. Feature metadata: which compounds to report',
'    feature_meta_tab: feature_metadata',
'    compound_col: Compound',
'    processing_name_col: Processing_name',
'    report_col: Report',
'    report_value: "YES"',
'    comment_col: Comment',
'',
'    # 4. Processing. Each step is off until you set it.',
'    blank_filter: False           # number (x blank mean) | False',
'    blank_head: Sample_type',
'    blank_name: Blank',
'    filter_features:',
'      enabled: False',
'      min_frac: 0.5               # 0-1, within the better-detected group',
'      method: within              # across | within | QC',
'      group_col: Treatment',
'    batch_correction: False       # True | False | feature_metadata YES/NO column',
'    bc_qc_label: "QC"',
'    bc_factor_name: "Sample_type"',
'    bc_header: "Chrom_Batch"',
'    normalize: False              # False | sample_metadata column to divide by',
'    replace_MVs: False            # False | fraction of the feature minimum',
'    impute_method: gaussian       # gaussian | constant',
'    impute_seed: 42',
'',
'  # ====================== REPORTS ======================',
'  data_quality_report:',
'    execute: True',
'    use_imputed: False            # False = CVs and plots on measured values',
'    sample_type_head: Sample_type',
'    qc_label: QC',
'    sample_labels: ["Sample"]',
'    injection_order_head: Injection_order',
'    sample_id_head: Sample_ID',
'    qc_pca:',
'      enabled: True',
'      color_by: [Injection_order, Sample_type]',
'',
'  stats_report:',
'    execute: True',
'    use_imputed: False',
'    sample_type_head: Sample_type',
'    sample_labels: ["Sample"]',
'    sample_id_head: Sample_ID',
'    global_summary_factor: "Treatment"',
'    p_adjust: BH                  # BH | bonferroni',
'    sig_threshold: 0.05',
'',
'    comparisons:',
'      enabled: True',
'      entries:',
'        - name: A_vs_B',
'          method: welch           # welch | student | wilcoxon | anova | kruskal',
'          compare: { factor: Treatment, levels: [A, B] }',
'',
'    boxplot:',
'      enabled: True',
'      factors: [Treatment]',
'',
'    pca:',
'      enabled: True',
'      components: [1, 2]',
'',
'    heatmap:',
'      enabled: True',
'      style: auto                 # auto | samples | group_means',
'      color_samples_by: Treatment',
'      group_features_by: ~        # feature_metadata column, or ~ for one block',
'      transform: log2             # none | log2 | sqrt',
'      scale: auto                 # auto | z | feature')
}

#' The short combine configuration
#' @return Character vector of lines.
#' @keywords internal
#' @noRd
.cfg_combine <- function(){
  c(
'# MRManalyzeR combine configuration: merge processed panels by sample, then',
'# run the statistics across them. Each input is the .RDS or .xlsx a normal',
'# run wrote. The fully commented reference:',
'#   newConfig("full_combine.yml", template = "combine", minimal = FALSE)',
'',
'combine:',
'  datasets:',
'    - path: "C:/path/to/results/panel_A_Area_.RDS"',
'      tag:  "A"                   # prefixes this panel\'s compounds',
'    - path: "C:/path/to/results/panel_B_Area_.RDS"',
'      tag:  "B"',
'  feature_meta_cols: ~            # ~ = the columns the panels share',
'  sample_meta_cols:  ~',
'  prefix_features:   True         # True | False | auto (only on collisions)',
'  sample_id_col:     Sample_ID    # what the panels are matched on',
'  output_stub:       "C:/path/to/results/combined/my_study"',
'  units:             "combined"',
'',
'stats_report:',
'  execute: True',
'  use_imputed: False              # False = each panel\'s measured values',
'  sample_type_head: Sample_type',
'  sample_labels: ["Sample"]',
'  global_summary_factor: "Treatment"',
'  p_adjust: BH',
'  sig_threshold: 0.05',
'',
'  comparisons:',
'    enabled: True',
'    entries:',
'      - name: A_vs_B',
'        method: welch             # welch | student | wilcoxon | anova | kruskal',
'        compare: { factor: Treatment, levels: [A, B] }',
'',
'  pca:',
'    enabled: True',
'    components: [1, 2]')
}
