# MRManalyzeR 0.99.7

* The `structToolbox` example in the getting-started vignette is guarded on
  `pls` and `structToolbox` being installed, and both are declared in
  Suggests. `structToolbox` suggests `pls` rather than requiring it, so a
  builder that installs our Suggests but not theirs had the model function
  without the engine behind it, and the vignette failed to build.

* The rename note in 0.99.6 said five functions where it listed seven.

# MRManalyzeR 0.99.6

Changes from the Bioconductor review.

* Seven functions take an `mrm` prefix, because their old names were also
  exported by packages a reader might have attached: `runPCA()`,
  `plotPCA()` and `plotHeatmap()` become `mrmRunPCA()`, `mrmPlotPCA()` and
  `mrmPlotHeatmap()` (`scater`); `runExample()` becomes `mrmRunExample()`
  (`shiny`); and the three validators become `mrmValidateInput()`,
  `mrmValidateConfig()` and `mrmValidateDesign()` (`validateConfig()` is a
  name quantMSImageR uses as well). A script that loaded two of these got
  whichever was attached last. The remaining exports were checked against
  every package installed here and keep their names. Scripts calling the old
  names need updating; the YAML configurations do not, as they name no
  functions.

* `newConfig()` writes a configuration to edit, so a first run does not start
  from a blank file. `minimal = TRUE` gives the keys a run cannot do without,
  `FALSE` copies the fully commented reference; `template = "combine"` writes
  the merge-several-panels form.

* Every exported step checks that it was handed a `DatasetExperiment`, and
  says which function was called with what, instead of failing further in on
  a slot the object does not have.

* A step that writes its record into `variable_meta` or `sample_meta` now
  warns when it is replacing columns from an earlier run -
  `addCVMetrics()`, `correctBatch()`, `imputeMissing()`, `filterBlanks()`,
  `filterFeatures()` and `filterSamples()`. The second set of numbers
  describes different data from the first.

* The package uses base R where it used to carry a dependency: `|>` instead
  of `magrittr`'s `%>%`, and base R's `%||%` instead of its own. `magrittr`
  is no longer a dependency, and `stats` and `utils` are declared in
  DESCRIPTION, which they always should have been.

* `=` is no longer used for assignment in package code.

* `mrm_validation`, what the three validation functions return, is
  documented: its three elements and what `print()` does with them.

* `data-raw/` is now `inst/script/`, with a note recording where everything
  in `inst/extdata` came from - the real workbook's publication, and the
  script that generates the simulated one.

* The README gives Bioconductor installation instructions.

Also in this release.

* `plotGroupHeatmap()` draws a key for its opacity. The cell opacity carries
  meaning - whether a sample moved with the direction its group's mean moved
  in, or against it - and until now nothing on the figure said so.

* Both vignettes and the two report templates were revised after an
  independent read. The corrections of substance: the detection-filter
  example demonstrated the opposite of what its text claimed, the
  blank-filter example ran on a matrix whose blanks were already masked, the
  walkthrough switched to the imputed dataset where it said it was using
  measured values, `CV_sample_vs_QC` was described with a threshold that
  belongs to a different statistic, and the stats report's Methods text still
  described the pre-0.99.3 engine (two tests per comparison, adjustment
  across the whole workbook).

* `feature_metadata` in the results workbook leads with `Lowest_measured`
  after the compound name: the smallest value the instrument actually
  reported for that compound in the **biological samples**, read before
  imputation, and `NA` for a compound never measured there. QC and blank
  injections are excluded - a pooled QC is not study material, and a blank is
  what the signal is measured against. It is the floor to read a borderline
  concentration against, and is rounded to the nearest whole number - the
  column is there to be read, and a peak area means nothing at the eighth
  decimal place. It is an integer whatever the datatype, so on a run
  reporting concentrations a floor below 0.5 reads as zero.

* `removed_features` leads with `Failed_step`, naming what took each
  compound out - `User exclude`, `Detection filter`, `Signal filter (LOD)`,
  `Never measured`, `Not in source data` or `No metadata entry`. The earliest
  decision wins, so a compound the curator excluded says so whatever else is
  true of it.

* Compounds dropped by the detection filter are listed in `removed_features`.
  They used to leave the matrix without leaving a row behind, so a panel
  member could disappear from the report with nothing recorded.

# MRManalyzeR 0.99.5

* `DT` moves from Suggests to Imports. Both reports call `DT::datatable()`
  without a fallback, so a run without it installed stopped part-way through
  rendering.

* The group-mean heatmap takes a colour scale: `scale: z` keeps one colour
  limit for the whole panel, so compounds can be compared with each other and
  a pale row is one whose groups separated less; `scale: feature` gives each
  compound its own limit, so every row fills the colour range and only the
  shape along the groups can be read; `scale: auto`, the example config's
  setting, picks shared at two groups and per compound at three or more,
  where the question is usually which group is high rather than which
  compound moved most. `plotGroupHeatmap(group_scale = )` is the same choice,
  without `auto`. `scale:` used to be ignored by that style.

* New package logo, drawn top right in both reports and used by the README.

* `colours:` in either report block takes a `levels:` map that pins a class to
  a colour - `COX: "#B03A2E"` - wherever it is drawn: the heatmap annotation
  bar, the PCA loadings and the correlation blocks. Levels not listed keep the
  palette. Names are matched ignoring whitespace, so a heading a spreadsheet
  stored with non-breaking spaces (`COX || CYP`) can be written with ordinary
  ones.

* `plotGroupHeatmap()` draws the compound names between the class colour bar
  and the cells whenever `group_features_by` is set - block, then names, then
  data, which is how a pathway figure is read. They used to sit on the axis,
  outside the bar. Nothing to configure; without a class column they stay on
  the axis, and `feature_labels: hide` still drops them. `mrmPlotHeatmap()` is
  unchanged: pheatmap draws its row names at the edge and cannot move them.

# MRManalyzeR 0.99.4

* New processing order: blank filter, detection filters, batch correction,
  normalisation, concentration adjustment, and imputation last. Batch
  correction now works on the instrument signal - normalising first rescaled
  the QC reference by per-sample amounts that pooled QCs usually do not have -
  and on measured values only, since imputing first let invented values move a
  batch median. Imputation is still applied per `processing_batch`. Results can
  differ for runs with `batch_correction:` on together with `normalize:`,
  `adjust_conc:` or `replace_MVs:`, and, with the median change below, for any
  run whose reference injections have missing values.

* `correctBatch()` counts a missing value as zero when it takes a median, so
  `c(NA, NA, 1, 2, 3)` has median 1 rather than 2. A compound whose reference
  median is zero or negative in any batch, or across all reference injections,
  is now corrected in no batch at all - where a batch with no measured
  reference value used to turn it into `NA`. Each such compound is listed on
  screen with the batches that failed, and recorded in the new
  `bc_failed_batches` column of `variable_meta`, next to `batch_corrected`.

* `batch_correction:` accepts the name of a `feature_metadata` column as well
  as `True` / `False`. The column holds `YES` or `NO` per compound and only the
  `YES` compounds are corrected; a blank cell counts as `NO`, with a warning.
  The column is checked before the data sheets are read, under its heading
  as typed or as `assembleDataset()` renames it. The same selection is
  available directly as `correctBatch(feature_col = )`.

* The internal `batch_correct` struct model is removed and its arithmetic moved
  into `correctBatch()`, which already wrapped it. Its training and prediction
  methods were identical, so every correction was computed twice, and its
  `qc_frac` and `sample_frac` settings were never used.

* An injection with an empty batch label now stops `correctBatch()` with an
  error that names it, instead of the unclear "batch(es) 'NA' have fewer than
  2 reference samples". `validateDesign()` no longer counts an empty label as
  a batch.

* An entry switched off with its own `enabled: False` in `comparisons:`,
  `correlations:` or `linear_models:` is now skipped everywhere. `runStats()`
  used to run it anyway.

* The reason recorded for a feature dropped after reading now says whether it
  had no values at all in the included samples, or had every value removed by
  the signal filter. Both used to read "All values below <filter> threshold"
  (or "All-NA after processing" with no filter).

* Injection order must be numeric. `validateInput()` gains
  `injection_order_head` and `include_value`, and `runMRManalyzeR()` stops
  before processing when an included injection's order is missing or not a
  number - checked when the data-quality report will be drawn, which is the
  only thing that reads the column. `plotDrift()` stops rather than drawing an
  empty plot.

* A `stats_report:` feature column the data does not have
  (`pca.loadings.color_by`, `heatmap.group_features_by` and the heatmap
  sub-block equivalents) now raises a warning, and the report says the
  grouping was not applied instead of describing it. A heading renamed on
  import (`S-group` stored as `S.group`) is now matched in every heatmap
  section, as it already was for the loadings and the general heatmap.

* Missing values are imputed by random draws below each compound's lowest
  measurement - a normal centred at half the ceiling (fraction x minimum), SD a
  sixth of it, truncated to lie between 0 and the ceiling - rather than all
  set to one constant, which shrank each compound's variance and tied its
  ranks. `imputeMissing()` gains `method = "gaussian"` (default) or
  `"constant"` (the old fill) and `seed` (default 42, so a rerun is
  identical); `mrmRunPCA()` gains `impute = "gaussian"`, now its default, and
  `seed`. In YAML: `impute_method:` and `impute_seed:` under
  `PeakMatrixProcessing:`, and `seed:` in the PCA `preprocess:` blocks. A
  report PCA whose `preprocess:` omits `impute:` now uses `gaussian` rather
  than `min`, and `mrmRunPCA()` counts zeros as missing, as `imputeMissing()`
  already did.

* The reports use measured values. With `replace_MVs:` set, the xlsx
  `matrix` sheet and the RDS still hold the imputed matrix, and a
  `matrix_measured` sheet and `<name>_measured.RDS` are added beside them. The
  data-quality report and the statistics (tests, summaries, boxplots, volcano)
  use the measured values and CVs are computed on them. The PCAs use the
  values `replace_MVs:` imputed, the same as the xlsx matrix, and impute for
  themselves only what is still missing - everything, when `replace_MVs:` is
  off. `use_imputed: True` in either report block switches that report's
  tests and plots to the imputed matrix. `processDataset()` returns the measured dataset as
  `measured`. Combine mode does the same: the merged matrix is built from
  what each panel exported, and where panels stored a `_measured.RDS` a
  merge of those is written beside it (`matrix_measured`,
  `<output_stub>_measured.RDS`) and used for the statistics and report.
  Each report says which it was given, judged from the data rather than the
  setting. A dataset stored before this version with no imputation record
  (multi-batch runs kept none) is taken as imputed when `replace_MVs:` is
  set, and a run with `execute: False` on it warns that its reports will use
  imputed values.

* A zero in the data is read as a non-detect and set to `NA` straight after
  reading, before any step sees it. TargetLynx reports an area of 0 when it
  finds no peak, and that 0 used to count as a measurement: the detection
  filters took it as detected, and in the measured matrix it pulled group
  means, CVs and tests down. It is now a gap wherever the matrix is not
  imputed, as imputation and the PCA already treated it. The blank filter
  leaves such a blank out of its mean, as it does a value the signal filter
  masked, and a compound that is zero in every included sample is dropped as
  "No values in the included samples".

* On measured data, gaps show as gaps: the heatmap draws them grey instead of
  as the compound's mean, and a correlation resting on fewer than three shared
  samples is left blank instead of drawn as a perfect +/-1.

* `imputeMissing()` splits by processing batch itself (`by =`), so the
  per-compound record survives multi-batch runs: `n_imputed`, `frac_imputed`,
  and `impute_max` - renamed from `impute_fill`, since random fills have no
  single value - with one `impute_max_<batch>` column per batch when there are
  several.

* A `tl_data:` block without `data_tab_names:` now reads every sheet whose
  name contains `lcms_data`, as the example configuration documents.
  `runMRManalyzeR()` substituted `"skyline_data"` whatever the data source, so
  such a config stopped with "sheet(s) not found in workbook: skyline_data".
  `processDataset()`'s `data_tab_names` now defaults to `NULL` and resolves to
  each reader's usual sheets.

* A `PeakMatrixProcessing:` block that leaves out `normalize:`,
  `blank_filter:` or `replace_MVs:` now treats that step as off. The absent key
  was passed on as `NULL`, and an absent `normalize:` stopped the run inside
  `normaliseMatrix()` with "argument is of length zero". This can change what a
  config that omits `replace_MVs:` writes: the imputation step it used to enter
  with `NULL` imputed nothing, but still turned every zero in the matrix into
  `NA`, and added an all-`NA` `impute_fill` column plus all-zero `n_imputed`
  and `frac_imputed` columns to `variable_meta`. An omitted `blank_filter:` no
  longer adds `blank_threshold` and `n_masked_blank` either. Set the keys
  explicitly to keep either step.

* `validateConfig()` and `validateDesign()` now check a stats report's
  comparisons only when the report will run, which means `execute: True` - the
  same condition `runMRManalyzeR()` renders it under. A report that is switched
  off, or that omits `execute:`, can no longer stop a processing run or warn
  about levels it will never test.

* A metabolite ratio whose numerator or denominator is not in the data now
  raises a warning naming the missing feature, instead of being skipped
  silently - typically a combine prefix that does not match.

* `runMRManalyzeRCombine()` runs the statistics when `ion_ratios:` is the only
  section with entries. Its gate tested the other three sections only, so a
  combine config whose stats block held nothing but ratios produced no stats at
  all.

* `combineDatasets()` names any input file it cannot find.

* Reading a TargetLynx export no longer raises a tidyselect deprecation
  warning about an external vector in a selection.

* In combine mode, `qc_remap:` and `feature_meta_rename:` on a tagged panel
  now apply when another panel has no tag. They were looked up under the tag,
  which `combineDatasets()` only uses when every panel has one.

* A batch or processing-batch heading can be given in the config either as
  written in the workbook or as stored after import (`Chrom-Batch` or
  `Chrom.Batch`). The checks made before import accepted only the first, and
  `validateDesign()` only the second, warning that the column was missing.
  `filterBlanks()` and `imputeMissing()` find `blank_head` (and
  `imputeMissing(by =)`) the same way: a `Sample-type` heading stopped
  `filterBlanks()`, and `imputeMissing()` then filled the blanks as samples.

# MRManalyzeR 0.99.3

* **Function names are now camelCase**, following the Bioconductor style
  guide (`runMRManalyzeR()`, `processDataset()`, `mrmPlotHeatmap()`, and so on).
  This renames every exported function; there are no back-compatible aliases,
  so calling code must be updated. Acronyms stay upper case - `mrmPlotPCA()`,
  `addCVMetrics()`, `readTargetLynx()`.

* New `filterFeatures()`: drops features by detection frequency. `method =
  "within"` computes the detection fraction separately within each level of a
  grouping factor and keeps a feature that clears the threshold in any one of
  them, so a compound present in only one treatment survives; `"across"` pools
  every study sample and `"QC"` uses the QC injections. QC and blank
  injections never count towards the fraction, though they stay in the
  dataset. The fraction used is recorded per feature as `detected_frac`.

* New `filterSamples()`: drops study injections whose missing-value fraction
  exceeds a threshold, recording `na_frac` for every row. Only study samples
  are eligible - a QC over the threshold is reported but kept, since removing
  one can leave a batch with nothing for `correctBatch()` to anchor to.

* Both are off by default and are wired into `processDataset()` through new
  `filter_features:` and `filter_samples:` YAML blocks. They run once on the
  whole dataset, after blank masking and before normalisation and imputation:
  a detection fraction is a property of the study rather than of a batch, and
  imputation destroys the missingness both filters measure. Features are
  filtered before samples, so an injection is judged against features that
  are actually detectable.

* `processDataset()`'s per-batch helper was split accordingly, into a
  blank-filter pass and a normalise / concentration / impute pass, with the
  global filters in between. Per-batch imputation behaviour is unchanged.

# MRManalyzeR 0.99.2

* New `plotGroupHeatmap()`: a group-mean heatmap where the hue is the group
  mean and the opacity of each sample's sub-cell is that sample's contribution
  to it, so a mean resting on a single sample is visible rather than hidden.
  `mrmPlotHeatmap()` is unchanged and remains the right choice while the sample
  axis is still legible - it never averages.
* The colour limit defaults from the group means themselves rather than from a
  per-sample scale, since averaging replicates shrinks a scaled value by
  roughly sqrt(n) and an inherited limit leaves the panel washed out.

# MRManalyzeR 0.99.1

* Vignette and README reworked following review: repositioned around flexible
  quantitative outputs (raw areas through to calibrated concentrations) rather
  than internal-standard quantification alone, and around defined target lists
  rather than hypothesis testing alone.
* Vignette retitled, and validation promoted to its own stage in the workflow
  table and diagram.
* Corrected the description of `runStats()`: the test is named explicitly per
  comparison, not inferred from the design.
* Blanks are now excluded from the PCA examples, and the imputation fraction
  lowered to 0.2.
* "Ion ratios" renamed to "metabolite ratios" to avoid confusion with
  qualifier/quantifier transition ratios.
* Added Antonio Checa as an author, and a package logo.

# MRManalyzeR 0.99.0

* Initial Bioconductor submission.
* Post-acquisition processing of targeted LC-MS/MS lipidomics and metabolomics
  results exported from Waters TargetLynx or Skyline, or supplied as a plain
  sample-by-analyte matrix: peak-matrix assembly, S/N and LOD/LOQ filtering,
  blank filtering, missing-value imputation, normalisation, internal-standard
  and volume adjustment of vendor-reported concentrations, and batch
  correction. Chromatographic peak detection, integration and
  calibration-curve fitting are outside its scope.
* Every step is an exported function operating on a
  `struct::DatasetExperiment`; a whole workflow can additionally be driven
  from a single YAML config via `runMRManalyzeR()`.
* Self-contained HTML data-quality and statistical reports: group summaries,
  PCA, sample x feature heatmaps, volcano plots, feature-feature correlations,
  linear models, and enzyme-activity ion ratios.
* `runMRManalyzeRCombine()` merges multiple acquisition panels into a single
  analysis.
* `runExample()` runs the full workflow on the bundled example dataset
  (a subset of Kolmert et al. 2018, doi:10.1016/j.prostaglandins.2018.05.005).
* Every processing step is also an exported function operating on a
  `struct::DatasetExperiment`: `readTargetLynx()` / `readSkyline()`,
  `assembleDataset()`, `filterBlanks()`, `normaliseMatrix()`,
  `adjustConcentration()`, `imputeMissing()` and `correctBatch()`, composed
  by `processDataset()`.
* Per-compound quality metrics (`CV_QC`, `CV_sample` and their ratio) are
  written into the dataset's `variable_meta`.
* A second bundled dataset, `example_synthetic.xlsx` (plus the processed
  `example_synthetic.RDS`), is entirely simulated from a fixed seed: two
  chromatographic batches, per-sample protein amounts and reconstitution
  volumes, and a planted treatment effect on 5 of its 20 analytes, so
  normalisation, concentration adjustment, batch correction and the statistics
  can each be demonstrated against a known truth.
