#' MRManalyzeR: targeted LC-MS lipidomics and metabolomics processing
#'
#' Processes and analyses targeted lipidomics and metabolomics data exported from
#' Waters TargetLynx or Skyline: builds a peak-area / concentration matrix with
#' signal-to-noise or LOD/LOQ filtering, blank filtering, batch correction,
#' normalisation, calibration/internal-standard concentration adjustment and
#' missing-value imputation, then runs the statistical analyses and
#' renders two self-contained HTML reports.
#'
#' The workflow falls into four stages. Everything from stage 1 onwards is
#' carried in a single `struct::DatasetExperiment`, so each function takes that
#' object and returns it.
#'
#' @section Entry points:
#' [runMRManalyzeR()] drives all four stages from a single YAML config;
#' [runMRManalyzeRCombine()] merges several acquisition panels into one
#' analysis; [mrmRunExample()] runs the bundled example dataset end to end.
#'
#' @section 1. Data parse:
#' What did the instrument report? [readTargetLynx()] and [readSkyline()]
#' read a workbook into a wide sample x compound matrix, which
#' [assembleDataset()] turns into a `struct::DatasetExperiment`.
#' [loadConfig()], [loadDataset()] and [combineDatasets()] read a YAML
#' config, a stored dataset, and several datasets merged into one.
#'
#' @section 2. Peak-matrix processing:
#' What is the true concentration? [filterBlanks()], [correctBatch()],
#' [normaliseMatrix()], [adjustConcentration()] and [imputeMissing()], with
#' [subsetDataset()] to slice a dataset by sample metadata.
#' [processDataset()] composes stages 1 and 2 in one call.
#'
#' @section 3. QC check:
#' Can these numbers be trusted? [mrmRunPCA()], plus the per-compound CV metrics
#' added to `variable_meta` by [runMRManalyzeR()] and the bundled
#' data-quality report.
#'
#' @section 4. Stats:
#' What do these numbers mean? [runStats()] runs the comparisons,
#' correlations and linear models described by the YAML `stats_report:` block
#' and returns the same tables the statistics report and the output xlsx use.
#'
#' The `plot_*()` functions are exported rather than hidden inside the report
#' templates so that a figure can be reproduced, restyled or subset in one
#' line, without editing an R Markdown file or rebuilding it from the matrix.
#' The bundled reports call the same functions.
#'
#' See `vignette("MRManalyzeR")` for a step-by-step walkthrough.
#'
#' @keywords internal
#' @importFrom rlang .data
#' @importFrom DT datatable
#' @importFrom utils head
#' @importFrom stats median
#' @importFrom struct DatasetExperiment
"_PACKAGE"

# dplyr NSE column references used across the readers / assemble step.
utils::globalVariables(c("ID", "ID2", "Name", "S/N", "Report", "Compound"))

#' Resolve a configured metadata column name against the stored names
#'
#' `struct::DatasetExperiment` runs `make.names()` over the metadata it is
#' handed, so a column headed `S-group` in the source workbook is stored as
#' `S.group`, and `pubchem/KEGG_id` as `pubchem.KEGG_id`. A config naming the
#' original then matches nothing, and a plot that colours by it falls back to
#' no colouring at all - which is how a hyphen in a spreadsheet heading turns
#' into a missing legend with no error anywhere.
#'
#' Try the name exactly as written first, so a workbook that genuinely has
#' both `S-group` and `S.group` is not silently redirected, then its
#' `make.names()` form, and last a column whose own `make.names()` form
#' matches - so the stored spelling (`Chrom.Batch`) also finds a heading as
#' written in the workbook (`Chrom-Batch`), which the checks made before
#' import see.
#'
#' @param name Column name from the config, or `NULL`.
#' @param df Data frame whose columns are being matched.
#' @return The matching column name, or `NULL` if neither form is present.
#' @keywords internal
#' @noRd
.resolve_meta_col <- function(name, df){
  if(is.null(name)) return(NULL)
  nm <- as.character(name)[1]
  if(is.na(nm) || !nzchar(nm))  return(NULL)
  if(nm %in% colnames(df))      return(nm)
  alt <- make.names(nm)
  if(alt %in% colnames(df))     return(alt)
  hit <- colnames(df)[make.names(colnames(df)) == alt]
  if(length(hit) == 1L)         return(hit)
  NULL
}

#' Shared plot theme
#'
#' One theme for every plot the package draws, so the two HTML reports and
#' anything a user builds by calling the `plot_*()` functions directly look
#' like the same suite.
#'
#' @return A `ggplot2` theme object.
#' @keywords internal
#' @noRd
.mrm_theme <- function(){
  ggplot2::theme_classic() +
    ggplot2::theme(
      legend.position = "bottom",
      axis.title = ggplot2::element_text(face = "bold", color = "black",
                                         size = 12),
      axis.text  = ggplot2::element_text(face = "bold", color = "black",
                                         size = 8))
}

#' Vendored ltc colour palettes
#'
#' Hex values taken from the `ltc` palettes of Loukas Theodosiou
#' (https://github.com/loukesio/ltc-color-palettes). They are inlined rather
#' than depended on because `ltc` is distributed only on GitHub, and a
#' Bioconductor package cannot require a remote install.
#'
#' `hat` is the qualitative palette used for categorical annotations - feature
#' classes, sample types, batches - and `heatmap2` the diverging scheme used
#' for z-scores and correlations.
#'
#' @keywords internal
#' @noRd
.ltc <- list(
  # Reordered from the ltc original so the strong, maximally separable hues
  # come first - blue, red, green, yellow, purple - and the off-shades only
  # appear once a variable has more than five levels. A two- or three-level
  # factor should not be drawn in two shades of amber.
  #
  # The off-shades are ordered so each sits as far as possible from the strong
  # hue it resembles: teal, then crimson (5 places after red), light green (5
  # after green) and orange (5 after yellow), so no legend pairs yellow with
  # orange until it has nine levels.
  #
  # ltc's black is replaced by a mid grey and moved last. Black reads as a
  # rule or an axis rather than as a category, and it dominates a heatmap
  # annotation strip; grey recedes, which is what a ninth-most-important
  # level should do. It is darker than the grey70 used for NA in
  # .group_fill_colors(), so a real level is never confused with a missing one.
  hat = c("#4e54ac", "#e8351e", "#17a769", "#efb306", "#852f88",
          "#0f8096", "#cd023d", "#7db954", "#eb990c", "#5a5a5a"),
  reading = c("#EFBC68", "#919F89", "#EDBDAE", "#57717C",
              "#5F97A4", "#CAEAC8", "#95A1AE", "#C8CFD6"),
  heatmap0 = c("#001219", "#005F73", "#0A9396", "#94D2BD", "#E9D8A6",
               "#EE9B00", "#CA6702", "#AE2012", "#9B2226"),
  heatmap2 = c("#ca0020", "#f4a582", "#f7f7f7", "#92c5de", "#0571b0")
)

#' The palette currently selected for categorical variables
#'
#' Read from `getOption("MRManalyzeR.palette")` so a report can set it once
#' from YAML rather than every plotting function taking a palette argument.
#'
#' @return Character vector of hex colours.
#' @keywords internal
#' @noRd
.mrm_palette <- function(){
  nm <- getOption("MRManalyzeR.palette", "hat")
  .ltc[[nm]] %||% .ltc$hat
}

#' Diverging colour ramp for z-scores and correlations
#'
#' Reversed so it runs low-blue to high-red, the direction these figures are
#' read. The scheme is `getOption("MRManalyzeR.diverging")`, default
#' `heatmap2`.
#'
#' @param n Number of steps.
#' @return Character vector of `n` hex colours.
#' @keywords internal
#' @noRd
.diverging_pal <- function(n = 101){
  nm <- getOption("MRManalyzeR.diverging", "heatmap2")
  pal <- .ltc[[nm]] %||% .ltc$heatmap2
  grDevices::colorRampPalette(rev(pal))(n)
}

#' Qualitative palette for n categorical levels
#'
#' The vendored `ltc` `hat` palette, which carries ten distinguishable hues -
#' more than Brewer's Set1 - and is interpolated beyond that.
#'
#' @param n Number of levels.
#' @return Character vector of `n` colours.
#' @keywords internal
#' @noRd
.qual_pal <- function(n){
  if(n < 1) return(character(0))
  pal <- .mrm_palette()
  if(n <= length(pal)) return(pal[seq_len(n)])
  grDevices::colorRampPalette(pal)(n)
}

#' Named colour vector for the sorted levels of a categorical variable
#'
#' Keyed on sorted levels so a given level maps to the same colour wherever it
#' appears - the heatmap feature annotation and the PCA loadings bars must
#' agree, or the two figures cannot be read against each other. `NA` maps to
#' grey; a two-level variable uses a fixed blue/red pair, which reads as an
#' ordered contrast where an arbitrary pair would not.
#'
#' @param x Vector of level labels.
#' @return Named character vector of colours.
#' @keywords internal
#' @noRd
.group_fill_colors <- function(x){
  x <- as.character(x)
  lv <- sort(unique(x[!is.na(x) & x != "NA"]))
  m <- if(length(lv) == 0)      character(0)
       else if(length(lv) == 2) stats::setNames(
                                  .mrm_palette()[c(1L, 2L)], lv)
       else                     stats::setNames(.qual_pal(length(lv)), lv)
  if(any(is.na(x) | x == "NA")) m <- c(m, "NA" = "grey70")

  # A configured colour outranks the palette, so a class keeps one colour
  # across studies and matches figures drawn outside the package.
  fixed <- .fixed_level_colours()
  if(length(fixed) && length(m)){
    at <- match(.colour_key(names(m)), .colour_key(names(fixed)))
    m[!is.na(at)] <- unname(fixed[at[!is.na(at)]])
  }
  m
}

#' Level name reduced to a matching key
#'
#' Whitespace is dropped, ordinary and non-breaking alike, so a heading
#' written `COX || CYP` with non-breaking spaces - which is what a
#' spreadsheet usually stores - is configured with ordinary ones.
#' @keywords internal
#' @noRd
.colour_key <- function(x) gsub("[[:space:]\u00a0]+", "", as.character(x))

#' Colours pinned to particular levels by the configuration
#'
#' `getOption("MRManalyzeR.level_colours")`, set by each report from its
#' `colours: levels:` block. A named character vector, or the named list a
#' YAML map arrives as.
#' @return Named character vector, empty when nothing is configured.
#' @keywords internal
#' @noRd
.fixed_level_colours <- function(){
  v <- getOption("MRManalyzeR.level_colours", NULL)
  if(!length(v)) return(character(0))
  v <- unlist(v, use.names = TRUE)
  v <- stats::setNames(as.character(v), names(v))
  v[!is.na(names(v)) & nzchar(names(v)) & !is.na(v)]
}

#' Annotation colour list for pheatmap
#'
#' One named colour vector per annotation column, via
#' [.group_fill_colors()] so a class keeps its colour across the heatmap
#' annotation, the PCA loadings and the correlation strips.
#'
#' @param ann Annotation data frame, or `NULL`.
#' @return A named list of colour vectors, or `NULL`.
#' @keywords internal
#' @noRd
.ann_colours <- function(ann){
  if(is.null(ann)) return(NULL)
  out <- lapply(ann, .group_fill_colors)
  names(out) <- colnames(ann)
  out
}

#' Move columns to the front of a metadata frame
#'
#' The identifier stays first and the named columns follow it, so the two
#' things a reader looks for - what the compound is, and what happened to it -
#' are the first columns of the sheet rather than the last.
#'
#' @param df A data frame.
#' @param cols Columns to pull forward, in the order wanted.
#' @param after Column to keep in front of them.
#' @return `df` with its columns reordered; missing names are ignored.
#' @keywords internal
#' @noRd
.cols_after <- function(df, cols, after = "Compound"){
  cols <- intersect(cols, colnames(df))
  if(!length(cols)) return(df)
  lead <- intersect(after, colnames(df))
  rest <- setdiff(colnames(df), c(lead, cols))
  df[, c(lead, cols, rest), drop = FALSE]
}

#' The lowest value actually measured for each compound
#'
#' Read off the matrix as measured, before any imputation, over the
#' **biological samples only** - so a compound's floor is the smallest number
#' the instrument reported for it in the study material, and a compound with
#' nothing measured there is `NA` rather than `Inf`. QC and blank injections
#' are excluded: a pooled QC is not a study sample, and a blank is the thing
#' the signal is measured against.
#'
#' Rounded to the nearest whole number, because the column is there to be
#' read: a peak area carries no meaning at the eighth decimal place. Note
#' that this is an integer whatever the datatype, so on a run reporting
#' concentrations a floor below 0.5 rounds to zero.
#'
#' @param de A `DatasetExperiment`.
#' @param sample_col,sample_labels `sample_meta` column and values marking the
#'   biological samples.
#' @param exclude Values of `sample_col` that are never samples (the QC and
#'   blank labels), used only when no row matches `sample_labels`.
#' @return Named numeric vector, one entry per compound, rounded to the
#'   nearest whole number.
#' @keywords internal
#' @noRd
.lowest_measured <- function(de, sample_col = "Sample_type",
                             sample_labels = "Sample",
                             exclude = character()){
  x <- as.data.frame(de$data)
  if(!ncol(x)) return(stats::setNames(numeric(0), character(0)))

  sm  <- as.data.frame(de$sample_meta)
  col <- .resolve_meta_col(sample_col, sm)
  rows <- rep(TRUE, nrow(x))
  if(!is.null(col) && nrow(sm) == nrow(x)){
    lv   <- as.character(sm[[col]])
    rows <- lv %in% sample_labels
    # A study that labels its samples something else would otherwise record
    # nothing at all, so fall back to "whatever is not a QC or a blank".
    if(!any(rows)) rows <- !(lv %in% exclude)
  }
  if(!any(rows))
    return(stats::setNames(rep(NA_real_, ncol(x)), colnames(x)))

  v <- suppressWarnings(vapply(x[rows, , drop = FALSE],
                               function(col) min(col, na.rm = TRUE),
                               numeric(1)))
  v[!is.finite(v)] <- NA_real_
  round(v)
}

#' Stop unless this is a DatasetExperiment
#'
#' Every exported step takes the dataset as its first argument, so the common
#' mistake is handing one a plain data frame or a matrix. Checking here names
#' the function and what arrived, instead of failing further in on a slot the
#' object does not have.
#'
#' @param de The object to check.
#' @param fn Name of the calling function, for the message.
#' @return `TRUE`, invisibly; called for the error.
#' @keywords internal
#' @noRd
.check_de <- function(de, fn){
  if(!methods::is(de, "DatasetExperiment"))
    stop(sprintf(
      "[%s] `de` must be a struct::DatasetExperiment, not %s. Build one with assembleDataset(), or read a stored one with loadDataset().",
      fn, paste(class(de), collapse = "/")), call. = FALSE)
  invisible(TRUE)
}

#' Warn before replacing columns a previous run wrote
#'
#' Each step records what it did in `variable_meta` or `sample_meta`. Running
#' the same step twice overwrites that record, and the second set of numbers
#' describes different data from the first - so say so rather than let a
#' recorded value quietly change meaning.
#'
#' @param df The metadata frame about to be written to.
#' @param cols Column names this step is about to set.
#' @param fn Name of the calling function, for the message.
#' @param where `"variable_meta"` or `"sample_meta"`, for the message.
#' @return `invisible(NULL)`; called for the warning.
#' @keywords internal
#' @noRd
.warn_overwrite <- function(df, cols, fn, where = "variable_meta"){
  hit <- intersect(cols, colnames(df))
  if(length(hit))
    warning(sprintf(
      "[%s] overwriting %s in %s: %s. Those values described the previous run.",
      fn, if(length(hit) == 1) "a column" else "columns", where,
      paste(hit, collapse = ", ")), call. = FALSE)
  invisible(NULL)
}
