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
#' analysis; [runExample()] runs the bundled example dataset end to end.
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
#' Can these numbers be trusted? [runPCA()], plus the per-compound CV metrics
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
#' @importFrom utils head
#' @importFrom stats median
#' @importFrom methods new
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
.resolve_meta_col = function(name, df){
  if(is.null(name)) return(NULL)
  nm = as.character(name)[1]
  if(is.na(nm) || !nzchar(nm))  return(NULL)
  if(nm %in% colnames(df))      return(nm)
  alt = make.names(nm)
  if(alt %in% colnames(df))     return(alt)
  hit = colnames(df)[make.names(colnames(df)) == alt]
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
.mrm_theme = function(){
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
.ltc = list(
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
.mrm_palette = function(){
  nm = getOption("MRManalyzeR.palette", "hat")
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
.diverging_pal = function(n = 101){
  nm  = getOption("MRManalyzeR.diverging", "heatmap2")
  pal = .ltc[[nm]] %||% .ltc$heatmap2
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
.qual_pal = function(n){
  if(n < 1) return(character(0))
  pal = .mrm_palette()
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
.group_fill_colors = function(x){
  x  = as.character(x)
  lv = sort(unique(x[!is.na(x) & x != "NA"]))
  m  = if(length(lv) == 0)      character(0)
       else if(length(lv) == 2) stats::setNames(
                                  .mrm_palette()[c(1L, 2L)], lv)
       else                     stats::setNames(.qual_pal(length(lv)), lv)
  if(any(is.na(x) | x == "NA")) m = c(m, "NA" = "grey70")
  m
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
.ann_colours = function(ann){
  if(is.null(ann)) return(NULL)
  out = lapply(ann, .group_fill_colors)
  names(out) = colnames(ann)
  out
}
