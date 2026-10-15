#' Read a plain sample x compound matrix
#'
#' The vendor-neutral route in. `readTargetLynx()` and `readSkyline()` exist
#' because those two exports are awkward shapes; anything that can produce a
#' rectangular table of samples against analytes needs no parser, only a
#' documented schema. That covers Sciex MultiQuant, Agilent MassHunter, Thermo
#' TraceFinder, and the hand-curated spreadsheets that most laboratories end up
#' with at least once.
#'
#' The schema is deliberately minimal: one identifier column, then one column
#' per analyte (or the transpose, with `orientation = "samples_cols"`). The
#' identifiers must match `sample_metadata`'s name column, and the analyte
#' names must match `feature_metadata`; everything else - filtering,
#' normalisation, statistics - is the same from here on.
#'
#' No signal filtering is applied, because a generic matrix carries no S/N or
#' LOD information. Mask values before import, or supply an `LOD`/`LOQ` column
#' in `feature_meta` and use [processDataset()] with `signal_filter`.
#'
#' @param x Path to an xlsx workbook, or a data frame / matrix already in
#'   memory.
#' @param data_tab_names Sheet name(s) to read when `x` is a path. Multiple
#'   sheets are row-bound, so an acquisition split across sheets can be given
#'   as a vector; sample identifiers must be unique across all of them.
#' @param id_col The label column - the one holding identifiers rather than
#'   measurements. It is dropped from the numeric matrix and used as row names.
#'   `NULL` uses the first column. What it identifies follows `orientation`:
#'   sample names under `"samples_rows"`, compound names under
#'   `"samples_cols"`, with the other axis coming from the header row.
#' @param orientation `"samples_rows"` (default: one row per sample) or
#'   `"samples_cols"` (one row per analyte, transposed on read).
#' @return A data frame, rows = samples, columns = compounds, numeric - the
#'   same contract as [readTargetLynx()] and [readSkyline()].
#' @examples
#' m <- data.frame(Name = c("S1", "S2", "S3"),
#'                 PGE2 = c(120, 95, 140),
#'                 PGD2 = c(45, 51, 38))
#' readPeakMatrix(m)
#' @family data parse
#' @export
readPeakMatrix <- function(x, data_tab_names = "matrix_data",
                            id_col = NULL,
                            orientation = c("samples_rows",
                                            "samples_cols")){

  orientation <- match.arg(orientation)

  raw <- if(is.character(x) && length(x) == 1L){
    parts <- lapply(data_tab_names, function(sh)
      openxlsx::read.xlsx(x, sheet = sh))
    do.call(rbind, parts)
  } else {
    as.data.frame(x, check.names = FALSE)
  }

  if(ncol(raw) < 2)
    stop("[readPeakMatrix] need an identifier column plus at least one measurement column.")

  # check.names = FALSE throughout: compound names routinely contain brackets,
  # parentheses and commas that make.names() would mangle, breaking every
  # downstream match against feature_metadata.
  raw <- as.data.frame(raw, check.names = FALSE)
  idc <- if(is.null(id_col)) 1L else match(id_col, colnames(raw))
  if(is.na(idc))
    stop(sprintf("[readPeakMatrix] id_col '%s' is not a column. Present: %s.",
                 id_col, paste(colnames(raw), collapse = ", ")))

  ids <- as.character(raw[[idc]])
  mat <- as.matrix(raw[, -idc, drop = FALSE])
  rownames(mat) <- ids

  if(orientation == "samples_cols"){
    mat <- t(mat)
    ids <- rownames(mat)
  }

  if(anyDuplicated(ids))
    stop(sprintf("[readPeakMatrix] duplicate sample identifier(s): %s",
                 paste(unique(ids[duplicated(ids)]), collapse = ", ")))

  # Sentinel non-detect tokens vary by vendor; resolve the common ones to NA
  # rather than letting them coerce the whole column to character.
  sentinel <- c("#N/A", "N/A", "NA", "n.d.", "ND", "nd", "-", "", "NaN",
               "Inf", "-Inf")
  mat[mat %in% sentinel] <- NA
  suppressWarnings(mode(mat) <- "numeric")

  as.data.frame(mat, stringsAsFactors = FALSE, check.names = FALSE)
}
