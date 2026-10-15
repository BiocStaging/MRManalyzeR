#' Read Skyline sheet(s) into a wide sample x compound matrix
#'
#' Reads the Skyline "molecule x sample" sheet(s), renaming molecule columns to
#' the canonical `Processing_name` where matched in `feature_meta`, and
#' optionally masks values below the per-compound `LOD`/`LOQ` threshold.
#'
#' @param xlsx_path Path to the Skyline xlsx workbook.
#' @param data_tab_names Sheet name(s) holding the molecule x sample matrix.
#' @param feature_meta Feature metadata (canonical `Processing_name` / `Report`
#'   columns; plus the `LOD`/`LOQ` column when `signal_filter` is set).
#' @param signal_filter `"LOD"`, `"LOQ"`, or `FALSE` - the per-compound floor
#'   below which values are set to `NA`.
#' @return A data frame, rows = samples, columns = compounds, numeric.
#' @examples
#' sky <- data.frame(Molecule = c("PGE2", "PGD2"),
#'                   S1 = c(100, 50), S2 = c(120, 55))
#' f <- tempfile(fileext = ".xlsx")
#' openxlsx::write.xlsx(list(skyline_data = sky), f)
#' fm <- data.frame(Processing_name = c("PGE2", "PGD2"), Report = "YES",
#'                  row.names = c("PGE2", "PGD2"))
#' readSkyline(f, "skyline_data", fm)
#' @family data parse
#' @export
readSkyline <- function(xlsx_path, data_tab_names = "skyline_data",
                        feature_meta, signal_filter = FALSE){

  out_table <- .read_skyline_matrix(xlsx_path, data_tab_names, feature_meta)

  if(!isFALSE(signal_filter) && signal_filter %in% c("LOD", "LOQ"))
    out_table <- .apply_lod_loq_mask(out_table, feature_meta, floor_col = signal_filter)

  out_table
}

#' Read one or more Skyline "molecule x sample" sheets into a
#' sample x compound matrix
#'
#' Skyline exports are already a matrix (rows = molecules, columns =
#' samples) rather than the stacked long-format TargetLynx layout, so this
#' reader does not need [`extractTable()`]'s block-parsing logic - it reads,
#' transposes, coerces sentinel non-detect tokens to `NA`, and (for multiple
#' sheets, e.g. one per acquisition batch) validates and row-binds into a
#' single matrix.
#'
#' The returned matrix uses the *same contract* as
#' `.build_wide_matrix()` + `.apply_snr_mask()` produce for the TargetLynx
#' path: rows named by sample, columns named by `Processing_name` wherever a
#' match was found (case-insensitively) in `fdata`, left as the raw Skyline
#' `Molecule` name otherwise. That means the existing
#' `generateDataMatrix()` mismatch-diagnostics / rename-to-Compound /
#' drop-all-NA-columns logic downstream applies unchanged to both readers.
#'
#' @param xlsx_path Path to the Skyline workbook.
#' @param data_tab_names Character vector of sheet name(s) holding the
#'   molecule x sample matrix. Multiple names are treated as separate
#'   acquisition batches and row-bound together; each must report the same
#'   set of `Report == "YES"` compounds (case-insensitively), and sample
#'   names must be unique across all listed sheets. Sheet-membership does
#'   *not* imply batch-correction grouping - that remains driven entirely by
#'   `sample_metadata` (same as multiple `lcms_data*` sheets on the
#'   TargetLynx path).
#' @param fdata Feature metadata; used only for the case-insensitive
#'   `Processing_name` lookup and the multi-sheet coverage check.
#' @return A data frame: rows = samples, columns = compounds (named by
#'   canonical `Processing_name` where matched, raw `Molecule` name
#'   otherwise), numeric values with sentinel tokens already resolved to `NA`.
#' @keywords internal
#' @noRd
.read_skyline_matrix <- function(xlsx_path, data_tab_names, fdata){

  fnames_yes <- fdata$Processing_name[fdata$Report == "YES"]
  # lower(Processing_name) -> canonical-cased Processing_name, for a
  # case-insensitive join against Skyline's Molecule names.
  ci_lookup <- stats::setNames(fdata$Processing_name, tolower(fdata$Processing_name))

  sheet_dfs <- vector("list", length(data_tab_names))
  covered_yes <- vector("list", length(data_tab_names))

  for(i in seq_along(data_tab_names)){
    sh <- data_tab_names[i]
    raw <- openxlsx::read.xlsx(xlsx_path, sheet = sh)

    if(ncol(raw) < 2)
      stop(sprintf("[readSkyline] sheet '%s' has no sample columns.", sh))

    molecules <- as.character(raw[[1]])
    if(anyDuplicated(molecules))
      stop(sprintf(
        "[readSkyline] sheet '%s' has duplicate Molecule entries: %s",
        sh, paste(unique(molecules[duplicated(molecules)]), collapse = ", ")))

    mat <- as.matrix(raw[, -1, drop = FALSE])
    rownames(mat) <- molecules

    # Sentinel non-numeric tokens -> NA before numeric coercion. Skyline
    # writes literal "#N/A" for non-detects; "221e" (infinity) shows up
    # for the occasional divide-by-zero ratio. The rest are defensive.
    sentinel <- c("#N/A", "221e", "-221e", "NaN", "Inf", "-Inf", "")
    mat[mat %in% sentinel] <- NA
    suppressWarnings(mode(mat) <- "numeric")

    mat_t <- t(mat)   # -> samples x molecules

    # check.names = FALSE: molecule names routinely contain characters
    # (brackets, parens, commas - e.g. "[D4]-PGE2", "11(12)-EpETrE") that
    # as.data.frame()'s default make.names() would otherwise mangle,
    # breaking every downstream name match.
    df <- as.data.frame(mat_t, stringsAsFactors = FALSE, check.names = FALSE)
    df$.sample_name <- rownames(mat_t)
    sheet_dfs[[i]] <- df

    covered_yes[[i]] <- sort(intersect(tolower(molecules), tolower(fnames_yes)))
  }
  names(sheet_dfs) <- data_tab_names

  # --- Validate: all sheets cover the same Report=="YES" compound set -----
  if(length(data_tab_names) > 1){
    ref <- covered_yes[[1]]
    bad <- vapply(covered_yes[-1], function(x) !identical(x, ref), logical(1))
    if(any(bad)){
      bad_sheet <- data_tab_names[-1][which(bad)[1]]
      stop(sprintf(
        "[readSkyline] data_tab_names sheets disagree on which Report=='YES' compounds they cover ('%s' vs '%s'). All listed sheets must report the same reportable-compound set.",
        data_tab_names[1], bad_sheet))
    }
  }

  # --- Validate: sample names unique across sheets ------------------------
  all_samples <- unlist(lapply(sheet_dfs, function(d) d$.sample_name))
  if(anyDuplicated(all_samples))
    stop(sprintf(
      "[readSkyline] duplicate sample name(s) across data_tab_names sheets: %s",
      paste(unique(all_samples[duplicated(all_samples)]), collapse = ", ")))

  # --- Merge: bind_rows unions columns (differing non-reportable/IS
  #     columns per batch are filled NA rather than erroring) --------------
  merged <- dplyr::bind_rows(sheet_dfs)
  merged <- as.data.frame(merged, check.names = FALSE)
  rownames(merged) <- merged$.sample_name
  merged$.sample_name <- NULL

  # --- Case-insensitive rename of matched columns to canonical Processing_name
  orig_cols <- colnames(merged)
  matched  <- tolower(orig_cols) %in% names(ci_lookup)
  new_cols <- orig_cols
  new_cols[matched] <- ci_lookup[tolower(orig_cols[matched])]
  colnames(merged) <- new_cols

  merged
}

#' Mask a wide matrix using a per-compound LOD or LOQ threshold
#'
#' Distinct from [`.apply_snr_mask()`]: SNR masks per-injection (the same
#' compound can be masked in one sample but not another, because S/N is
#' measured per peak). LOD/LOQ masks per-compound (one threshold per
#' compound from `feature_metadata`, applied against every sample's value
#' for that compound) - the threshold is constant, but which cells end up
#' `NA` still varies per sample because the *measured value* does.
#'
#' @param out_table Wide matrix, samples x compounds (canonical
#'   `Processing_name` columns).
#' @param fdata Feature metadata containing the threshold column.
#' @param floor_col Column name in `fdata` to use as the per-compound
#'   threshold - `"LOD"` or `"LOQ"`.
#' @return `out_table` with sub-threshold cells set to `NA`.
#' @keywords internal
#' @noRd
.apply_lod_loq_mask <- function(out_table, fdata, floor_col){

  thresh_lookup <- stats::setNames(
    suppressWarnings(as.numeric(fdata[[floor_col]])),
    fdata$Processing_name
  )

  shared <- intersect(colnames(out_table), names(thresh_lookup))
  thresh <- thresh_lookup[shared]
  has_thresh <- shared[!is.na(thresh)]

  for(cn in has_thresh){
    below <- !is.na(out_table[[cn]]) & out_table[[cn]] < thresh_lookup[[cn]]
    out_table[[cn]][below] <- NA
  }

  out_table
}
