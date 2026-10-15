#' Load a `DatasetExperiment` from `.RDS` or `.xlsx`
#'
#' Used by [`combineDatasets()`] so the combine workflow can mix RDS and
#' xlsx inputs (xlsx loading reconstructs a `DatasetExperiment` from the
#' standard sheets `feature_metadata`, `sample_metadata`, `matrix` written
#' by [`runMRManalyzeR()`]).
#'
#' @param path Path to either an `.RDS` containing a
#'   `struct::DatasetExperiment`, or an `.xlsx` produced by
#'   [`runMRManalyzeR()`] (with `feature_metadata`, `sample_metadata`,
#'   `matrix` tabs).
#' @param sample_id_col Column in `sample_metadata` used as the row-name
#'   key when loading from xlsx. Default `"Sample_ID"`.
#' @return A `struct::DatasetExperiment`.
#' @examples
#' de <- struct::DatasetExperiment(
#'   data          = data.frame(A = c(1, 2), B = c(3, 4),
#'                              row.names = c("S1", "S2")),
#'   sample_meta   = data.frame(Sample_ID = c("S1", "S2"),
#'                              row.names = c("S1", "S2")),
#'   variable_meta = data.frame(Compound = c("A", "B"),
#'                              row.names = c("A", "B")))
#' f <- tempfile(fileext = ".RDS"); saveRDS(de, f)
#' loadDataset(f)
#' @family data parse
#' @export
loadDataset <- function(path, sample_id_col = "Sample_ID"){

  if(!file.exists(path)) stop("[loadDataset] file not found: ", path)

  ext <- tolower(tools::file_ext(path))

  if(ext == "rds"){
    de <- readRDS(path)
    # Validate variable_meta rownames against data colnames. Old RDS files
    # can have drifted because Compound values may not match the matrix
    # column headers exactly. If they don't match, rebuild the DE so the
    # variable_meta is in lockstep with the data.
    dat <- as.data.frame(de$data)
    vm <- as.data.frame(de$variable_meta)
    sm <- as.data.frame(de$sample_meta)
    dat_cn <- colnames(dat)
    if(!is.null(dat_cn) && !identical(as.character(rownames(vm)),
                                      as.character(dat_cn))){
      if(!"Compound" %in% colnames(vm))
        stop(sprintf("[loadDataset] %s: variable_meta has no Compound column to align on.", path))
      missing_in_vm <- setdiff(dat_cn, as.character(vm$Compound))
      if(length(missing_in_vm)){
        stop(sprintf(
          "[loadDataset] %s: %d data column(s) have no matching Compound row in variable_meta (e.g. %s). The RDS appears to have drifted between data and variable_meta - re-export from the source xlsx.",
          path, length(missing_in_vm),
          paste(head(missing_in_vm, 5), collapse = ", ")))
      }
      rownames(vm) <- as.character(vm$Compound)
      vm          <- vm[dat_cn, , drop = FALSE]
      de <- struct::DatasetExperiment(
        name          = de$name,
        data          = dat,
        sample_meta   = sm,
        variable_meta = vm
      )
    }
    return(de)
  }

  if(ext == "xlsx"){
    sheets <- openxlsx::getSheetNames(path)
    need  <- c("feature_metadata", "sample_metadata", "matrix")
    miss  <- setdiff(need, sheets)
    if(length(miss))
      stop(sprintf("[loadDataset] %s is missing required sheet(s): %s",
                   path, paste(miss, collapse = ", ")))

    fmeta <- openxlsx::read.xlsx(path, sheet = "feature_metadata")
    smeta <- openxlsx::read.xlsx(path, sheet = "sample_metadata")
    M    <- openxlsx::read.xlsx(path, sheet = "matrix", rowNames = TRUE)
    M    <- as.matrix(M)

    if(!sample_id_col %in% colnames(smeta))
      stop(sprintf("[loadDataset] %s: sample_metadata lacks column '%s'.",
                   path, sample_id_col))

    # Align matrix rows to sample_meta in declared order, keyed by sample_id_col.
    sid <- as.character(smeta[[sample_id_col]])
    if(all(sid %in% rownames(M))){
      M <- M[sid, , drop = FALSE]
    } else if(nrow(M) == nrow(smeta)){
      # Fall back to row-order alignment.
      rownames(M) <- sid
    } else {
      stop(sprintf(
        "[loadDataset] %s: cannot align matrix rows to sample_metadata (%d vs %d).",
        path, nrow(M), nrow(smeta)))
    }
    rownames(smeta) <- sid

    # Variable_meta keyed by Compound (the convention used elsewhere in the package).
    if("Compound" %in% colnames(fmeta))
      rownames(fmeta) <- as.character(fmeta$Compound)

    # Ensure variable_meta rows are in the same order as matrix columns.
    # Reorder by name when possible; otherwise error out loudly so the user
    # can fix the input rather than getting silent misalignment downstream.
    if(all(colnames(M) %in% rownames(fmeta))){
      fmeta <- fmeta[colnames(M), , drop = FALSE]
    } else {
      missing_in_fmeta <- setdiff(colnames(M), rownames(fmeta))
      stop(sprintf(
        "[loadDataset] %s: %d matrix column(s) have no matching row in feature_metadata$Compound (e.g. %s). Fix the xlsx so matrix headers and feature_metadata$Compound use identical names.",
        path, length(missing_in_fmeta),
        paste(head(missing_in_fmeta, 5), collapse = ", ")))
    }

    return(struct::DatasetExperiment(
      name          = tools::file_path_sans_ext(basename(path)),
      data          = as.data.frame(M),
      sample_meta   = smeta,
      variable_meta = fmeta
    ))
  }

  stop("[loadDataset] unsupported extension '.", ext, "' for: ", path)
}
