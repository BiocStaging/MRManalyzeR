#' Subset a `DatasetExperiment` by sample-metadata conditions
#'
#' Filters `sample_meta` to rows matching all `conditions` (a named list of
#' `column = value` pairs, joined by AND) and slices `data` to the matching
#' samples. Optionally restricts to a feature subset.
#'
#' @param de A `struct::DatasetExperiment`.
#' @param conditions A named list, e.g. `list(Treatment = "HDM", Sex = "Male")`.
#'   Values may be a scalar or vector (interpreted as `%in%`). `NULL` or
#'   empty list returns the input unchanged on the sample dimension.
#' @param features Optional character vector of feature (column) names to
#'   retain. `NULL` keeps all features.
#' @param drop_empty_features If `TRUE` (default), features that are entirely
#'   `NA` after subsetting are dropped from both `data` and `variable_meta`.
#' @return A new `DatasetExperiment` with `data`, `sample_meta`, and
#'   `variable_meta` consistently sliced.
#' @examples
#' m  <- data.frame(A = c(10, 20, 30), B = c(1, 2, 3),
#'                  row.names = c("S1", "S2", "S3"))
#' fm <- data.frame(Compound = c("A", "B"), row.names = c("A", "B"))
#' sm <- data.frame(Sample_ID = c("S1", "S2", "S3"), Group = c("x", "y", "x"),
#'                  row.names = c("S1", "S2", "S3"))
#' de <- struct::DatasetExperiment(data = m, sample_meta = sm, variable_meta = fm)
#' subsetDataset(de, conditions = list(Group = "x"))
#' @family peak-matrix processing
#' @export
subsetDataset <- function(de, conditions = NULL, features = NULL,
                     drop_empty_features = TRUE){
  .check_de(de, "subsetDataset")

  smeta <- de$sample_meta
  data <- de$data
  vmeta <- de$variable_meta

  # --- Sample filter ------------------------------------------------------
  keep_row <- rep(TRUE, nrow(smeta))
  if(length(conditions)){
    for(col in names(conditions)){
      if(!col %in% colnames(smeta))
        stop(sprintf("subsetDataset: column '%s' not in sample_meta.", col))
      keep_row <- keep_row & (smeta[[col]] %in% conditions[[col]])
    }
  }

  smeta_out <- smeta[keep_row, , drop = FALSE]
  if(nrow(smeta_out) == 0)
    warning("subsetDataset: subset produced 0 samples.")

  data_out <- data[keep_row, , drop = FALSE]

  # --- Feature filter -----------------------------------------------------
  if(!is.null(features)){
    miss <- setdiff(features, colnames(data_out))
    if(length(miss))
      warning("subsetDataset: features not in matrix: ",
              paste(head(miss, 5), collapse = ", "),
              if(length(miss) > 5) " ..." else "")
    keep_feat <- intersect(features, colnames(data_out))
    data_out <- data_out[, keep_feat, drop = FALSE]
    vmeta_out <- vmeta[match(keep_feat, vmeta$Compound), , drop = FALSE]
  } else {
    vmeta_out <- vmeta
  }

  # --- Drop all-NA features ----------------------------------------------
  if(isTRUE(drop_empty_features) && ncol(data_out) > 0){
    has_data <- colSums(!is.na(data_out)) > 0
    data_out <- data_out[, has_data, drop = FALSE]
    vmeta_out <- vmeta_out[match(colnames(data_out), vmeta_out$Compound), , drop = FALSE]
  }

  struct::DatasetExperiment(
    data          = data_out,
    sample_meta   = smeta_out,
    variable_meta = vmeta_out
  )
}
