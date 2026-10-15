#' Assemble a struct::DatasetExperiment from a matrix + metadata
#'
#' Aligns `feature_meta` and `sample_meta` to the matrix (features to columns,
#' samples to rows, both by name) and builds the `DatasetExperiment` so
#' `variable_meta`/`sample_meta` rownames match the data assay exactly.
#'
#' @param x Wide numeric matrix / data frame, rows = samples, columns =
#'   compounds (canonical `Compound` names).
#' @param feature_meta Feature metadata; rows with `Report == "YES"` whose
#'   `Compound` is in `x` are kept, in matrix-column order.
#' @param sample_meta Sample metadata; rows whose `name_col` value is a row of
#'   `x` are kept.
#' @param name_col `sample_meta` column matching the matrix row names.
#' @return A `struct::DatasetExperiment`.
#' @examples
#' m  <- data.frame(A = c(10, 20, 30), B = c(1, 2, 3),
#'                  row.names = c("S1", "S2", "S3"))
#' fm <- data.frame(Compound = c("A", "B"), Report = "YES",
#'                  row.names = c("A", "B"))
#' sm <- data.frame(Name = c("S1", "S2", "S3"), Group = c("x", "y", "x"),
#'                  row.names = c("S1", "S2", "S3"))
#' assembleDataset(m, fm, sm)
#' @family data parse
#' @export
assembleDataset <- function(x, feature_meta, sample_meta, name_col = "Name"){

  # Check what was handed over before using it: the failures otherwise surface
  # deep inside dplyr, naming a column the caller never mentioned.
  if(!is.data.frame(x) && !is.matrix(x))
    stop(sprintf(
      "[assembleDataset] `x` must be a matrix or data frame of samples x compounds, not %s.",
      paste(class(x), collapse = "/")), call. = FALSE)
  if(is.null(rownames(x)) || is.null(colnames(x)))
    stop("[assembleDataset] `x` needs row names (injections) and column names (compounds).",
         call. = FALSE)
  for(.a in c("feature_meta", "sample_meta")){
    .v <- get(.a)
    if(!is.data.frame(.v))
      stop(sprintf("[assembleDataset] `%s` must be a data frame, not %s.",
                   .a, paste(class(.v), collapse = "/")), call. = FALSE)
  }
  .need <- function(df, cols, what){
    miss <- setdiff(cols, colnames(df))
    if(length(miss))
      stop(sprintf("[assembleDataset] %s is missing column(s): %s. Present: %s.",
                   what, paste(miss, collapse = ", "),
                   paste(colnames(df), collapse = ", ")), call. = FALSE)
  }
  .need(feature_meta, c("Compound", "Report"), "feature_meta")
  .need(sample_meta,  name_col,               "sample_meta")

  # Feature side: keep only feature_meta rows whose Compound is in the matrix,
  # in matrix-column order, so variable_meta rownames match assay colnames.
  fdata_output <- feature_meta |>
    dplyr::filter(Report == "YES", Compound %in% colnames(x)) |>
    dplyr::arrange(match(Compound, colnames(x)))
  rownames(fdata_output) <- fdata_output$Compound

  # Sample side: intersect first so dropped samples (blank filter, etc.)
  # don't get reintroduced as all-NA rows.
  metadata_ar <- sample_meta[sample_meta[[name_col]] %in% rownames(x), , drop = FALSE]
  metadata_ar <- metadata_ar[order(metadata_ar[[name_col]]), , drop = FALSE]
  rownames(metadata_ar) <- metadata_ar[[name_col]]
  out_matrix_ar <- x[metadata_ar[[name_col]], fdata_output$Compound, drop = FALSE]

  # struct::DatasetExperiment runs make.names() over the metadata it stores, so
  # a heading like `S-group` becomes `S.group`. Said out loud because the
  # symptom otherwise appears far away and looks like nothing: a config asking
  # to colour by `S-group` matches no column, and the plot quietly draws
  # everything in one colour with no error to trace back.
  .report_renamed(fdata_output, "feature_metadata")
  .report_renamed(metadata_ar,  "sample_metadata")

  struct::DatasetExperiment(
    data          = out_matrix_ar,
    sample_meta   = metadata_ar,
    variable_meta = fdata_output
  )
}

#' Announce columns whose names will not survive storage
#' @keywords internal
#' @noRd
.report_renamed <- function(df, what){
  from <- colnames(df)
  to  <- make.names(from, unique = TRUE)
  chg <- from != to
  if(any(chg))
    message(sprintf(
      "[assembleDataset] %s column(s) renamed for storage: %s. Use the new name in the config.",
      what, paste(sprintf("'%s' -> '%s'", from[chg], to[chg]),
                  collapse = ", ")))
}
