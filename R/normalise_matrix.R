#' Normalise a dataset by a per-sample metadata column
#'
#' Divides each sample (row) of the dataset by that sample's value in
#' `sample_meta[[column]]` -- for example protein content, tissue weight, or
#' volume of biofluid -- putting all samples on a per-unit basis.
#'
#' @param de A `struct::DatasetExperiment`.
#' @param column `sample_meta` column to divide by.
#' @return `de` with each row of `data` divided by its per-sample divisor.
#' @family peak-matrix processing
#' @examples
#' de <- struct::DatasetExperiment(
#'   data = data.frame(PGE2 = c(10, 20), PGD2 = c(30, 40),
#'                     row.names = c("S1", "S2")),
#'   sample_meta = data.frame(protein_ug = c(2, 4), row.names = c("S1", "S2")),
#'   variable_meta = data.frame(Compound = c("PGE2", "PGD2"),
#'                              row.names = c("PGE2", "PGD2")))
#' normaliseMatrix(de, column = "protein_ug")$data
#' @export
normaliseMatrix <- function(de, column){
  .check_de(de, "normaliseMatrix")

  smeta <- as.data.frame(de$sample_meta)
  if(!column %in% colnames(smeta))
    stop(sprintf("[normaliseMatrix] sample_meta has no '%s' column.", column))

  # data rows and sample_meta rows are aligned in a DatasetExperiment, so the
  # divisor recycles down each feature column.
  de$data <- as.data.frame(de$data) / smeta[[column]]
  de
}
