#' Read a TargetLynx workbook into a wide sample x compound matrix
#'
#' Reads the TargetLynx sheet(s) with [extractTable()], pivots the long table
#' to a wide samples x compound matrix of the requested `datatype`, and
#' optionally masks values whose per-injection S/N is below `snr`. Columns are
#' the raw TargetLynx compound (`Processing_name`) identifiers; the
#' rename-to-`Compound` and feature/sample filtering happen later in
#' [processDataset()].
#'
#' @param xlsx_path Path to the TargetLynx xlsx workbook.
#' @param datatype TargetLynx column to report (e.g. "Area", "Response",
#'   "ng/mL").
#' @param tl_headers TargetLynx column headers to extract.
#' @param snr S/N threshold; values below it are set to `NA`. `FALSE` to skip.
#' @param data_tab_names Sheet name(s) to read; `NULL` auto-matches sheets
#'   whose name contains "lcms_data".
#' @return A data frame, rows = samples, columns = compounds
#'   (`Processing_name`), numeric.
#' @examples
#' xlsx <- system.file("extdata", "example_data.xlsx", package = "MRManalyzeR")
#' m <- readTargetLynx(xlsx, datatype = "Area", snr = 3)
#' dim(m)
#' @family data parse
#' @export
readTargetLynx = function(xlsx_path, datatype = "Area",
                           tl_headers = c("ID", "Name", "Area", "ng/mL", "Response", "S/N"),
                           snr = FALSE, data_tab_names = NULL){

  .read_targetlynx(xlsx_path, datatype = datatype, tl_headers = tl_headers,
                   snr = snr, data_tab_names = data_tab_names)$masked
}

#' Read a TargetLynx workbook, keeping the matrix before and after S/N masking
#'
#' [processDataset()] needs both: the masked matrix to process, and the
#' unmasked one to tell a compound that was never measured from one the S/N
#' filter removed.
#' @return A list with `raw` (before masking) and `masked` (after; the same
#'   object when `snr` is `FALSE`).
#' @keywords internal
#' @noRd
.read_targetlynx = function(xlsx_path, datatype, tl_headers, snr,
                            data_tab_names){
  lcms_table = extractTable(xlsx_path, tl_headers = tl_headers,
                            data_tab_names = data_tab_names) %>%
    subset(Name != "")

  raw = .build_wide_matrix(lcms_table, datatype)
  list(raw    = raw,
       masked = if(isFALSE(snr)) raw
                else .apply_snr_mask(raw, lcms_table, snr))
}
