#' Add per-compound CV metrics to a dataset
#'
#' Computes the per-feature coefficient of variation (100 * sd / mean) over the
#' QC injections (`CV_QC`) and over the biological samples (`CV_sample`), plus
#' their ratio (`CV_sample_vs_QC`), and adds all three as columns of
#' `variable_meta`, matched to each feature by row. QC and sample rows are
#' identified from `sample_type_head` in `sample_meta`.
#'
#' `CV_QC` is technical variability: every QC injection is the same material,
#' so its spread is the measurement noise of that compound in that run, and it
#' is the usual quantification filter. `CV_sample` carries the biological
#' spread *plus* that same noise. Their ratio says whether there is biology
#' left to find: comfortably above 1 and the compound varies more between
#' samples than between replicate injections of identical material; near 1 and
#' it does not.
#'
#' The column is named `CV_sample_vs_QC` rather than `CV_sample/CV_QC` because
#' assigning into a `DatasetExperiment` passes the frame through
#' [make.names()], which would rewrite the `/` to a `.` anyway.
#'
#' Calling it twice is harmless - the columns are overwritten, not duplicated.
#' [runMRManalyzeR()] applies it to every run, so a dataset read back from the
#' output RDS already carries these columns.
#'
#' @param de A `struct::DatasetExperiment`.
#' @param sample_type_head `sample_meta` column classifying sample type.
#' @param qc_label Value(s) in `sample_type_head` marking QC injections.
#' @param sample_labels Value(s) marking biological samples.
#' @return `de` with `CV_QC`, `CV_sample` and `CV_sample_vs_QC` added to
#'   `variable_meta`. `CV_QC` is `NA` when the run holds fewer than two QC
#'   injections.
#' @examples
#' de <- loadDataset(system.file("extdata", "example_synthetic.RDS",
#'                                package = "MRManalyzeR"))
#' de <- addCVMetrics(de, qc_label = "QC")
#' head(as.data.frame(de$variable_meta)[, c("Compound", "CV_QC",
#'                                          "CV_sample", "CV_sample_vs_QC")])
#' @family QC check
#' @export
addCVMetrics <- function(de, sample_type_head = "Sample_type",
                          qc_label = "QC", sample_labels = "Sample"){
  .check_de(de, "addCVMetrics")
  dm   <- as.data.frame(de$data)          # samples x features
  smeta <- as.data.frame(de$sample_meta)
  vm   <- as.data.frame(de$variable_meta)
  if(ncol(dm) == 0 || !sample_type_head %in% colnames(smeta)) return(de)

  cv_fn <- function(x) stats::sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE) * 100
  cv_col <- function(idx){
    if(length(idx) < 2){
      v <- rep(NA_real_, ncol(dm)); names(v) <- colnames(dm); return(v)
    }
    round(vapply(dm, function(x) cv_fn(x[idx]), numeric(1)), 1)
  }
  cv_qc    <- cv_col(which(smeta[[sample_type_head]] %in% qc_label))
  cv_sample <- cv_col(which(smeta[[sample_type_head]] %in% sample_labels))

  .warn_overwrite(vm, c("CV_QC", "CV_sample", "CV_sample_vs_QC"),
                  "addCVMetrics")
  key <- rownames(vm)                       # Compound == colnames(dm)
  vm$CV_QC          <- unname(cv_qc[key])
  vm$CV_sample      <- unname(cv_sample[key])
  vm$CV_sample_vs_QC <- round(unname(cv_sample[key]) / unname(cv_qc[key]), 2)

  de$variable_meta <- vm
  de
}
