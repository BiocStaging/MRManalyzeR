#' Per-group summary statistics for every compound
#'
#' The numbers behind a boxplot: n, mean, standard deviation and standard error
#' for each compound in each group. [plotBoxplot()] draws bars and error bars
#' from exactly these values, and [runMRManalyzeR()] writes them to the
#' `summary` sheet of the stats workbook - so what a reader measures off the
#' figure and what they read in the spreadsheet are the same numbers.
#'
#' Standard deviation describes the spread of the animals; standard error
#' describes how well the mean is pinned down, and shrinks with n. They answer
#' different questions and the choice belongs to the reader, so both are
#' returned and `plotBoxplot(error_bar =)` picks one.
#'
#' @param de A `struct::DatasetExperiment`.
#' @param group_by `sample_meta` column defining the groups.
#' @param compounds Features to summarise; `NULL` uses all of them.
#' @return A data frame with one row per compound per group: `Compound`,
#'   `group`, `n` (non-missing values), `mean`, `sd`, `se`.
#' @examples
#' de <- loadDataset(system.file("extdata", "example_synthetic.RDS",
#'                                package = "MRManalyzeR"))
#' de <- subsetDataset(de, conditions = list(Sample_type = "Sample"))
#' head(summariseGroups(de, "Treatment"))
#' @family stats
#' @export
summariseGroups <- function(de, group_by, compounds = NULL){
  .check_de(de, "summariseGroups")

  dm   <- as.data.frame(de$data)
  smeta <- as.data.frame(de$sample_meta)

  if(!group_by %in% colnames(smeta))
    stop(sprintf("[summariseGroups] '%s' is not a sample_meta column.",
                 group_by))

  if(is.null(compounds)) compounds <- colnames(dm)
  miss <- setdiff(compounds, colnames(dm))
  if(length(miss))
    stop(sprintf("[summariseGroups] not features in this dataset: %s",
                 paste(utils::head(miss, 5), collapse = ", ")))

  grp <- as.character(smeta[[group_by]])

  rows <- lapply(compounds, function(f){
    v <- dm[[f]]
    sp <- split(v, grp)
    data.frame(
      Compound = f,
      group    = names(sp),
      n        = vapply(sp, function(x) sum(!is.na(x)), integer(1)),
      mean     = vapply(sp, function(x) mean(x, na.rm = TRUE), numeric(1)),
      sd       = vapply(sp, function(x) stats::sd(x, na.rm = TRUE), numeric(1)),
      stringsAsFactors = FALSE)
  })

  out <- do.call(rbind, rows)
  out$se <- out$sd / sqrt(pmax(out$n, 1))
  rownames(out) <- NULL
  out
}
