#' Summarise which compounds were reported and which were dropped
#'
#' Returns one row per compound the run started with, marking whether it
#' survived into the reported matrix and why not if it did not. Exclusions
#' accumulate from several places - signal-to-noise or LOD/LOQ masking that
#' left a feature entirely missing, features whose `Processing_name` had no
#' rows in the source data, and source columns with no `feature_metadata` entry
#' - so the reason recorded in `Comment` is often the only trace of what
#' happened to a compound.
#'
#' This is the table the data-quality report shows as "compounds included" and
#' "compounds excluded"; keeping it as one frame with a `Status` column makes
#' it sortable and searchable in one place, and writable to a single sheet.
#'
#' @param de A `struct::DatasetExperiment` - the reported dataset.
#' @param removed_features Optional data frame of dropped features, as returned
#'   in `[[2]]` by [processDataset()]. `NULL` (or zero rows) yields an
#'   included-only table, which is what happens when a run reuses a stored RDS
#'   rather than reprocessing.
#' @return A data frame with columns `Compound`, `Status` (`"included"` /
#'   `"excluded"`) and `Comment`.
#' @examples
#' de <- loadDataset(system.file("extdata", "example_synthetic.RDS",
#'                                package = "MRManalyzeR"))
#' head(summariseFeatures(de))
#' @family QC check
#' @export
summariseFeatures <- function(de, removed_features = NULL){
  .check_de(de, "summariseFeatures")

  vm <- as.data.frame(de$variable_meta)

  pull <- function(df, status){
    if(is.null(df) || nrow(df) == 0)
      return(data.frame(Compound = character(0), Status = character(0),
                        Comment = character(0), stringsAsFactors = FALSE))
    data.frame(
      Compound = as.character(df$Compound),
      Status   = status,
      Comment  = if("Comment" %in% colnames(df))
        as.character(df$Comment) else NA_character_,
      stringsAsFactors = FALSE)
  }

  out <- rbind(pull(vm, "included"), pull(removed_features, "excluded"))
  rownames(out) <- NULL
  out
}
