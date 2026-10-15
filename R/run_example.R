#' Render the bundled example MRManalyzeR reports
#'
#' Convenience wrapper that drives the full workflow against the bundled
#' oxylipin example dataset (a subset of Kolmert et al. 2018,
#' \doi{10.1016/j.prostaglandins.2018.05.005}). It loads
#' `inst/extdata/example_data.xlsx` and `inst/extdata/example_config.yml`,
#' rewrites the YAML's `paths:` block to point at a writable results
#' directory + the bundled xlsx, and calls [`runMRManalyzeR()`]. After the
#' run the generated HTML reports are opened in interactive sessions.
#'
#' @param results_dir Directory to write outputs (xlsx, RDS, HTML reports)
#'   to. Defaults to a fresh `tempdir()` so repeated calls do not collide.
#' @param open Logical. Open the HTML reports in the browser / RStudio
#'   viewer after rendering? Default `TRUE` in interactive sessions.
#' @param render Logical. Render the HTML reports? Default `TRUE`. Set
#'   `FALSE` to build the peak matrix + xlsx/RDS outputs only (no pandoc) -
#'   a quick smoke test of the processing workflow.
#' @return Invisibly, the list returned by [`runMRManalyzeR()`], whose
#'   `output_directory` element is the resolved results directory.
#' @examples
#' # Light run: build the peak matrix + outputs, skip HTML report rendering.
#' res <- mrmRunExample(render = FALSE, open = FALSE)
#' res$output_directory       # where the xlsx + rds landed
#' @family entry points
#' @export
mrmRunExample <- function(results_dir = tempfile("MRManalyzeR_example_"),
                       open        = interactive(),
                       render      = TRUE){

  data_xlsx <- system.file("extdata", "example_data.xlsx",
                          package = "MRManalyzeR")
  base_yaml <- system.file("extdata", "example_config.yml",
                          package = "MRManalyzeR")
  if(!nzchar(data_xlsx) || !file.exists(data_xlsx))
    stop("[mrmRunExample] bundled example_data.xlsx not found ",
         "- reinstall MRManalyzeR.")
  if(!nzchar(base_yaml) || !file.exists(base_yaml))
    stop("[mrmRunExample] bundled example_config.yml not found ",
         "- reinstall MRManalyzeR.")

  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

  # `loadConfig()` already resolves !concat_path so the returned `cfg` has
  # plain string paths. We override the data + result locations to point at
  # the bundled xlsx and a writable temp folder, then write the modified
  # YAML alongside the outputs so the user can inspect / edit it.
  cfg <- loadConfig(base_yaml)
  cfg$project$paths$data_dir   <- dirname(data_xlsx)
  cfg$project$paths$result_dir <- results_dir
  cfg$project$paths$project_dir <- results_dir
  cfg$project$paths$fn         <- tools::file_path_sans_ext(basename(data_xlsx))

  # Light mode: skip the (pandoc-heavy) report renders, keep the matrix build.
  if(!isTRUE(render)){
    cfg$project$data_quality_report$execute <- FALSE
    cfg$project$stats_report$execute        <- FALSE
  }

  out_yaml <- file.path(results_dir, "example_config.yml")
  yaml::write_yaml(cfg, out_yaml)

  message("[mrmRunExample] outputs -> ", results_dir)
  res <- runMRManalyzeR(out_yaml)

  if(isTRUE(open)){
    htmls <- list.files(results_dir, pattern = "\\.html$", full.names = TRUE)
    for(h in htmls) utils::browseURL(h)
  }

  res$output_directory <- results_dir
  invisible(res)
}
