#' Normal Q-Q plot for one compound
#'
#' Plots the sample quantiles of one compound against theoretical normal
#' quantiles, with the reference line, and reports the Shapiro-Wilk statistic
#' in the title. Points falling along the line mean the parametric tests in
#' [runStats()] are on safe ground for that compound; a curve means they are
#' not, and the Wilcoxon or Kruskal-Wallis result reported alongside is the one
#' to trust.
#'
#' Targeted panels are frequently right-skewed, so the same compound is worth
#' looking at under more than one `transform` - which is why the data-quality
#' report renders a tab per transform rather than picking one.
#'
#' Only the biological samples should normally be included: blanks and QCs are
#' different populations and would distort the distribution. Restrict with
#' `levels`.
#'
#' @param de A `struct::DatasetExperiment`.
#' @param compound Feature (column) name to plot.
#' @param transform One of `"none"`, `"log2"` (log2(x + 1)) or `"sqrt"`.
#' @param sample_type_head `sample_meta` column classifying sample type.
#' @param levels Value(s) of `sample_type_head` to include; `NULL` keeps all.
#' @param value_label Axis label for the measured value, e.g. the datatype
#'   units.
#' @return A `ggplot`. When fewer than three non-missing values are available
#'   the plot is an empty panel whose title says so, so a report loop does not
#'   have to special-case sparse compounds.
#' @examples
#' de <- loadDataset(system.file("extdata", "example_synthetic.RDS",
#'                                package = "MRManalyzeR"))
#' plotQQ(de, "Analyte_02", transform = "log2", levels = "Sample")
#' @family QC check
#' @export
plotQQ <- function(de, compound, transform = c("none", "log2", "sqrt"),
                   sample_type_head = "Sample_type",
                   levels           = NULL,
                   value_label      = "value"){
  .check_de(de, "plotQQ")

  transform <- match.arg(transform)
  dm   <- as.data.frame(de$data)
  smeta <- as.data.frame(de$sample_meta)
  if(!compound %in% colnames(dm))
    stop(sprintf("[plotQQ] '%s' is not a feature in this dataset.", compound))

  keep <- if(is.null(levels) || !sample_type_head %in% colnames(smeta))
    rep(TRUE, nrow(dm)) else smeta[[sample_type_head]] %in% levels

  v <- dm[keep, compound]
  v <- v[!is.na(v)]
  v <- switch(transform,
             log2 = suppressWarnings(log2(v + 1)),
             sqrt = suppressWarnings(sqrt(pmax(v, 0))),
             v)

  if(length(v) < 3)
    return(ggplot2::ggplot() + ggplot2::theme_void() +
             ggplot2::labs(title = sprintf(
               "%s - too few non-NA values (n = %d)", compound, length(v))))

  # Shapiro-Wilk is only defined for 3 <= n <= 5000.
  sw_lbl <- ""
  if(length(v) <= 5000){
    sw <- tryCatch(stats::shapiro.test(v), error = function(e) NULL)
    if(!is.null(sw))
      sw_lbl <- sprintf("  Shapiro-Wilk W = %.3f, p = %s", sw$statistic,
                       format.pval(sw$p.value, digits = 3, eps = 1e-4))
  }

  y_lab <- switch(transform,
                 log2 = sprintf("log2(%s + 1) sample quantiles", value_label),
                 sqrt = sprintf("sqrt(%s) sample quantiles", value_label),
                 sprintf("%s sample quantiles", value_label))

  ggplot2::ggplot(data.frame(value = v), ggplot2::aes(sample = .data$value)) +
    ggplot2::stat_qq(size = 2, alpha = 0.7) +
    ggplot2::stat_qq_line(color = "#d62728", linewidth = 0.6) +
    ggplot2::labs(x = "Theoretical normal quantiles", y = y_lab,
                  title = sprintf("%s   (n = %d)%s",
                                  compound, length(v), sw_lbl)) +
    .mrm_theme()
}
