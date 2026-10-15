#' Plot PCA scores, coloured by a sample-metadata column
#'
#' Draws the scores from [mrmRunPCA()] with the percentage of variance explained
#' on each axis. The colouring is the whole point: the same scores answer a
#' different question depending on what you colour them by. Colour by sample
#' type and the question is whether the QC injections cluster tightly in the
#' middle of the samples. Colour by chromatographic batch and it is whether the
#' variance driving component 1 is analytical rather than biological - if the
#' batches separate, correct before interpreting anything. Colour by the study
#' factor and it is whether there is any separation to find at all.
#'
#' Which is why the reports draw one plot per column listed under
#' `qc_pca.color_by` rather than choosing one.
#'
#' Numeric colourings (injection order, a numeric batch id) get a continuous
#' viridis scale; categorical ones get a discrete palette. Pass `numeric`
#' explicitly to override the automatic choice - useful when a column is stored
#' as character but is really a sequence.
#'
#' @param pca The list returned by [mrmRunPCA()].
#' @param sample_meta Sample metadata, or a `DatasetExperiment` to take it
#'   from. Rows are matched to the PCA scores by name, so samples dropped by
#'   `sample_na_max` do not have to be removed by the caller.
#' @param colour_by `sample_meta` column to colour points by.
#' @param components Length-2 integer vector: which principal components to
#'   plot.
#' @param label `"none"`, `"all"`, or `"outliers"` (beyond 2.5 SD on either
#'   component).
#' @param label_factor `sample_meta` column supplying point labels and the
#'   tooltip id; `NULL` uses row names.
#' @param ellipse `"none"` or `"group"` - a normal-probability ellipse per
#'   level of `colour_by`, drawn only for categorical colourings with more than
#'   one level.
#' @param ellipse_level Confidence level for the ellipse.
#' @param numeric `NA` (default) to decide from the column's type, or
#'   `TRUE`/`FALSE` to force a continuous or discrete colour scale.
#' @return A `ggplot`, or `NULL` if the PCA could not be fit (`pca$pr` is
#'   `NULL`), which lets a report skip the section without erroring.
#' @examples
#' de  <- loadDataset(system.file("extdata", "example_synthetic.RDS",
#'                                 package = "MRManalyzeR"))
#' pca <- mrmRunPCA(de, transform = "log2")
#' mrmPlotPCA(pca, de, colour_by = "Chrom_Batch")
#' @family QC check
#' @export
mrmPlotPCA <- function(pca, sample_meta, colour_by,
                    components    = c(1, 2),
                    label         = c("none", "all", "outliers"),
                    label_factor  = NULL,
                    ellipse       = c("none", "group"),
                    ellipse_level = 0.95,
                    numeric       = NA){

  label  <- match.arg(label)
  ellipse <- match.arg(ellipse)

  if(is.null(pca$pr)) return(NULL)

  if(methods::is(sample_meta, "DatasetExperiment"))
    sample_meta <- sample_meta$sample_meta
  smeta <- as.data.frame(sample_meta)

  x <- pca$pr$x
  pc <- as.integer(components)[c(1L, 2L)]
  if(any(pc > ncol(x)))
    stop(sprintf("[mrmPlotPCA] requested component %d but the fit has %d.",
                 max(pc), ncol(x)))

  # Re-align: mrmRunPCA() may drop samples via sample_na_max.
  smeta <- smeta[rownames(x), , drop = FALSE]
  colour_by <- .resolve_meta_col(colour_by, smeta) %||%
    stop(sprintf("[mrmPlotPCA] '%s' is not a sample_meta column. Present: %s.",
                 colour_by, paste(colnames(smeta), collapse = ", ")))
  label_factor <- .resolve_meta_col(label_factor, smeta)

  vexp <- (pca$pr$sdev^2) / sum(pca$pr$sdev^2) * 100

  df <- data.frame(PCx = x[, pc[1]], PCy = x[, pc[2]],
                  colour = smeta[[colour_by]], stringsAsFactors = FALSE)
  df$label <- if(!is.null(label_factor))
    as.character(smeta[[label_factor]]) else rownames(x)

  is_num <- if(is.na(numeric)) is.numeric(df$colour) else isTRUE(numeric)
  if(is_num) df$colour <- suppressWarnings(as.numeric(df$colour))

  df$tooltip <- sprintf("%s<br>%s: %s<br>PC%d: %.3g<br>PC%d: %.3g",
                       df$label, colour_by, as.character(df$colour),
                       pc[1], df$PCx, pc[2], df$PCy)

  p <- ggplot2::ggplot(df, ggplot2::aes(.data$PCx, .data$PCy,
                                       text = .data$tooltip)) +
    ggplot2::geom_point(ggplot2::aes(color = .data$colour),
                        size = 3, alpha = 0.85) +
    ggplot2::labs(x = sprintf("PC%d (%.1f%%)", pc[1], vexp[pc[1]]),
                  y = sprintf("PC%d (%.1f%%)", pc[2], vexp[pc[2]]),
                  color = colour_by) +
    # Scores panels are square and legend-heavy, so they use a plainer,
    # larger-based theme than the per-compound plots.
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::theme(
      legend.position   = "bottom",
      legend.box.margin = ggplot2::margin(t = 12),
      axis.title.x      = ggplot2::element_text(
                            margin = ggplot2::margin(t = 6, b = 0)),
      plot.margin       = ggplot2::margin(4, 6, 4, 6))

  p <- if(is_num) p + ggplot2::scale_color_viridis_c()
      else       p + ggplot2::scale_color_manual(
                       values = .group_fill_colors(df$colour),
                       na.value = "grey60")

  if(ellipse == "group" && !is_num && length(unique(df$colour)) > 1)
    p <- p + ggplot2::stat_ellipse(ggplot2::aes(group = .data$colour),
                                  type = "norm", level = ellipse_level,
                                  linewidth = 0.4)

  if(label == "all"){
    p <- p + ggplot2::geom_text(ggplot2::aes(label = .data$label),
                               size = 2.5, vjust = -1, show.legend = FALSE)
  } else if(label == "outliers"){
    out <- abs(scale(df$PCx)[, 1]) > 2.5 | abs(scale(df$PCy)[, 1]) > 2.5
    if(any(out))
      p <- p + ggplot2::geom_text(data = df[out, , drop = FALSE],
                                 ggplot2::aes(label = .data$label),
                                 size = 2.5, vjust = -1, show.legend = FALSE)
  }

  p
}
