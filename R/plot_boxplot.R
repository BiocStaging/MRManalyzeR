#' Per-compound box or bar plot by group
#'
#' One compound's values split by a grouping factor - the plot that actually
#' gets looked at when a statistical test comes back significant, because it
#' shows whether the difference rests on the whole group or on one animal.
#'
#' `type = "box"` draws a box with the individual injections jittered over it,
#' which is the honest view for the small groups a targeted study usually has:
#' the reader can count the points. `type = "bar"` draws group means with SD or
#' SE error bars, matching [summariseGroups()] exactly. Reports commonly show
#' both, the box for distribution and the bar for the summary the statistics
#' were computed on.
#'
#' The y-axis is capped at `quantile(value, yaxis_quant) * yaxis_scalar` so one
#' extreme injection cannot flatten the groups.
#'
#' @param de A `struct::DatasetExperiment`.
#' @param compound Feature (column) name to plot.
#' @param group_by `sample_meta` column defining the x-axis groups.
#' @param facet_by Optional `sample_meta` column to facet across, e.g. one
#'   panel per comparison.
#' @param type `"box"` or `"bar"`.
#' @param error_bar `"SD"`, `"SE"` or `"none"` - bars only.
#' @param value_label Y-axis label, e.g. the datatype units.
#' @param sample_id_head `sample_meta` column identifying points in the
#'   tooltip; `NULL` uses row names.
#' @param tooltip_factor Additional `sample_meta` column shown in the tooltip,
#'   typically the primary study factor.
#' @param yaxis_quant,yaxis_scalar Quantile used to cap the y-axis, and the
#'   headroom multiplier above it.
#' @param fill_brewer Use the package's qualitative palette for group fills;
#'   `FALSE` keeps ggplot2's default hue scale.
#' @return A `ggplot`.
#' @examples
#' de <- loadDataset(system.file("extdata", "example_synthetic.RDS",
#'                                package = "MRManalyzeR"))
#' de <- subsetDataset(de, conditions = list(Sample_type = "Sample"))
#' plotBoxplot(de, "Analyte_02", "Treatment", value_label = "ng/mL")
#' plotBoxplot(de, "Analyte_02", "Treatment", type = "bar",
#'              error_bar = "SE")
#' @family stats
#' @seealso [summariseGroups()] for the numbers the bars are drawn from.
#' @export
plotBoxplot <- function(de, compound, group_by, facet_by = NULL,
                        type           = c("box", "bar"),
                        error_bar      = c("SD", "SE", "none"),
                        value_label    = "value",
                        sample_id_head = NULL,
                        tooltip_factor = NULL,
                        yaxis_quant    = 0.95,
                        yaxis_scalar   = 1.3,
                        fill_brewer    = TRUE){
  .check_de(de, "plotBoxplot")

  type <- match.arg(type)

  # YAML supplies these uppercase ("SD" / "SE" / "NONE"); accept either case
  # rather than make every caller normalise.
  error_bar <- toupper(error_bar[1])
  if(!error_bar %in% c("SD", "SE", "NONE"))
    stop("[plotBoxplot] error_bar must be one of: SD, SE, none.")

  dm   <- as.data.frame(de$data)
  smeta <- as.data.frame(de$sample_meta)

  if(!compound %in% colnames(dm))
    stop(sprintf("[plotBoxplot] '%s' is not a feature in this dataset.",
                 compound))
  group_by <- .resolve_meta_col(group_by, smeta) %||%
    stop(sprintf("[plotBoxplot] '%s' is not a sample_meta column. Present: %s.",
                 group_by, paste(colnames(smeta), collapse = ", ")))
  facet_by      <- .resolve_meta_col(facet_by, smeta)
  sample_id_head <- .resolve_meta_col(sample_id_head, smeta)
  tooltip_factor <- .resolve_meta_col(tooltip_factor, smeta)

  df <- data.frame(group = as.character(smeta[[group_by]]),
                  value = dm[[compound]],
                  stringsAsFactors = FALSE)

  sid <- if(!is.null(sample_id_head))
    as.character(smeta[[sample_id_head]]) else rownames(dm)
  extra <- if(!is.null(tooltip_factor))
    as.character(smeta[[tooltip_factor]]) else NA_character_

  df$tooltip <- sprintf("%s: %s<br>%s: %.3g<br>%s: %s",
                       sample_id_head %||% "sample", sid,
                       value_label, df$value,
                       tooltip_factor %||% "group", extra)

  if(!is.null(facet_by)) df$facet <- as.character(smeta[[facet_by]])

  cap <- stats::quantile(df$value, probs = yaxis_quant, na.rm = TRUE) *
        yaxis_scalar

  p <- ggplot2::ggplot(df, ggplot2::aes(.data$group, .data$value,
                                       fill = .data$group)) +
    ggplot2::coord_cartesian(ylim = c(0, cap)) +
    ggplot2::labs(x = NULL, y = value_label, title = compound) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1),
      legend.position = "none")

  if(isTRUE(fill_brewer))
    p <- p + ggplot2::scale_fill_manual(
              values = .group_fill_colors(df$group), na.value = "grey60")

  if(!is.null(facet_by))
    p <- p + ggplot2::facet_wrap(~ facet, scales = "free_x", nrow = 1)

  if(type == "box"){
    # The tooltip is bound only to the jitter layer. Putting `text` on the
    # parent ggplot makes ggplotly() emit one trace per row across every
    # layer, which draws horizontal lines through the boxes.
    p + ggplot2::geom_boxplot(alpha = 0.3, outlier.shape = NA) +
      ggplot2::geom_jitter(ggplot2::aes(text = .data$tooltip),
                           width = 0.15, alpha = 0.7, size = 1.4)
  } else {
    err_fn <- switch(
      error_bar,
      NONE = NULL,
      SD = function(x){
        m <- mean(x, na.rm = TRUE); s <- stats::sd(x, na.rm = TRUE)
        c(y = m, ymin = m - s, ymax = m + s)
      },
      SE = function(x){
        m <- mean(x, na.rm = TRUE); n <- sum(!is.na(x))
        s <- stats::sd(x, na.rm = TRUE) / sqrt(max(1, n))
        c(y = m, ymin = m - s, ymax = m + s)
      })

    # Jitter is deliberately omitted on bars: ggplotly() renders jittered
    # points under a bar as horizontal segments. The box view shows them.
    p <- p + ggplot2::stat_summary(geom = "bar", fun = "mean", alpha = 0.4,
                                  colour = "black", width = 0.7)
    if(!is.null(err_fn))
      # ggplotly does not translate errorbar width 1:1 with bar width; a very
      # small absolute width is what renders as a cap narrower than the bar.
      p <- p + ggplot2::stat_summary(geom = "errorbar", fun.data = err_fn,
                                    width = 0.05, linewidth = 0.4)
    p
  }
}
