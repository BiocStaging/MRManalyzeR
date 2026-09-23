#' Plot one compound's intensity against injection order
#'
#' The core analytical-drift check: a compound's response plotted against the
#' order its injections were acquired in, coloured by sample type so the pooled
#' QCs can be tracked through the sequence. QCs are identical material, so a
#' trend, a step or a widening spread in them is instrumental - loss of
#' sensitivity, a column change, a batch boundary - and not biology.
#'
#' The y-axis is capped at `quantile(value, yaxis_quant) * yaxis_scalar` so a
#' single large value cannot flatten the rest of the series; points above the
#' cap are drawn outside the panel rather than dropped.
#'
#' The returned plot carries a `text` aesthetic holding the per-point tooltip,
#' which [plotly::ggplotly()] picks up with `tooltip = "text"`; it is ignored
#' by a plain `print()`.
#'
#' @param de A `struct::DatasetExperiment`.
#' @param compound Feature (column) name to plot.
#' @param sample_type_head `sample_meta` column classifying sample type.
#' @param levels Value(s) of `sample_type_head` to include; `NULL` keeps all.
#' @param injection_order_head `sample_meta` column holding acquisition order.
#'   Must be numeric; falls back to row position when the column is absent.
#' @param sample_id_head `sample_meta` column used to label points in the
#'   tooltip; `NULL` uses the row names.
#' @param value_label Axis / tooltip label for the measured value, e.g. the
#'   datatype units ("ng/mL", "Area").
#' @param yaxis_quant,yaxis_scalar Quantile used to cap the y-axis, and the
#'   headroom multiplier applied above it.
#' @return A `ggplot`.
#' @examples
#' de <- loadDataset(system.file("extdata", "example_synthetic.RDS",
#'                                package = "MRManalyzeR"))
#' plotDrift(de, "Analyte_02", value_label = "ng/mL")
#' @family QC check
#' @export
plotDrift = function(de, compound,
                      sample_type_head     = "Sample_type",
                      levels               = NULL,
                      injection_order_head = "Injection_order",
                      sample_id_head       = NULL,
                      value_label          = "value",
                      yaxis_quant          = 0.95,
                      yaxis_scalar         = 1.3){

  dm    = as.data.frame(de$data)
  smeta = as.data.frame(de$sample_meta)
  if(!compound %in% colnames(dm))
    stop(sprintf("[plotDrift] '%s' is not a feature in this dataset.",
                 compound))

  keep = if(is.null(levels) || !sample_type_head %in% colnames(smeta))
    rep(TRUE, nrow(dm)) else smeta[[sample_type_head]] %in% levels

  # A text injection order (a formula returning "_48", say) would become NA
  # and every point would silently drop out of the plot.
  ord = which(keep)
  if(injection_order_head %in% colnames(smeta)){
    raw = smeta[[injection_order_head]][keep]
    ord = suppressWarnings(as.numeric(as.character(raw)))
    bad = is.na(ord) & !is.na(raw)
    if(any(bad))
      stop(sprintf("[plotDrift] '%s' must be numeric; found: %s.",
                   injection_order_head,
                   paste(sprintf("'%s'", utils::head(
                           unique(as.character(raw[bad])), 5)),
                         collapse = ", ")), call. = FALSE)
  }

  df = data.frame(
    value = dm[keep, compound],
    order = ord,
    type  = if(sample_type_head %in% colnames(smeta))
      as.character(smeta[[sample_type_head]][keep]) else "sample",
    id    = if(!is.null(sample_id_head) && sample_id_head %in% colnames(smeta))
      as.character(smeta[[sample_id_head]][keep]) else rownames(dm)[keep],
    stringsAsFactors = FALSE)

  df$tooltip = sprintf("%s: %s<br>%s: %.3g<br>%s: %s",
                       sample_id_head %||% "sample", df$id,
                       value_label, df$value, sample_type_head, df$type)

  cap = stats::quantile(df$value, probs = yaxis_quant, na.rm = TRUE) *
        yaxis_scalar

  ggplot2::ggplot(df, ggplot2::aes(.data$order, .data$value)) +
    ggplot2::geom_point(ggplot2::aes(color = .data$type, text = .data$tooltip),
                        size = 3, alpha = 0.6, na.rm = TRUE) +
    ggplot2::scale_color_manual(values = .group_fill_colors(df$type),
                                na.value = "grey60") +
    ggplot2::coord_cartesian(ylim = c(0, cap)) +
    ggplot2::labs(x = "Injection order", y = value_label, title = compound,
                  color = NULL) +
    .mrm_theme()
}
