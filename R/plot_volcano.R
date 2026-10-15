#' Volcano plot for one comparison
#'
#' Effect size against significance for every compound in a comparison: log2
#' fold-change on x, `-log10(p)` on y. Together they say what neither says
#' alone - a large fold-change on three noisy animals is not a finding, and a
#' tiny but exquisitely reproducible shift usually is not interesting either.
#' The compounds worth attention sit in the upper corners.
#'
#' Rows are taken from the `stats` table returned by [runStats()], so the plot
#' and the `stats` sheet of the output xlsx are the same numbers. Only rows
#' carrying a fold change are eligible, which excludes the global tests
#' (`anova`, `kruskal`) - they ask whether *any* group differs and have no
#' direction. A two-level comparison contributes its single test row; a
#' three-or-more-level comparison contributes one volcano per post-hoc pair,
#' selected with `pair`.
#'
#' Every significant feature is labelled, up to `top_n_label`; beyond that the
#' most significant are kept, since a volcano labelled with eighty compounds
#' communicates nothing.
#'
#' @param stats The `stats` data frame from [runStats()], or the whole list.
#' @param comparison Name of the comparison to plot. `NULL` uses the only one
#'   present, and errors if the table holds more than one.
#' @param pair For post-hoc rows, the group pair as `"A vs B"`. `NULL` uses the
#'   two-group test rows.
#' @param use_adjusted Plot BH-adjusted p on the y-axis instead of raw p.
#' @param sig_threshold Significance cut-off; drawn as a dashed line and used
#'   to colour and label points.
#' @param top_n_label Maximum number of significant features to label.
#' @return A `ggplot`, or `NULL` when the comparison has no plottable rows.
#' @examples
#' de <- loadDataset(system.file("extdata", "example_synthetic.RDS",
#'                                package = "MRManalyzeR"))
#' de <- subsetDataset(de, conditions = list(Sample_type = "Sample"))
#' params <- list(comparisons = list(enabled = TRUE, entries = list(
#'   list(name = "PBS_vs_HDM", method = "welch",
#'        compare = list(factor = "Treatment",
#'                       levels = c("PBS", "HDM"))))))
#' st <- runStats(de, params)
#' plotVolcano(st, "PBS_vs_HDM")
#' @family stats
#' @export
plotVolcano <- function(stats, comparison = NULL, pair = NULL,
                        use_adjusted  = FALSE,
                        sig_threshold = 0.05,
                        top_n_label   = 10){

  if(is.list(stats) && !is.data.frame(stats) && !is.null(stats$stats))
    stats <- stats$stats
  stats <- as.data.frame(stats)

  p_col <- if(isTRUE(use_adjusted)) "p_adj" else "p_value"
  p_lbl <- if(isTRUE(use_adjusted)) "BH-adjusted p" else "raw p"

  posthoc <- c("tukey", "pairwise_wilcoxon")

  df <- stats[!is.na(stats$log2FC) & !is.na(stats[[p_col]]), , drop = FALSE]

  if(is.null(comparison)){
    comps <- unique(df$comparison)
    if(length(comps) > 1)
      stop("[plotVolcano] several comparisons present; name one: ",
           paste(comps, collapse = ", "))
    comparison <- comps[1]
  }
  df <- df[df$comparison %in% comparison, , drop = FALSE]

  if(!is.null(pair)){
    df <- df[df$method %in% posthoc &
              paste(df$group1, "vs", df$group2) %in% pair, , drop = FALSE]
  } else if(any(!df$method %in% posthoc)){
    df <- df[!df$method %in% posthoc, , drop = FALSE]
  }

  if(nrow(df) == 0) return(NULL)

  title <- if(is.null(pair)) comparison else paste0(comparison, " - ", pair)

  df$neglog10p <- -log10(df[[p_col]])
  df$sig      <- ifelse(df[[p_col]] < sig_threshold, "sig", "ns")
  df$tooltip  <- sprintf("%s<br>log2FC: %.3g<br>%s: %.3g",
                         df$feature, df$log2FC, p_lbl, df[[p_col]])

  sig_rows <- df[df$sig == "sig", , drop = FALSE]
  sig_rows <- sig_rows[order(sig_rows[[p_col]]), , drop = FALSE]
  to_label <- utils::head(sig_rows, top_n_label)

  p <- ggplot2::ggplot(df, ggplot2::aes(.data$log2FC, .data$neglog10p,
                                       color = .data$sig,
                                       text = .data$tooltip)) +
    ggplot2::geom_point(alpha = 0.7, size = 2) +
    # Significant points take the palette's red rather than a private hex, so
    # the volcano matches every other figure in the suite.
    ggplot2::scale_color_manual(
      values = c(sig = .mrm_palette()[2], ns = "grey60")) +
    ggplot2::geom_hline(yintercept = -log10(sig_threshold),
                        linetype = "dashed", color = "grey40") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                        color = "grey40") +
    ggplot2::labs(x = "log2 fold-change",
                  y = sprintf("-log10(%s)", p_lbl),
                  color = NULL, title = title) +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::theme(legend.position = "none")

  if(nrow(to_label) > 0){
    # ggrepel keeps labels off the points; fall back to plain text if the
    # (suggested) package is not installed.
    p <- if(requireNamespace("ggrepel", quietly = TRUE))
      p + ggrepel::geom_text_repel(data = to_label,
                                   ggplot2::aes(label = .data$feature),
                                   color = "#d62728", fontface = "bold",
                                   size = 3.2, max.overlaps = Inf, force = 4,
                                   min.segment.length = 0,
                                   segment.color = "#d62728",
                                   segment.size = 0.3, show.legend = FALSE)
    else
      p + ggplot2::geom_text(data = to_label,
                             ggplot2::aes(label = .data$feature),
                             color = "#d62728", fontface = "bold",
                             size = 3.2, vjust = -0.7, show.legend = FALSE)
  }

  p
}
