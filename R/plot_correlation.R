#' Feature-by-feature correlation heatmap
#'
#' Draws one triangle of the correlation matrix for a set of compounds. Where
#' the comparisons ask whether a compound differs between groups, this asks
#' which compounds move *together* - which in a targeted panel usually means
#' they share a branch of the pathway, so blocks of correlation along the
#' diagonal are the pathway structure showing itself.
#'
#' Only one triangle is drawn, because a correlation matrix is symmetric and
#' showing both halves doubles the ink for no information. The diagonal is left
#' grey unless self-correlations are present in the input.
#'
#' Features are ordered by their `group_by` class into contiguous blocks and
#' clustered within each block, so a class stays together and reads as a unit;
#' with no `group_by`, clustering is global. A coloured strip along each axis
#' carries the class, using the same palette as [plotLoadings()] and the
#' feature annotation in [plotHeatmap()] - so a class is one colour everywhere
#' in the report.
#'
#' @param correlations The `correlations` data frame from [runStats()], or the
#'   whole list. Long format: `feature_a`, `feature_b`, `estimate`, plus
#'   `correlation` / `subset` / `method` identifying each run.
#' @param variable_meta Feature metadata, or a `DatasetExperiment` to take it
#'   from. Only needed for `group_by`.
#' @param name,subset,method Select one run when the table holds several.
#'   `NULL` requires that only one is present.
#' @param group_by `variable_meta` column used to block and colour features.
#' @param cluster Cluster features within each block (globally when there is no
#'   `group_by`), using average linkage on `1 - |r|`.
#' @param show_values Print the coefficient in each cell. Ignored above 30
#'   features, where the text is unreadable anyway.
#' @return A `ggplot`, or `NULL` if the selection is empty.
#' @examples
#' de <- loadDataset(system.file("extdata", "example_synthetic.RDS",
#'                                package = "MRManalyzeR"))
#' de <- subsetDataset(de, conditions = list(Sample_type = "Sample"))
#' params <- list(correlations = list(enabled = TRUE, entries = list(
#'   list(name = "pairs", methods = "spearman", subsets = list(list())))))
#' st <- runStats(de, params)
#' plotCorrelation(st, de, group_by = "Enzymatic_pathway")
#' @family stats
#' @export
plotCorrelation = function(correlations, variable_meta = NULL,
                            name = NULL, subset = NULL, method = NULL,
                            group_by    = NULL,
                            cluster     = TRUE,
                            show_values = FALSE){

  if(is.list(correlations) && !is.data.frame(correlations) &&
     !is.null(correlations$correlations))
    correlations = correlations$correlations
  cdf = as.data.frame(correlations)
  if(nrow(cdf) == 0) return(NULL)

  pick = function(df, col, val){
    if(is.null(val)){
      u = unique(df[[col]])
      if(length(u) > 1)
        stop(sprintf("[plotCorrelation] several values of '%s'; name one: %s",
                     col, paste(u, collapse = ", ")))
      return(df)
    }
    df[df[[col]] %in% val, , drop = FALSE]
  }
  cdf = pick(cdf, "correlation", name)
  cdf = pick(cdf, "subset",      subset)
  cdf = pick(cdf, "method",      method)
  if(nrow(cdf) == 0) return(NULL)

  ttl = sprintf("%s - %s - %s", cdf$correlation[1], cdf$subset[1],
                cdf$method[1])

  # Long lower triangle -> full square, so ordering and clustering can work on
  # a proper matrix.
  feats = unique(c(cdf$feature_a, cdf$feature_b))
  M = matrix(NA_real_, length(feats), length(feats),
             dimnames = list(feats, feats))
  M[cbind(cdf$feature_a, cdf$feature_b)] = cdf$estimate
  M[cbind(cdf$feature_b, cdf$feature_a)] = cdf$estimate

  if(methods::is(variable_meta, "DatasetExperiment"))
    variable_meta = variable_meta$variable_meta
  vmeta = if(is.null(variable_meta)) NULL else as.data.frame(variable_meta)

  have_groups = !is.null(group_by) && !is.null(vmeta) &&
    group_by %in% colnames(vmeta) && "Compound" %in% colnames(vmeta)

  clust_within = function(fs){
    if(isTRUE(cluster) && length(fs) > 2){
      # A pair with no estimate counts as unrelated, so one gap does not
      # stop the whole block from being clustered.
      D = 1 - abs(M[fs, fs, drop = FALSE])
      D[is.na(D)] = 1
      d = tryCatch(stats::as.dist(D), error = function(e) NULL)
      h = if(!is.null(d)) tryCatch(stats::hclust(d, method = "average"),
                                   error = function(e) NULL) else NULL
      if(!is.null(h)) fs = fs[h$order]
    }
    fs
  }

  if(have_groups){
    g = as.character(vmeta[[group_by]][match(feats, vmeta$Compound)])
    g[is.na(g)] = "NA"
    feats = unlist(lapply(unique(g),
                          function(k) clust_within(feats[g == k])),
                   use.names = FALSE)
    grp_vals = as.character(vmeta[[group_by]][match(feats, vmeta$Compound)])
    grp_vals[is.na(grp_vals)] = "NA"
  } else {
    feats = clust_within(feats)
    grp_vals = rep("NA", length(feats))
  }

  pal = if(have_groups) .group_fill_colors(grp_vals) else NULL
  feat_cols = if(have_groups) unname(pal[grp_vals])
              else rep("grey70", length(feats))
  names(feat_cols) = feats

  # A dedicated axis level carries the class colour blocks, so the axis text
  # itself can stay black and readable.
  strip = "__grp__"
  lev = if(have_groups) c(strip, rev(feats)) else rev(feats)

  long = expand.grid(feature_a = feats, feature_b = feats,
                     stringsAsFactors = FALSE)
  long$r = M[cbind(long$feature_a, long$feature_b)]
  ia = match(long$feature_a, feats); ib = match(long$feature_b, feats)
  long = long[ia > ib | (ia == ib & !is.na(long$r)), , drop = FALSE]
  long$feature_a = factor(long$feature_a, levels = lev)
  long$feature_b = factor(long$feature_b, levels = lev)

  blank = function(l) ifelse(l == strip, "", l)

  p = ggplot2::ggplot(long, ggplot2::aes(.data$feature_a, .data$feature_b,
                                         fill = .data$r)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.1) +
    ggplot2::scale_fill_gradientn(colours = .diverging_pal(101),
                                  limits = c(-1, 1), na.value = "grey90") +
    ggplot2::scale_x_discrete(labels = blank, drop = FALSE,
                              na.translate = FALSE) +
    ggplot2::scale_y_discrete(labels = blank, drop = FALSE,
                              na.translate = FALSE) +
    ggplot2::coord_fixed() +
    ggplot2::theme_classic(base_size = 10) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5,
                                          hjust = 1, size = 6),
      axis.text.y = ggplot2::element_text(size = 6),
      axis.ticks  = ggplot2::element_blank(),
      legend.position = "right") +
    ggplot2::labs(x = NULL, y = NULL, fill = cdf$method[1], title = ttl)

  if(have_groups){
    x_strip = data.frame(feature_a = factor(feats, levels = lev),
                         feature_b = factor(strip, levels = lev))
    y_strip = data.frame(feature_a = factor(strip, levels = lev),
                         feature_b = factor(feats, levels = lev))
    # Invisible points carry only the discrete scale so the class legend
    # renders. Numeric NA is deliberate: a factor NA would make the discrete
    # axes draw a phantom "NA" row and column.
    legend_pts = data.frame(x = NA_real_, y = NA_real_,
                            grp = factor(names(pal), levels = names(pal)))
    p = p +
      ggplot2::geom_tile(data = x_strip,
                         ggplot2::aes(.data$feature_a, .data$feature_b),
                         fill = feat_cols[feats], inherit.aes = FALSE) +
      ggplot2::geom_tile(data = y_strip,
                         ggplot2::aes(.data$feature_a, .data$feature_b),
                         fill = feat_cols[feats], inherit.aes = FALSE) +
      ggplot2::geom_point(data = legend_pts,
                          ggplot2::aes(.data$x, .data$y, colour = .data$grp),
                          inherit.aes = FALSE, na.rm = TRUE, size = 3) +
      ggplot2::scale_colour_manual(values = pal, name = group_by,
                                   drop = FALSE) +
      ggplot2::guides(colour = ggplot2::guide_legend(
        override.aes = list(size = 3, shape = 15)))
  }

  if(isTRUE(show_values) && length(feats) <= 30)
    p = p + ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", .data$r)),
                               size = 2)

  p
}
