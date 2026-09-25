#' Group-mean heatmap that shows which samples produced the mean
#'
#' The same picture [plotHeatmap()] draws - compounds on rows, samples on
#' columns, a coloured annotation bar for the sample factor along the top and
#' for the compound class down the left, gaps between blocks - with one thing
#' changed: what a cell's colour means. The **hue** is the group's mean scaled
#' value, and the **opacity** of each sample's cell is how much that sample
#' contributed to it.
#'
#' A block that is uniformly solid is an effect the whole group shares. A block
#' carrying one opaque cell among faded ones is a mean resting on a single
#' sample - which is what a plain group-mean heatmap hides, and what a reader
#' needs in order to judge whether a colour contrast is worth believing.
#'
#' [plotHeatmap()] never averages, so it does not have this problem, and it is
#' the better choice while the sample axis stays legible - in practice, while
#' there are two groups. Past three the cells get narrow and the eye stops
#' resolving individual samples anyway, which is where this earns its place.
#'
#' @section How a cell is built:
#' Values are transformed, then scaled per compound across **all** samples, so
#' a cell reads "high or low for this compound" rather than "abundant or not".
#' Rows are therefore not comparable to each other in absolute terms.
#'
#' The hue is the mean of those scores within a group. The opacity is the
#' sample's own score multiplied by the sign of its group mean, so a sample
#' that moved with its group is opaque and one that sat near zero, or pulled
#' the other way, fades. `alpha_floor` stops a non-contributing sample becoming
#' invisible and reading as missing data; genuinely missing values are drawn in
#' `na.value` grey at full opacity, so the two cannot be confused.
#'
#' Note that "contributed" is measured against the mean of **all** samples, not
#' against the group's own members, so it is not the same as "has the largest
#' value". A group of four holding one very high sample and three low ones has
#' a *negative* mean if the three sit further from the overall centre than the
#' one sits above it - and then the three are what produced that mean, and the
#' high sample is the one opposing it and fading out. The opacity answers "did
#' this sample push its group where the group ended up", not "is this sample
#' large".
#'
#' @section Choosing the colour limit:
#' `cap` defaults to a quantile of the group means actually being drawn, and
#' that default matters more than it looks. Averaging *n* replicates shrinks a
#' scaled value by roughly `sqrt(n)`, so a limit chosen for per-sample scores -
#' the `+/-2.5` that suits [plotHeatmap()] - leaves a group-mean panel using a
#' small part of the colour range and looking washed out, whatever the opacity
#' is set to. Take the limit from the quantity being drawn.
#'
#' @param de A `struct::DatasetExperiment`.
#' @param color_samples_by `sample_meta` column defining the groups to average
#'   within, and the annotation bar along the top. Named to match
#'   [plotHeatmap()], where it does the annotation half of the same job.
#' @param features Compounds to include; `NULL` uses all of them.
#' @param group_features_by `variable_meta` column used to block compounds and
#'   colour the annotation bar down the left; `NULL` leaves them unblocked.
#' @param transform `"none"`, `"log2"` or `"sqrt"`, applied before scaling.
#'   Negative values are floored at 0 first.
#' @param cap Colour limit, in standard deviations of the group means. `NULL`
#'   takes the 97th percentile of `abs()` of the group means, so a few cells
#'   saturate and the rest get the full range.
#' @param alpha_floor Opacity given to a sample that did not contribute to its
#'   group's mean. Below about 0.4 those cells read as missing data.
#' @param feature_labels `"show"`, `"small"` or `"hide"`. With
#'   `group_features_by` set, the names are drawn between the class colour bar
#'   and the cells - block, then names, then data, which is how a pathway
#'   figure is read. Without it they sit on the axis as usual.
#' @param sample_labels Draw the sample axis. `FALSE` for wide cohorts, where
#'   the samples stay individually visible but not individually named.
#' @param sample_id_head `sample_meta` column supplying sample labels; `NULL`
#'   uses row names.
#' @param title Plot title.
#' @return A `ggplot`, or `NULL` when there is too little data to draw (fewer
#'   than 2 samples, or fewer than 2 compounds with any variance).
#' @examples
#' de <- loadDataset(system.file("extdata", "example_synthetic.RDS",
#'                                package = "MRManalyzeR"))
#' de <- subsetDataset(de, conditions = list(Sample_type = "Sample"))
#' plotGroupHeatmap(de, color_samples_by = "Treatment",
#'                    group_features_by = "Enzymatic_pathway")
#' @family stats
#' @seealso [plotHeatmap()], the same layout without the averaging.
#' @export
plotGroupHeatmap = function(de, color_samples_by,
                              features          = NULL,
                              group_features_by = NULL,
                              transform         = c("none", "log2", "sqrt"),
                              cap               = NULL,
                              alpha_floor       = 0.6,
                              feature_labels    = c("show", "small", "hide"),
                              sample_labels     = TRUE,
                              sample_id_head    = NULL,
                              title             = NULL){

  transform      = match.arg(transform)
  feature_labels = match.arg(feature_labels)

  # Three fill scales in one plot - the values, the sample annotation and the
  # compound annotation. ggplot allows one per plot; ggnewscale is what lets
  # the annotation bars use the same palette as plotHeatmap() instead of being
  # demoted to facet strips.
  if(!requireNamespace("ggnewscale", quietly = TRUE)){
    warning("[plotGroupHeatmap] ggnewscale is not installed.")
    return(NULL)
  }

  smeta = as.data.frame(de$sample_meta)
  vmeta = as.data.frame(de$variable_meta)
  M     = as.matrix(as.data.frame(de$data))
  if(!is.null(features)) M = M[, intersect(features, colnames(M)),
                               drop = FALSE]
  if(nrow(M) < 2 || ncol(M) < 2) return(NULL)

  color_samples_by = .resolve_meta_col(color_samples_by, smeta) %||%
    stop(sprintf(
      "[plotGroupHeatmap] '%s' is not a sample_meta column. Present: %s.",
      color_samples_by, paste(colnames(smeta), collapse = ", ")))
  group_features_by = .resolve_meta_col(group_features_by, vmeta)
  sample_id_head    = .resolve_meta_col(sample_id_head, smeta)

  if(!is.null(sample_id_head))
    rownames(M) = as.character(smeta[[sample_id_head]])

  if(transform == "log2"){
    M[M < 0 & !is.na(M)] = 0; M = log2(M + 1)
  } else if(transform == "sqrt"){
    M[M < 0 & !is.na(M)] = 0; M = sqrt(M)
  }

  # A zero-variance compound has no scaled value at all - it would be a row of
  # NaN rather than a row of zeros - so drop it instead of drawing it grey.
  vars = apply(M, 2, stats::var, na.rm = TRUE)
  M    = M[, is.finite(vars) & vars > 0, drop = FALSE]
  if(ncol(M) < 2) return(NULL)

  # Scaled per compound across ALL samples, not within group: scaling within
  # group would remove the very differences the figure is about.
  Z   = scale(M, center = TRUE, scale = TRUE)
  grp = factor(as.character(smeta[[color_samples_by]]))

  gm = vapply(levels(grp),
              function(g) colMeans(Z[grp == g, , drop = FALSE], na.rm = TRUE),
              numeric(ncol(Z)))
  dimnames(gm) = list(colnames(Z), levels(grp))

  # Each sample's contribution to its own group's mean: its score signed
  # towards that mean. A group mean of exactly zero gives sign() == 0, which
  # correctly leaves every sample in it with no contribution to make.
  contrib = Z
  for(g in levels(grp)){
    i = which(grp == g)
    contrib[i, ] = sweep(Z[i, , drop = FALSE], 2, sign(gm[, g]), `*`)
  }

  a_ref = stats::quantile(abs(Z), 0.95, na.rm = TRUE)
  if(!is.finite(a_ref) || a_ref <= 0) a_ref = 1
  alpha_floor = min(max(alpha_floor, 0), 1)
  alpha = pmin(pmax(contrib, 0) / as.numeric(a_ref), 1)
  alpha = alpha_floor + (1 - alpha_floor) * alpha

  if(is.null(cap)) cap = stats::quantile(abs(gm), 0.97, na.rm = TRUE)
  cap = as.numeric(cap)
  if(!is.finite(cap) || cap <= 0) cap = 1

  # ---- axis order, with a spacer level between blocks ----------------------
  # pheatmap separates blocks with a gap; here that is an axis level nothing is
  # drawn at, leaving the panel background showing through. Spacers are named
  # ".gap*" and blanked by the axis labeller, as is the ".ann" level carrying
  # the compound annotation bar.
  feats = colnames(Z)
  cls   = if(is.null(group_features_by)) rep("", length(feats))
          else as.character(vmeta[[group_features_by]][
            match(feats, vmeta$Compound)])
  cls[is.na(cls)] = "NA"
  ford  = order(match(cls, unique(cls)))
  feats = feats[ford]; cls = cls[ford]

  sord  = order(grp)
  samps = rownames(Z)[sord]
  sgrp  = as.character(grp)[sord]

  .with_gaps = function(x, blocks, tag){
    out = character(0); n = 0L
    for(b in unique(blocks)){
      if(n > 0L){ n = n + 1L; out = c(out, paste0(".gap", tag, n)) }
      out = c(out, x[blocks == b])
    }
    out
  }
  feat_disp = .with_gaps(feats, cls,  "F")
  samp_disp = .with_gaps(samps, sgrp, "S")

  # y is drawn bottom-up, so the display order is reversed and the sample
  # annotation bar becomes the topmost level.
  # Names between the class bar and the cells: reserve columns for them on
  # the x axis, then draw them as text rather than as axis labels. How many
  # columns is a guess from the longest name, because how many cells a name
  # spans depends on the panel size, which is not known until it is drawn.
  # Without a class bar there is nothing to sit inside of, so the axis keeps
  # them.
  lab_in = !is.null(group_features_by) && feature_labels != "hide"
  n_lab  = if(!lab_in) 0L else
    as.integer(min(max(ceiling(max(nchar(feats), 1L) / 6), 1), 20))

  y_lv = c(rev(feat_disp), ".gapTop", ".ann")
  x_lv = if(is.null(group_features_by)) samp_disp
         else c(".ann", if(n_lab) paste0(".lab", seq_len(n_lab)),
                ".gapLeft", samp_disp)

  fac = function(v, lv) factor(v, levels = lv)

  # ---- the cells ----------------------------------------------------------
  d = expand.grid(sample = samps, feature = feats, stringsAsFactors = FALSE)
  d$group = sgrp[match(d$sample, samps)]
  d$fill  = pmin(pmax(gm[cbind(d$feature, d$group)], -cap), cap)
  d$alpha = alpha[cbind(d$sample, d$feature)]

  # A missing value is grey at full opacity, so it cannot be mistaken for a
  # sample that simply did not contribute.
  na_cell = is.na(Z[cbind(d$sample, d$feature)])
  d$fill[na_cell]  = NA_real_
  d$alpha[na_cell] = 1
  d$sample  = fac(d$sample,  x_lv)
  d$feature = fac(d$feature, y_lv)

  # .group_fill_colors() keys on sorted levels, so a class gets the same colour
  # here as it does in plotHeatmap(), plotLoadings() and plotCorrelation().
  s_cols = .group_fill_colors(sgrp)
  top    = data.frame(sample  = fac(samps, x_lv),
                      feature = fac(".ann", y_lv),
                      ann     = factor(sgrp, levels = names(s_cols)),
                      stringsAsFactors = FALSE)

  fs = switch(feature_labels, show = 9, small = 6, 6)
  blank_dot = function(l) ifelse(substr(l, 1, 1) == ".", "", l)

  p = ggplot2::ggplot() +
    ggplot2::geom_tile(
      data = d,
      ggplot2::aes(.data$sample, .data$feature, fill = .data$fill,
                   alpha = .data$alpha)) +
    ggplot2::scale_fill_gradientn(colours = .diverging_pal(101),
                                  limits = c(-cap, cap), na.value = "grey88",
                                  name = "group mean") +
    ggplot2::scale_alpha_identity() +
    ggnewscale::new_scale_fill() +
    ggplot2::geom_tile(
      data = top,
      ggplot2::aes(.data$sample, .data$feature, fill = .data$ann)) +
    ggplot2::scale_fill_manual(values = s_cols, name = color_samples_by,
                               drop = FALSE)

  if(!is.null(group_features_by)){
    f_cols = .group_fill_colors(cls)
    left   = data.frame(sample  = fac(".ann", x_lv),
                        feature = fac(feats, y_lv),
                        ann     = factor(cls, levels = names(f_cols)),
                        stringsAsFactors = FALSE)
    p = p + ggnewscale::new_scale_fill() +
      ggplot2::geom_tile(
        data = left,
        ggplot2::aes(.data$sample, .data$feature, fill = .data$ann)) +
      ggplot2::scale_fill_manual(values = f_cols, name = group_features_by,
                                 drop = FALSE)
  }

  # Drawn at the left edge of the spacer column, running leftwards into the
  # reserved columns, so the names end where the cells begin.
  if(lab_in)
    p = p + ggplot2::geom_text(
      data = data.frame(sample  = fac(".gapLeft", x_lv),
                        feature = fac(feats, y_lv),
                        label   = feats, stringsAsFactors = FALSE),
      ggplot2::aes(.data$sample, .data$feature, label = .data$label),
      hjust = 1, nudge_x = -0.5, size = fs / ggplot2::.pt,
      inherit.aes = FALSE)

  p +
    ggplot2::scale_x_discrete(labels = blank_dot, drop = FALSE) +
    ggplot2::scale_y_discrete(labels = blank_dot, drop = FALSE) +
    ggplot2::labs(x = NULL, y = NULL, title = title) +
    .mrm_theme() +
    ggplot2::theme(
      legend.position = "right",
      axis.text.x = if(isTRUE(sample_labels))
        ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6)
      else ggplot2::element_blank(),
      axis.text.y = if(feature_labels == "hide" || lab_in)
        ggplot2::element_blank() else ggplot2::element_text(size = fs),
      axis.ticks  = ggplot2::element_blank(),
      axis.line   = ggplot2::element_blank())
}
