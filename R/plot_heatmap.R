#' Sample-by-feature heatmap with class annotations
#'
#' A z-scored heatmap of samples against compounds, with samples blocked by a
#' study factor and compounds blocked by their class. Where the comparisons
#' test one compound at a time, this shows the whole panel at once - which is
#' how you see that a treatment moved an entire pathway rather than a handful
#' of compounds that happened to clear a threshold.
#'
#' Values are transformed, then scaled per compound (`scale = "z"`) so that
#' compounds spanning three orders of magnitude are comparable in one picture -
#' without it the abundant compounds dominate and everything else is a flat
#' wash. Scores are capped at `cap` standard deviations so one outlier cannot
#' consume the whole colour range. Zero-variance compounds are dropped.
#' Missing values are drawn grey; only for ordering and clustering do they
#' count as the compound's mean (0 after scaling).
#'
#' Compounds are ordered into contiguous blocks by `group_features_by`, with
#' gaps drawn between blocks, and clustered *within* each block when
#' `cluster_within_groups` is set - so a class stays together and the
#' clustering tells you about structure inside it, not about which classes
#' happen to resemble each other. The class palette matches [plotLoadings()]
#' and [plotCorrelation()].
#'
#' `orientation` decides which axis carries samples: `"samples_y"` puts them on
#' rows, `"samples_x"` transposes, and `"auto"` puts whichever dimension is
#' larger on the y-axis so a wide panel does not become tall and unreadable.
#'
#' @param de A `struct::DatasetExperiment`.
#' @param features Compounds to include; `NULL` uses all of them.
#' @param title Plot title.
#' @param color_samples_by `sample_meta` column used to block and annotate
#'   samples.
#' @param group_features_by `variable_meta` column used to block and annotate
#'   compounds; `NULL` leaves them unblocked.
#' @param transform `"none"`, `"log2"` or `"sqrt"`, applied before scaling.
#'   Negative values are floored at 0 first.
#' @param scale `"z"` to scale each compound, or `"none"`.
#' @param cluster_rows Cluster the sample axis.
#' @param cluster_within_groups Cluster compounds inside each class block.
#' @param orientation `"samples_y"`, `"samples_x"` or `"auto"`.
#' @param feature_labels `"show"`, `"small"` or `"hide"` - the compound axis
#'   only; sample labels always stay readable.
#' @param sample_id_head `sample_meta` column supplying sample labels; `NULL`
#'   uses row names.
#' @param cap Colour-scale limit in standard deviations.
#' @param cell_size Cell size in points: one value for square cells, or
#'   `c(width, height)`. Because these are absolute units they survive being
#'   stretched to fill a canvas, which a proportional cell does not - so leave
#'   this set unless you want cells to take whatever aspect the figure
#'   dimensions imply. `NULL` restores that fill-the-canvas behaviour.
#' @return A `ggplot` wrapping the heatmap grob, or `NULL` when there is too
#'   little data to draw (fewer than 3 samples, 2 features, or no variable
#'   compound).
#' @examples
#' de <- loadDataset(system.file("extdata", "example_synthetic.RDS",
#'                                package = "MRManalyzeR"))
#' de <- subsetDataset(de, conditions = list(Sample_type = "Sample"))
#' plotHeatmap(de, color_samples_by = "Treatment",
#'              group_features_by = "Enzymatic_pathway",
#'              title = "All features")
#' @family stats
#' @export
plotHeatmap = function(de, features = NULL, title = NULL,
                        color_samples_by,
                        group_features_by     = NULL,
                        transform             = c("none", "log2", "sqrt"),
                        scale                 = c("z", "none"),
                        cluster_rows          = FALSE,
                        cluster_within_groups = TRUE,
                        orientation           = c("samples_y", "samples_x",
                                                  "auto"),
                        feature_labels        = c("show", "small", "hide"),
                        sample_id_head        = NULL,
                        cap                   = 2.5,
                        cell_size             = 10){

  transform      = match.arg(transform)
  scale          = match.arg(scale)
  orientation    = match.arg(orientation)
  feature_labels = match.arg(feature_labels)

  if(!requireNamespace("pheatmap", quietly = TRUE)){
    warning("[plotHeatmap] pheatmap is not installed.")
    return(NULL)
  }

  smeta = as.data.frame(de$sample_meta)
  vmeta = as.data.frame(de$variable_meta)
  M     = as.matrix(as.data.frame(de$data))
  if(!is.null(features)) M = M[, intersect(features, colnames(M)),
                               drop = FALSE]

  if(nrow(M) < 3 || ncol(M) < 2) return(NULL)
  color_samples_by = .resolve_meta_col(color_samples_by, smeta) %||%
    stop(sprintf("[plotHeatmap] '%s' is not a sample_meta column. Present: %s.",
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

  vars = apply(M, 2, stats::var, na.rm = TRUE)
  M    = M[, is.finite(vars) & vars > 0, drop = FALSE]
  if(ncol(M) == 0) return(NULL)

  if(scale == "z") M = scale(M, center = TRUE, scale = TRUE)
  M[M >  cap] =  cap
  M[M < -cap] = -cap

  ann_row = data.frame(row.names = rownames(M))
  ann_row[[color_samples_by]] = smeta[[color_samples_by]]

  ann_col = NULL; gaps_col = NULL
  if(!is.null(group_features_by)){
    feats = colnames(M)
    g = as.character(vmeta[[group_features_by]][match(feats,
                                                      vmeta$Compound)])
    g[is.na(g)] = "NA"
    levels_in_order = unique(g)

    ordered = unlist(lapply(levels_in_order, function(k){
      in_k = feats[g == k]
      if(isTRUE(cluster_within_groups) && length(in_k) > 1){
        Mk = M[, in_k, drop = FALSE]
        Mk[is.na(Mk)] = 0
        d = stats::dist(t(Mk))
        if(all(is.finite(d)))
          in_k = in_k[stats::hclust(d, method = "average")$order]
      }
      in_k
    }), use.names = FALSE)

    M = M[, ordered, drop = FALSE]
    g2 = as.character(vmeta[[group_features_by]][match(ordered,
                                                       vmeta$Compound)])
    g2[is.na(g2)] = "NA"
    ann_col = data.frame(row.names = ordered)
    ann_col[[group_features_by]] = factor(g2, levels = levels_in_order)
    gaps_col = utils::head(cumsum(rle(g2)$lengths), -1)
  }

  ord = order(as.character(ann_row[[color_samples_by]]))
  M       = M[ord, , drop = FALSE]
  ann_row = ann_row[ord, , drop = FALSE]
  gaps_row = utils::head(
    cumsum(rle(as.character(ann_row[[color_samples_by]]))$lengths), -1)

  ann_colors = c(.ann_colours(ann_row), .ann_colours(ann_col))

  # Samples are clustered on a copy with gaps at the compound mean, so a gap
  # cannot break the clustering; in the picture itself gaps stay grey.
  hc = FALSE
  if(isTRUE(cluster_rows) && nrow(M) > 1){
    M_fill = M
    M_fill[is.na(M_fill)] = 0
    hc = stats::hclust(stats::dist(M_fill), method = "complete")
  }

  flip = identical(orientation, "samples_x") ||
    (identical(orientation, "auto") && nrow(M) > ncol(M))

  if(flip){
    M_use = t(M)
    ann_row_use = ann_col;  ann_col_use = ann_row
    gaps_row_use = gaps_col; gaps_col_use = gaps_row
    cluster_rows_use = FALSE; cluster_cols_use = hc
  } else {
    M_use = M
    ann_row_use = ann_row;  ann_col_use = ann_col
    gaps_row_use = gaps_row; gaps_col_use = gaps_col
    cluster_rows_use = hc; cluster_cols_use = FALSE
  }

  feat_fs   = switch(feature_labels, show = 9, small = 4, 4)
  feat_show = feature_labels != "hide"
  if(flip){   # features on rows, samples on columns
    show_rn = feat_show; fs_row = feat_fs; show_cn = TRUE;      fs_col = 9
  } else {    # features on columns, samples on rows
    show_cn = feat_show; fs_col = feat_fs; show_rn = TRUE;      fs_row = 9
  }

  ph = pheatmap::pheatmap(
    M_use,
    cluster_rows      = cluster_rows_use,
    cluster_cols      = cluster_cols_use,
    gaps_row          = if(length(gaps_row_use)) gaps_row_use else NULL,
    gaps_col          = if(length(gaps_col_use)) gaps_col_use else NULL,
    annotation_row    = ann_row_use,
    annotation_col    = ann_col_use,
    annotation_colors = if(length(ann_colors)) ann_colors else NA,
    color  = .diverging_pal(101),
    breaks = seq(-cap, cap, length.out = 102),
    # pheatmap wants NA, not NULL, for "size it to the canvas". A length-1
    # cell_size gives square cells; length-2 is c(width, height).
    cellwidth     = .cell_dim(cell_size[1]),
    cellheight    = .cell_dim(cell_size[min(2, length(cell_size))]),
    show_rownames = show_rn, show_colnames = show_cn,
    fontsize_row  = fs_row,  fontsize_col  = fs_col,
    border_color  = NA,
    na_col        = "grey85",
    legend_breaks = c(-cap, -1, 0, 1, cap),
    legend_labels = c(sprintf("<= -%g", cap), "-1", "0", "1",
                      sprintf(">= %g", cap)),
    # pheatmap's own default is NA, and it tests `is.na(main)` - so passing the
    # documented title = NULL default straight through fails with "argument is
    # of length zero" rather than drawing an untitled heatmap.
    main   = title %||% NA,
    silent = TRUE)

  # Wrapping the pheatmap gtable in a ggplot makes it fill the knitr canvas and
  # be captured in results="asis" loops. print(ph$gtable) emits a text
  # description, grid.draw() is not captured, and cowplot::plot_grid() places
  # the grob at natural size leaving whitespace.
  ggplot2::ggplot() +
    ggplot2::annotation_custom(ph$gtable, xmin = -Inf, xmax = Inf,
                               ymin = -Inf, ymax = Inf) +
    ggplot2::theme_void()
}

#' Coerce a cell size to what pheatmap expects
#'
#' `NULL`, `NA` or a non-finite value all mean "let pheatmap size the cell to
#' the canvas", which pheatmap spells `NA`.
#'
#' @param x Cell edge in points, or `NULL`.
#' @return A finite numeric, or `NA`.
#' @keywords internal
#' @noRd
.cell_dim = function(x){
  if(is.null(x) || length(x) == 0) return(NA)
  x = suppressWarnings(as.numeric(x[1]))
  if(!is.finite(x) || x <= 0) return(NA)
  x
}
