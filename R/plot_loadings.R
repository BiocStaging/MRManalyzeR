#' Plot PCA loadings as ranked bars, one panel per component
#'
#' The scores plot says *whether* samples separate; the loadings say *which
#' compounds* drove it. Each feature's contribution to a component is drawn as
#' a bar, sorted within its colour group so the compounds pulling hardest in
#' each direction sit at the ends.
#'
#' Colouring by a `variable_meta` column - typically the enzymatic pathway -
#' is what turns this from a ranked list into something interpretable: if the
#' bars at one end of PC1 are predominantly one pathway, the separation in the
#' scores plot has a biochemical reading rather than just a statistical one.
#' The palette is shared with the heatmap feature annotation, so a class is the
#' same colour in both figures.
#'
#' @param pca The list returned by [mrmRunPCA()].
#' @param variable_meta Feature metadata, or a `DatasetExperiment` to take it
#'   from. Matched to the loadings by `Compound`.
#' @param components Length-2 integer vector: which components to show. Falls
#'   back to the first two if the fit has fewer.
#' @param colour_by `variable_meta` column used to colour and group the bars;
#'   `NULL` draws every bar in one colour.
#' @return A `ggplot`, or `NULL` if the PCA could not be fit.
#' @examples
#' de  <- loadDataset(system.file("extdata", "example_synthetic.RDS",
#'                                 package = "MRManalyzeR"))
#' pca <- mrmRunPCA(de, transform = "log2")
#' plotLoadings(pca, de, colour_by = "Enzymatic_pathway")
#' @family QC check
#' @seealso [mrmPlotPCA()] for the matching scores plot.
#' @export
plotLoadings <- function(pca, variable_meta, components = c(1, 2),
                         colour_by = NULL){

  if(is.null(pca$pr)) return(NULL)

  if(methods::is(variable_meta, "DatasetExperiment"))
    variable_meta <- variable_meta$variable_meta
  vmeta <- as.data.frame(variable_meta)

  pr <- pca$pr
  pc <- as.integer(components)[c(1L, 2L)]
  if(any(pc > ncol(pr$x))) pc <- c(1L, min(2L, ncol(pr$x)))
  pc_labs <- sprintf("PC%d", pc)

  feats <- rownames(pr$rotation)

  colour_by <- .resolve_meta_col(colour_by, vmeta)
  has_colour <- !is.null(colour_by)
  colour_vec <- if(has_colour){
    cv <- as.character(vmeta[[colour_by]][match(feats, vmeta$Compound)])
    cv[is.na(cv)] <- "NA"
    cv
  } else {
    rep("(all)", length(feats))
  }
  groups_in_order <- unique(colour_vec)

  # Bars are ordered within colour group, so each class reads as a block whose
  # strongest contributors are at its edges.
  build_pc <- function(loading, pc_lab){
    ord <- unlist(lapply(groups_in_order, function(g){
      idx <- which(colour_vec == g)
      if(!length(idx)) return(character(0))
      feats[idx[order(loading[idx], decreasing = TRUE)]]
    }))
    list(df = data.frame(feature = feats,
                         feat_pc = paste0(feats, "__", pc_lab),
                         loading = loading,
                         colour  = factor(colour_vec,
                                          levels = groups_in_order),
                         PC      = pc_lab,
                         stringsAsFactors = FALSE),
         levels = rev(paste0(ord, "__", pc_lab)))
  }

  a <- build_pc(pr$rotation[, pc[1]], pc_labs[1])
  b <- build_pc(pr$rotation[, pc[2]], pc_labs[2])

  long <- rbind(a$df, b$df)
  long$feat_pc <- factor(long$feat_pc, levels = c(a$levels, b$levels))
  long$PC     <- factor(long$PC, levels = pc_labs)

  ggplot2::ggplot(long, ggplot2::aes(.data$loading, .data$feat_pc,
                                     fill = .data$colour)) +
    ggplot2::geom_col() +
    ggplot2::geom_vline(xintercept = 0, color = "grey50", linewidth = 0.3) +
    ggplot2::facet_wrap(~ PC, scales = "free", nrow = 1) +
    ggplot2::scale_y_discrete(labels = function(x) sub("__[^_]+$", "", x)) +
    ggplot2::scale_fill_manual(values = .group_fill_colors(colour_vec),
                               na.value = "grey60", drop = FALSE) +
    ggplot2::labs(x = "Loading", y = NULL,
                  fill = if(has_colour) colour_by else "") +
    ggplot2::theme_classic(base_size = 10) +
    ggplot2::theme(axis.text.y     = ggplot2::element_text(size = 6),
                   strip.text      = ggplot2::element_text(face = "bold"),
                   legend.position = "bottom") +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 3, byrow = TRUE))
}
