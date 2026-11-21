#' scBrailleR: Accessible Text Descriptions of Single-Cell Plots
#'
#' Extends BrailleR with `VI()` and `Augment()` methods for plots
#' commonly produced in single-cell RNA-seq analysis. Currently
#' provides accessible wrappers for [Seurat::ElbowPlot()] and
#' [Seurat::DimPlot()].
#'
#' @keywords internal
#' @importFrom BrailleR VI Augment WhereXY
#' @importFrom stats dist mad median quantile
#' @importFrom utils data
"_PACKAGE"
