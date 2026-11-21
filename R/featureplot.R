### scBrailleR: accessible description of a FeaturePlot
###
### A FeaturePlot is the DimPlot scatter coloured by one gene's
### expression. A sighted user reads, per gene:
###   1. how many cells express it at all;
###   2. where in the embedding expression is high (the "hotspot");
###   3. which cluster that region belongs to.
### (1) and (3) are answered from the expression vector and the group
### labels directly; (2) with a mean-expression grid aligned to the
### DimPlot grids, so the two descriptions can be read side by side.

#' Accessible FeaturePlot for Seurat or SingleCellExperiment
#'
#' Wraps [Seurat::FeaturePlot()] (Seurat) or [scater::plotReducedDim()]
#' coloured by expression (SingleCellExperiment), and attaches the
#' embedding, the expression of each feature and the group labels so
#' that [BrailleR::VI()] can describe where each gene is expressed.
#'
#' @param object A Seurat or SingleCellExperiment object with at least
#'   one reduction computed.
#' @param features Character vector of genes (or other features) to plot.
#' @param reduction Reduction name (default: Seurat's default reduction,
#'   or the first `reducedDimNames()` entry for SCE).
#' @param group.by Grouping used to say which cluster expresses each gene.
#'   Defaults to active idents (Seurat) or `colLabels()` (SCE, if present;
#'   otherwise no per-group facts are reported).
#' @param dims Which two dimensions of the reduction to plot.
#' @param layer Expression values to describe: a Seurat layer (default
#'   `"data"`) or an SCE assay (default `"logcounts"`).
#' @param ... Passed to [Seurat::FeaturePlot()] or [scater::plotReducedDim()].
#' @return A `ggplot` (a `patchwork` for several features) with extra
#'   class `"sc_featureplot"`.
#' @examples
#' \donttest{
#' if (requireNamespace("SeuratObject", quietly = TRUE) &&
#'     requireNamespace("Seurat", quietly = TRUE)) {
#'   data("pbmc_small", package = "SeuratObject")
#'   p <- BR_FeaturePlot(pbmc_small, c("CD3E", "LYZ"), reduction = "tsne")
#'   BrailleR::VI(p)
#'   BrailleR::VI(p, Describe = TRUE)   # adds a per-group table per gene
#' }
#' }
#' @export
BR_FeaturePlot <- function(object,
                           features,
                           reduction = NULL,
                           group.by  = NULL,
                           dims      = c(1, 2),
                           layer     = NULL,
                           ...) {
  if (inherits(object, "SingleCellExperiment")) {
    if (!requireNamespace("scater", quietly = TRUE)) {
      stop("scater is required for BR_FeaturePlot() on a SingleCellExperiment")
    }
    reduction <- reduction %||% SingleCellExperiment::reducedDimNames(object)[1]
    if (is.null(reduction)) stop("object has no reducedDims")
    layer <- layer %||% "logcounts"
    missing_f <- setdiff(features, rownames(object))
    if (length(missing_f)) stop("features not found: ", paste(missing_f, collapse = ", "))
    plots <- lapply(features, function(f)
      scater::plotReducedDim(object, dimred = reduction, ncomponents = dims,
                             colour_by = f, by_exprs_values = layer, ...))
    p <- if (length(plots) == 1) plots[[1]] else patchwork::wrap_plots(plots)
    emb  <- SingleCellExperiment::reducedDim(object, reduction)[, dims, drop = FALSE]
    expr <- t(as.matrix(SummarizedExperiment::assay(object, layer)[features, , drop = FALSE]))
    cd <- SummarizedExperiment::colData(object)
    if (is.null(group.by) && "label" %in% colnames(cd)) group.by <- "label"
    ids <- if (is.null(group.by)) NULL else cd[[group.by]]
  } else if (inherits(object, "Seurat")) {
    if (!requireNamespace("Seurat", quietly = TRUE)) {
      stop("Seurat is required for BR_FeaturePlot()")
    }
    reduction <- reduction %||% Seurat::DefaultDimReduc(object)
    layer <- layer %||% "data"
    p <- Seurat::FeaturePlot(object, features = features,
                             reduction = reduction, dims = dims, ...)
    emb  <- Seurat::Embeddings(object, reduction = reduction)[, dims, drop = FALSE]
    expr <- as.matrix(Seurat::FetchData(object, vars = features, layer = layer))
    ids  <- if (is.null(group.by)) Seurat::Idents(object) else object[[group.by, drop = TRUE]]
  } else {
    stop("BR_FeaturePlot() supports Seurat and SingleCellExperiment objects, not ",
         class(object)[1])
  }
  if (is.null(colnames(emb))) colnames(emb) <- paste0(reduction, "_", dims)

  attr(p, "embedding") <- emb
  attr(p, "expr")      <- expr
  attr(p, "idents")    <- if (is.null(ids)) NULL else as.factor(ids)
  attr(p, "reduction") <- reduction
  attr(p, "dims")      <- dims
  attr(p, "group_by")  <- group.by
  attr(p, "layer")     <- layer
  class(p) <- c("sc_featureplot", class(p))
  p
}

# Mean of v per (x-bin, y-bin), same orientation as .modal_grid:
# rows = y descending, cols = x ascending. Empty bins are NA.
.mean_grid <- function(emb, v, grid = c(4, 4)) {
  xi <- cut(emb[, 1], .even_breaks(emb[, 1], grid[1]), labels = FALSE, include.lowest = TRUE)
  yi <- cut(emb[, 2], .even_breaks(emb[, 2], grid[2]), labels = FALSE, include.lowest = TRUE)
  m <- tapply(v, list(factor(yi, seq_len(grid[2])), factor(xi, seq_len(grid[1]))), mean)
  m <- m[rev(seq_len(grid[2])), , drop = FALSE]
  dimnames(m) <- list(as.character(rev(seq_len(grid[2]))),
                      as.character(seq_len(grid[1])))
  m
}

# Facts for one feature. `ids` may be NULL.
.feature_facts <- function(emb, v, ids, grid, modal) {
  expressing <- v > 0
  mg <- .mean_grid(emb, v, grid)
  hot <- which(mg == max(mg, na.rm = TRUE), arr.ind = TRUE)[1, ]
  out <- list(
    n_expressing   = sum(expressing),
    pct_expressing = 100 * mean(expressing),
    max            = max(v),
    median_pos     = if (any(expressing)) median(v[expressing]) else NA_real_,
    mean_grid      = mg,
    hotspot        = list(row = rownames(mg)[hot[1]], col = colnames(mg)[hot[2]],
                          mean = mg[hot[1], hot[2]],
                          group = if (is.null(modal)) NULL else modal[hot[1], hot[2]])
  )
  if (!is.null(ids)) {
    ids <- droplevels(ids)
    by_group <- data.frame(
      group          = levels(ids),
      n              = as.integer(table(ids)),
      mean           = as.numeric(tapply(v, ids, mean)),
      pct_expressing = as.numeric(100 * tapply(expressing, ids, mean))
    )
    out$by_group <- by_group[order(by_group$mean, decreasing = TRUE), , drop = FALSE]
  }
  out
}

#' @exportS3Method BrailleR::Augment sc_featureplot
Augment.sc_featureplot <- function(x, grid = c(4L, 4L)) {
  emb  <- attr(x, "embedding")
  expr <- attr(x, "expr")
  ids  <- attr(x, "idents")
  modal <- if (is.null(ids)) NULL else .modal_grid(emb, ids, grid = grid)
  x$sc <- list(
    reduction  = attr(x, "reduction"),
    group_by   = attr(x, "group_by"),
    layer      = attr(x, "layer"),
    n_cells    = nrow(emb),
    axis_names = colnames(emb),
    grid_size  = grid,
    modal_grid = modal,
    features   = lapply(stats::setNames(colnames(expr), colnames(expr)),
                        function(f) .feature_facts(emb, expr[, f], ids, grid, modal))
  )
  BrailleR::Augment.ggplot(x)
}

#' VI() method for sc_featureplot
#'
#' For each feature: how many cells express it, the bin of the embedding
#' with the highest mean expression (and which group dominates that bin),
#' the top groups by mean expression, and a mean-expression grid aligned
#' with the [BR_DimPlot()] grids. `Describe = TRUE` adds the full
#' per-group table.
#'
#' @param x A `sc_featureplot`.
#' @param Describe If `TRUE`, also print the per-group table per feature.
#' @param ... Ignored.
#' @return Invisibly, the augmented plot.
#' @exportS3Method BrailleR::VI sc_featureplot
VI.sc_featureplot <- function(x, Describe = FALSE, ...) {
  x <- BrailleR::Augment(x)
  S <- x$sc
  cat(sprintf(
    "Feature plot of %d feature%s (%s) on the %s embedding of %d cells, values from %s.\n",
    length(S$features), if (length(S$features) == 1) "" else "s",
    paste(names(S$features), collapse = ", "),
    toupper(S$reduction), S$n_cells, S$layer
  ))
  grid_note <- sprintf(
    "%dx%d bins; rows = %s descending, cols = %s ascending; '.' = empty bin",
    S$grid_size[2], S$grid_size[1], S$axis_names[2], S$axis_names[1])

  for (f in names(S$features)) {
    F <- S$features[[f]]
    cat(sprintf("\n%s: expressed in %d cells (%.1f%%), maximum %.2f",
                f, F$n_expressing, F$pct_expressing, F$max))
    if (!is.na(F$median_pos)) cat(sprintf(", median %.2f among expressing cells", F$median_pos))
    cat(".\n")
    if (F$n_expressing == 0) next

    h <- F$hotspot
    cat(sprintf("Highest mean expression in bin row %s, column %s (mean %.2f)%s.\n",
                h$row, h$col, h$mean,
                if (is.null(h$group)) "" else sprintf(", where group %s is most common", h$group)))

    if (!is.null(F$by_group)) {
      top <- utils::head(F$by_group, 3)
      cat(sprintf("Top groups by mean: %s.\n", paste(sprintf(
        "%s (mean %.2f, %.0f%% expressing)", top$group, top$mean, top$pct_expressing),
        collapse = "; ")))
    }

    cat(sprintf("Mean-expression grid (%s):\n", grid_note))
    g <- formatC(F$mean_grid, format = "f", digits = 2)
    g[is.na(F$mean_grid)] <- "."
    print(noquote(g))

    if (isTRUE(Describe) && !is.null(F$by_group)) {
      cat("Per-group table (sorted by mean):\n")
      tab <- F$by_group
      tab$mean <- round(tab$mean, 2)
      tab$pct_expressing <- round(tab$pct_expressing, 1)
      print(tab, row.names = FALSE)
    }
  }
  invisible(x)
}
