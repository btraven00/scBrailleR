### scBrailleR: accessible description of Seurat::DimPlot
###
### A DimPlot is a 2D scatter of a reduction (UMAP, t-SNE, PCA), with
### each cell coloured by an identity / cluster label. A sighted user
### reads three things from it:
###   1. how many clusters, how big each is;
###   2. where each cluster sits and how spread out it is;
###   3. which clusters are isolated and which overlap.
### A second layer of information is the *shape* of the embedding: the
### overall density landscape, and which region of the plot is "owned"
### by which cluster. We answer both layers numerically, the second one
### using a BrailleR::WhereXY-style grid summary.

#' Accessible Seurat DimPlot
#'
#' Wraps [Seurat::DimPlot()] and attaches a minimal facts list (per-
#' cluster size, centroid, spread, plus the raw embedding and idents)
#' so that [BrailleR::VI()] can describe the embedding in words.
#'
#' Note: the full Seurat object is *not* attached to the plot — only
#' the two embedding columns and the identity factor are kept.
#'
#' @param object A Seurat object with at least one reduction computed.
#' @param reduction Reduction name (default: Seurat's default reduction).
#' @param group.by Metadata column to colour by; defaults to active idents.
#' @param dims Which two dimensions of the reduction to plot (default
#'   `c(1, 2)`).
#' @param ... Passed to [Seurat::DimPlot()].
#' @return A `ggplot` with extra class `"seurat_dimplot"` and the
#'   attributes consumed by `Augment.seurat_dimplot()`.
#' @examples
#' \donttest{
#' if (requireNamespace("SeuratObject", quietly = TRUE) &&
#'     requireNamespace("Seurat", quietly = TRUE)) {
#'   data("pbmc_small", package = "SeuratObject")
#'
#'   # Default output: per-cluster facts + overall and modal-cluster grids.
#'   p <- BR_DimPlot(pbmc_small, reduction = "tsne")
#'   BrailleR::VI(p)
#'
#'   # Long form: additionally prints a per-cluster density grid.
#'   BrailleR::VI(p, Describe = TRUE)
#'
#'   # Colour by a different metadata column instead of active idents.
#'   BrailleR::VI(BR_DimPlot(pbmc_small, reduction = "tsne",
#'                           group.by = "groups"))
#' }
#' }
#' @export
BR_DimPlot <- function(object,
                       reduction = NULL,
                       group.by  = NULL,
                       dims      = c(1, 2),
                       ...) {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat is required for BR_DimPlot()")
  }
  reduction <- if (is.null(reduction)) Seurat::DefaultDimReduc(object) else reduction
  p <- Seurat::DimPlot(object,
                       reduction = reduction,
                       group.by  = group.by,
                       dims      = dims,
                       ...)

  emb <- Seurat::Embeddings(object, reduction = reduction)[, dims, drop = FALSE]
  ids <- if (is.null(group.by)) Seurat::Idents(object) else object[[group.by, drop = TRUE]]

  attr(p, "embedding") <- emb
  attr(p, "idents")    <- as.factor(ids)
  attr(p, "reduction") <- reduction
  attr(p, "dims")      <- dims
  attr(p, "group_by")  <- group.by
  class(p) <- c("seurat_dimplot", class(p))
  p
}

# Per-cluster centroid + spread + isolation.
# isolation = (distance to nearest other centroid) / (geom. mean of own MADs).
# Higher = more isolated. NA when only one cluster.
.cluster_facts <- function(emb, idents) {
  levs <- levels(droplevels(idents))
  if (!length(levs)) return(list())

  rows <- lapply(levs, function(lev) {
    keep <- idents == lev
    if (!any(keep)) return(NULL)
    e <- emb[keep, , drop = FALSE]
    list(
      label    = lev,
      n        = sum(keep),
      centroid = c(median(e[, 1]), median(e[, 2])),
      spread   = c(mad(e[, 1]),    mad(e[, 2]))
    )
  })
  rows <- Filter(Negate(is.null), rows)
  total <- sum(vapply(rows, `[[`, numeric(1), "n"))
  for (i in seq_along(rows)) rows[[i]]$pct <- 100 * rows[[i]]$n / total

  if (length(rows) > 1) {
    cents <- do.call(rbind, lapply(rows, `[[`, "centroid"))
    d <- as.matrix(dist(cents))
    diag(d) <- Inf
    for (i in seq_along(rows)) {
      own_spread <- sqrt(prod(pmax(rows[[i]]$spread, 1e-9)))
      rows[[i]]$nearest_d   <- min(d[i, ])
      rows[[i]]$isolation   <- rows[[i]]$nearest_d / own_spread
      rows[[i]]$nearest_lab <- rows[[which.min(d[i, ])]]$label
    }
  }
  rows
}

# Equal-width bin breaks over [min, max], with a tiny pad so the
# extreme points fall inside the outermost bin (matches WhereXY's cut()).
.even_breaks <- function(v, n) {
  rng <- range(v)
  pad <- 1e-6 * diff(rng) + .Machine$double.eps
  seq(rng[1] - pad, rng[2] + pad, length.out = n + 1)
}

# Modal-cluster grid: for each (x-bin, y-bin) cell, the cluster label
# with the most cells in that bin. Empty cells are "." and ties are
# resolved by taking the first level. Rows are flipped so that higher
# y appears at the top of the printed matrix (matches WhereXY).
.modal_grid <- function(emb, idents, grid = c(4, 4)) {
  xb <- .even_breaks(emb[, 1], grid[1])
  yb <- .even_breaks(emb[, 2], grid[2])
  xi <- cut(emb[, 1], breaks = xb, labels = FALSE, include.lowest = TRUE)
  yi <- cut(emb[, 2], breaks = yb, labels = FALSE, include.lowest = TRUE)
  out <- matrix(".", nrow = grid[2], ncol = grid[1])
  for (cx in seq_len(grid[1])) {
    for (cy in seq_len(grid[2])) {
      sel <- xi == cx & yi == cy
      if (!any(sel, na.rm = TRUE)) next
      tab <- sort(table(droplevels(idents[sel])), decreasing = TRUE)
      out[cy, cx] <- names(tab)[1]
    }
  }
  out <- out[rev(seq_len(grid[2])), , drop = FALSE]
  dimnames(out) <- list(as.character(rev(seq_len(grid[2]))),
                        as.character(seq_len(grid[1])))
  out
}

# Density grid restricted to one cluster, using global breaks so it
# aligns with the overall and modal grids.
.cluster_density_grid <- function(emb, idents, label, grid = c(4, 4)) {
  sel <- idents == label
  e   <- emb[sel, , drop = FALSE]
  xb  <- .even_breaks(emb[, 1], grid[1])
  yb  <- .even_breaks(emb[, 2], grid[2])
  xi  <- cut(e[, 1], breaks = xb, labels = FALSE, include.lowest = TRUE)
  yi  <- cut(e[, 2], breaks = yb, labels = FALSE, include.lowest = TRUE)
  xi  <- factor(xi, levels = seq_len(grid[1]))
  yi  <- factor(yi, levels = seq_len(grid[2]))
  out <- table(yi, xi)
  out <- out[rev(seq_len(grid[2])), , drop = FALSE]
  storage.mode(out) <- "integer"
  dimnames(out) <- list(as.character(rev(seq_len(grid[2]))),
                        as.character(seq_len(grid[1])))
  out
}

#' Augment a seurat_dimplot with the facts needed for description
#'
#' Populates `x$Seurat` with per-cluster summaries and three grids:
#' an overall density grid (via [BrailleR::WhereXY]), a modal-cluster
#' grid, and per-cluster density grids (only built for the `Describe`
#' path of `VI()`).
#'
#' @param x A `seurat_dimplot`.
#' @param grid Two-integer vector `c(nx, ny)` giving grid resolution.
#' @return The input with class `Augmented` appended and `x$Seurat`
#'   populated.
#' @exportS3Method BrailleR::Augment seurat_dimplot
Augment.seurat_dimplot <- function(x, grid = c(4L, 4L)) {
  emb <- attr(x, "embedding")
  ids <- attr(x, "idents")
  per <- .cluster_facts(emb, ids)

  density_grid <- BrailleR::WhereXY(emb[, 1], emb[, 2],
                                    grid = grid, addmargins = TRUE)
  modal_grid <- .modal_grid(emb, ids, grid = grid)

  x$Seurat <- list(
    reduction    = attr(x, "reduction"),
    dims         = attr(x, "dims"),
    group_by     = attr(x, "group_by"),
    n_cells      = nrow(emb),
    n_clusters   = length(per),
    axis_names   = colnames(emb),
    axis_ranges  = list(x = range(emb[, 1]), y = range(emb[, 2])),
    per_cluster  = per,
    grid_size    = grid,
    density_grid = density_grid,
    modal_grid   = modal_grid,
    # closures so we can build per-cluster grids on demand without
    # bloating the augmented object up front.
    cluster_density_grid = function(label)
      .cluster_density_grid(emb, ids, label, grid = grid)
  )
  BrailleR::Augment.ggplot(x)
}

#' VI() method for seurat_dimplot
#'
#' Default output reports overall counts, per-cluster size/centroid/
#' spread, pairwise isolation/closest-pair, the overall density grid,
#' and the modal-cluster grid. `Describe = TRUE` additionally prints a
#' per-cluster density grid for each cluster.
#'
#' @param x A `seurat_dimplot`.
#' @param Describe If `TRUE`, also print per-cluster density grids.
#' @param ... Ignored.
#' @return Invisibly, the augmented plot.
#' @exportS3Method BrailleR::VI seurat_dimplot
VI.seurat_dimplot <- function(x, Describe = FALSE, ...) {
  x <- BrailleR::Augment(x)
  S <- x$Seurat

  group_lab <- if (is.null(S$group_by)) "active identity" else S$group_by
  cat(sprintf(
    "%s embedding of %d cells, coloured by %s, across %d groups.\n",
    toupper(S$reduction), S$n_cells, group_lab, S$n_clusters
  ))
  cat(sprintf(
    "Axes: %s in [%.2f, %.2f]; %s in [%.2f, %.2f].\n",
    S$axis_names[1], S$axis_ranges$x[1], S$axis_ranges$x[2],
    S$axis_names[2], S$axis_ranges$y[1], S$axis_ranges$y[2]
  ))

  per <- S$per_cluster
  ord <- order(vapply(per, `[[`, numeric(1), "n"), decreasing = TRUE)
  for (i in ord) {
    cl <- per[[i]]
    cat(sprintf(
      "Group %s: n=%d (%.1f%%), centred near (%.2f, %.2f), spread (%.2f, %.2f).\n",
      cl$label, cl$n, cl$pct,
      cl$centroid[1], cl$centroid[2],
      cl$spread[1],   cl$spread[2]
    ))
  }

  if (S$n_clusters > 1) {
    isol <- vapply(per, function(cl) cl$isolation %||% NA_real_, numeric(1))
    most_iso  <- per[[which.max(isol)]]
    least_iso <- per[[which.min(isol)]]
    cat(sprintf(
      "Most isolated group: %s (nearest centroid is %s at distance %.2f, %.1fx its own spread).\n",
      most_iso$label, most_iso$nearest_lab,
      most_iso$nearest_d, most_iso$isolation
    ))
    cat(sprintf(
      "Closest pair: %s and %s (centroid distance %.2f, %.1fx %s's spread).\n",
      least_iso$label, least_iso$nearest_lab,
      least_iso$nearest_d, least_iso$isolation, least_iso$label
    ))
  }

  cat(sprintf(
    "\nOverall density grid (%dx%d bins of equal width; rows = %s descending, cols = %s ascending; row/col sums in margins):\n",
    S$grid_size[2], S$grid_size[1], S$axis_names[2], S$axis_names[1]
  ))
  print(S$density_grid)

  cat("\nModal-cluster grid (label of the most-populous group in each bin; '.' = empty bin):\n")
  print(noquote(S$modal_grid))

  if (isTRUE(Describe)) {
    cat("\nPer-cluster density grids:\n")
    for (i in ord) {
      cl <- per[[i]]
      cat(sprintf("\nGroup %s (n=%d):\n", cl$label, cl$n))
      print(S$cluster_density_grid(cl$label))
    }
  }

  invisible(x)
}

`%||%` <- function(a, b) if (is.null(a)) b else a
