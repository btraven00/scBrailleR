### scBrailleR: accessible description of Seurat::ElbowPlot
###
### An ElbowPlot shows standard deviation of each principal component
### (y-axis) vs PC index (x-axis). The sighted-user task is "where does
### the curve flatten?" — we answer it numerically.

#' Accessible ElbowPlot for Seurat or SingleCellExperiment
#'
#' Wraps [Seurat::ElbowPlot()] (Seurat), or draws the equivalent ggplot
#' from the `varExplained` attribute that [scater::runPCA()] leaves on
#' the PCA reducedDim (SingleCellExperiment), and attaches the per-PC standard
#' deviations and a detected elbow index so that [BrailleR::VI()]
#' can describe the plot in words.
#'
#' @param object A Seurat object after `RunPCA()`, or a
#'   SingleCellExperiment after [scater::runPCA()].
#' @param ndims Number of PCs to plot.
#' @param reduction Name of the reduction (default `"pca"` for Seurat,
#'   `"PCA"` for SCE).
#' @param ... Passed to [Seurat::ElbowPlot()]; ignored for SCE.
#' @return A ggplot object with extra class `"sc_elbowplot"` and
#'   attributes used by `Augment`/`VI`.
#' @examples
#' \donttest{
#' if (requireNamespace("SeuratObject", quietly = TRUE) &&
#'     requireNamespace("Seurat", quietly = TRUE)) {
#'   data("pbmc_small", package = "SeuratObject")
#'
#'   # Default: short summary with a sharpness judgement.
#'   p <- BR_ElbowPlot(pbmc_small, ndims = 15)
#'   BrailleR::VI(p)
#'
#'   # Long form: variance milestones + a per-PC table.
#'   BrailleR::VI(p, Describe = TRUE)
#' }
#' }
#' @export
BR_ElbowPlot <- function(object, ndims = 20, reduction = NULL, ...) {
  if (inherits(object, "SingleCellExperiment")) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      stop("ggplot2 is required for BR_ElbowPlot() on a SingleCellExperiment")
    }
    reduction <- reduction %||% "PCA"
    pcs <- SingleCellExperiment::reducedDim(object, reduction)
    # scater::runPCA / scran::fixedPCA store per-PC variance here.
    ve <- attr(pcs, "varExplained")
    if (is.null(ve)) stop(sprintf(
      "reducedDim '%s' has no 'varExplained' attribute; run scater::runPCA()",
      reduction))
    ndims <- min(ndims, length(ve))
    sdev <- sqrt(ve[seq_len(ndims)])
    # ponytail: mimics Seurat::ElbowPlot; scater has no elbow plot.
    p <- ggplot2::ggplot(data.frame(dims = seq_len(ndims), stdev = sdev),
                         ggplot2::aes(dims, stdev)) +
      ggplot2::geom_point() +
      ggplot2::labs(x = reduction, y = "Standard Deviation")
  } else if (inherits(object, "Seurat")) {
    if (!requireNamespace("Seurat", quietly = TRUE)) {
      stop("Seurat is required for BR_ElbowPlot()")
    }
    reduction <- reduction %||% "pca"
    p <- Seurat::ElbowPlot(object, ndims = ndims, reduction = reduction, ...)
    sdev <- Seurat::Stdev(object, reduction = reduction)[seq_len(ndims)]
  } else {
    stop("BR_ElbowPlot() supports Seurat and SingleCellExperiment objects, not ",
         class(object)[1])
  }
  attr(p, "sdev")      <- sdev
  attr(p, "reduction") <- reduction
  class(p) <- c("sc_elbowplot", class(p))
  p
}

# Detect the elbow as the point of maximum perpendicular distance from
# the line joining the first and last (PC, sdev) points. Simple and
# scale-invariant; matches the heuristic blind users would ask us to apply.
.find_elbow <- function(y) {
  n <- length(y)
  if (n < 3) return(NA_integer_)
  x <- seq_along(y)
  p1 <- c(x[1], y[1]); p2 <- c(x[n], y[n])
  v  <- p2 - p1
  vn <- sqrt(sum(v * v))
  dists <- vapply(seq_len(n), function(i) {
    w <- c(x[i], y[i]) - p1
    abs(w[1] * v[2] - w[2] * v[1]) / vn
  }, numeric(1))
  which.max(dists)
}

#' @exportS3Method BrailleR::Augment sc_elbowplot
Augment.sc_elbowplot <- function(x) {
  sdev  <- attr(x, "sdev")
  elbow <- .find_elbow(sdev)
  var <- sdev^2
  total_var <- sum(var)
  drops <- c(NA_real_, -diff(sdev))   # positive when sdev decreases
  table <- data.frame(
    PC        = seq_along(sdev),
    sdev      = round(sdev, 3),
    drop      = round(drops, 3),
    var_pct   = round(100 * var / total_var, 2),
    cum_pct   = round(100 * cumsum(var) / total_var, 2)
  )
  milestones <- vapply(c(0.5, 0.75, 0.9), function(t) {
    which(cumsum(var) / total_var >= t)[1]
  }, integer(1))
  names(milestones) <- c("p50", "p75", "p90")

  # Sharpness: drop at the elbow vs median drop elsewhere.
  sharpness <- if (!is.na(elbow) && length(sdev) >= 4) {
    other <- drops[-c(1, elbow)]
    other <- other[!is.na(other)]
    if (length(other) && median(other) > 0) drops[elbow] / median(other) else NA_real_
  } else NA_real_

  x$sc <- list(
    reduction = attr(x, "reduction"),
    ndims     = length(sdev),
    sdev      = sdev,
    drops     = drops,
    var_pct   = 100 * var / total_var,
    cum_pct   = 100 * cumsum(var) / total_var,
    elbow     = elbow,
    sharpness = sharpness,
    milestones = milestones,
    table     = table,
    var_explained_cum = if (!is.na(elbow))
      sum(var[seq_len(elbow)]) / total_var else NA_real_
  )
  BrailleR::Augment.ggplot(x)
}

#' @exportS3Method BrailleR::VI sc_elbowplot
VI.sc_elbowplot <- function(x, Describe = FALSE, ...) {
  x <- BrailleR::Augment(x)
  S <- x$sc

  cat(sprintf(
    "Elbow plot of %d principal components from reduction %s.\n",
    S$ndims, S$reduction
  ))
  cat(sprintf(
    "Standard deviations range from %.2f (PC1) to %.2f (PC%d).\n",
    S$sdev[1], S$sdev[S$ndims], S$ndims
  ))
  if (!is.na(S$elbow)) {
    sharp_lab <- if (is.na(S$sharpness)) {
      ""
    } else if (S$sharpness >= 2) {
      " The elbow is sharp."
    } else if (S$sharpness >= 1.2) {
      " The elbow is moderate."
    } else {
      " The curve decays gradually; the elbow position is approximate."
    }
    cat(sprintf(
      "Estimated elbow at PC%d (sdev %.2f). PCs 1-%d together account for %.1f%% of the variance among the %d PCs shown.%s\n",
      S$elbow, S$sdev[S$elbow], S$elbow,
      100 * S$var_explained_cum, S$ndims, sharp_lab
    ))
  } else {
    cat("Too few PCs to estimate an elbow.\n")
  }

  if (isTRUE(Describe)) {
    cat("\nVariance milestones:\n")
    for (nm in names(S$milestones)) {
      pc <- S$milestones[[nm]]
      pct <- c(p50 = 50, p75 = 75, p90 = 90)[[nm]]
      if (is.na(pc)) {
        cat(sprintf("  %d%% of variance: not reached within the %d PCs shown.\n",
                    pct, S$ndims))
      } else {
        cat(sprintf("  %d%% of variance reached at PC%d (cumulative %.1f%%).\n",
                    pct, pc, S$cum_pct[pc]))
      }
    }
    cat("\nPer-PC table (sdev, drop from previous PC, variance %, cumulative %):\n")
    print(S$table, row.names = FALSE)
  }

  invisible(x)
}
