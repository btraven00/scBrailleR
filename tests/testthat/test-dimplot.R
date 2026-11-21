### Tests for BR_DimPlot and its VI/Augment methods.

# Synthetic 2D data with three well-separated clusters used for
# Seurat-free unit tests of the grid helpers.
make_synth <- function(seed = 1) {
  set.seed(seed)
  emb <- rbind(
    cbind(rnorm(30, -3, 0.3), rnorm(30, -3, 0.3)),
    cbind(rnorm(30,  3, 0.3), rnorm(30,  3, 0.3)),
    cbind(rnorm(30,  3, 0.3), rnorm(30, -3, 0.3))
  )
  colnames(emb) <- c("D1", "D2")
  ids <- factor(rep(c("A", "B", "C"), each = 30))
  list(emb = emb, ids = ids)
}

# Helpers are unexported; reach in.
cluster_facts        <- scBrailleR:::.cluster_facts
modal_grid           <- scBrailleR:::.modal_grid
cluster_density_grid <- scBrailleR:::.cluster_density_grid


test_that(".cluster_facts: counts, percentages, centroids, isolation", {
  d <- make_synth()
  cf <- cluster_facts(d$emb, d$ids)

  expect_length(cf, 3)
  labs <- vapply(cf, `[[`, character(1), "label")
  expect_setequal(labs, c("A", "B", "C"))

  ns <- vapply(cf, `[[`, numeric(1), "n")
  expect_equal(sum(ns), 90)
  pcts <- vapply(cf, `[[`, numeric(1), "pct")
  expect_equal(sum(pcts), 100, tolerance = 1e-8)

  # Centroid of A should be near (-3, -3).
  cA <- cf[[which(labs == "A")]]$centroid
  expect_lt(abs(cA[1] - (-3)), 0.3)
  expect_lt(abs(cA[2] - (-3)), 0.3)

  # All three groups are well separated relative to their tiny spread.
  isol <- vapply(cf, `[[`, numeric(1), "isolation")
  expect_true(all(isol > 5))
})

test_that(".cluster_facts handles a single cluster (no isolation field)", {
  d <- make_synth()
  ids <- factor(rep("solo", nrow(d$emb)))
  cf <- cluster_facts(d$emb, ids)
  expect_length(cf, 1)
  expect_null(cf[[1]]$isolation)
  expect_equal(cf[[1]]$n, 90)
  expect_equal(cf[[1]]$pct, 100)
})

test_that(".modal_grid returns a character matrix of the right shape", {
  d <- make_synth()
  m <- modal_grid(d$emb, d$ids, grid = c(4, 4))
  expect_true(is.matrix(m))
  expect_equal(dim(m), c(4, 4))
  expect_true(is.character(m))
  # Every non-empty cell is one of A/B/C, empty cells are ".".
  expect_true(all(as.vector(m) %in% c("A", "B", "C", ".")))
  # With three well-separated clusters, at least 3 distinct labels appear.
  expect_gte(length(setdiff(unique(as.vector(m)), ".")), 3)
})

test_that(".modal_grid orientation: top row has highest y", {
  # Place a single point in the extreme upper-right corner only.
  emb <- rbind(c(0, 0), c(0.01, 0.01), c(1, 1))
  colnames(emb) <- c("D1", "D2")
  ids <- factor(c("low", "low", "high"))
  m <- modal_grid(emb, ids, grid = c(2, 2))
  expect_equal(m[1, 2], "high")   # top-right
  expect_equal(m[2, 1], "low")    # bottom-left
})

test_that(".cluster_density_grid counts add up to that cluster's n", {
  d <- make_synth()
  for (lab in levels(d$ids)) {
    g <- cluster_density_grid(d$emb, d$ids, lab, grid = c(4, 4))
    expect_equal(sum(g), sum(d$ids == lab))
    expect_equal(dim(g), c(4, 4))
  }
})


test_that("BR_DimPlot attaches embedding + idents and dispatches VI", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("SeuratObject")
  skip_if_not_installed("BrailleR")

  data("pbmc_small", package = "SeuratObject")
  p <- BR_DimPlot(pbmc_small, reduction = "tsne")

  expect_s3_class(p, "seurat_dimplot")
  expect_s3_class(p, "ggplot")

  emb <- attr(p, "embedding")
  expect_equal(ncol(emb), 2)
  expect_equal(nrow(emb), ncol(pbmc_small))   # one row per cell
  expect_s3_class(attr(p, "idents"), "factor")
  expect_length(attr(p, "idents"), ncol(pbmc_small))
  expect_identical(attr(p, "reduction"), "tsne")
})

test_that("Augment.seurat_dimplot populates Seurat facts and grids", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("SeuratObject")
  skip_if_not_installed("BrailleR")

  data("pbmc_small", package = "SeuratObject")
  a <- BrailleR::Augment(BR_DimPlot(pbmc_small, reduction = "tsne"))
  S <- a$Seurat

  expect_true(all(c("n_cells", "n_clusters", "per_cluster",
                    "density_grid", "modal_grid",
                    "cluster_density_grid") %in% names(S)))
  expect_equal(S$n_cells, ncol(pbmc_small))
  expect_gte(S$n_clusters, 1)

  # density_grid is a count matrix with margins: dims = grid + 1
  expect_equal(dim(S$density_grid), S$grid_size + 1L)
  expect_equal(dim(S$modal_grid),   rev(S$grid_size))

  # cluster_density_grid is a closure that returns a count matrix.
  first_lab <- S$per_cluster[[1]]$label
  g <- S$cluster_density_grid(first_lab)
  expect_equal(dim(g), rev(S$grid_size))
  expect_equal(sum(g), S$per_cluster[[1]]$n)
})

test_that("VI() prints per-cluster facts and both grids by default", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("SeuratObject")
  skip_if_not_installed("BrailleR")

  data("pbmc_small", package = "SeuratObject")
  out <- capture.output(BrailleR::VI(BR_DimPlot(pbmc_small, reduction = "tsne")))
  txt <- paste(out, collapse = "\n")

  expect_match(txt, "TSNE embedding of \\d+ cells")
  expect_match(txt, "Group .+: n=\\d+ \\(.+%\\)")
  expect_match(txt, "Most isolated group:")
  expect_match(txt, "Closest pair:")
  expect_match(txt, "Overall density grid")
  expect_match(txt, "Modal-cluster grid")
})

test_that("VI(Describe = TRUE) additionally prints a per-cluster density grid", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("SeuratObject")
  skip_if_not_installed("BrailleR")

  data("pbmc_small", package = "SeuratObject")
  out <- capture.output(
    BrailleR::VI(BR_DimPlot(pbmc_small, reduction = "tsne"), Describe = TRUE)
  )
  txt <- paste(out, collapse = "\n")
  expect_match(txt, "Per-cluster density grids:")
  # There must be at least one "Group X (n=...)" header in the per-cluster block.
  expect_match(txt, "Group .+ \\(n=\\d+\\):")
})
