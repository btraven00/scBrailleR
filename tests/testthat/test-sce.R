### SingleCellExperiment support for BR_DimPlot / BR_ElbowPlot.

skip_sce <- function() {
  skip_if_not_installed("SingleCellExperiment")
  skip_if_not_installed("scater")
  skip_if_not_installed("scuttle")
  skip_if_not_installed("BrailleR")
}

mock_sce <- function() {
  set.seed(1)
  sce <- scuttle::logNormCounts(scuttle::mockSCE(ncells = 150, ngenes = 300))
  sce <- scater::runPCA(sce, ncomponents = 10)
  SingleCellExperiment::colLabels(sce) <- factor(sample(c("a", "b", "c"), 150, TRUE))
  sce
}

test_that("BR_DimPlot on SCE uses colLabels and the first reducedDim", {
  skip_sce()
  sce <- mock_sce()
  p <- BR_DimPlot(sce)
  expect_s3_class(p, "sc_dimplot")
  expect_identical(attr(p, "reduction"), "PCA")
  expect_identical(attr(p, "group_by"), "label")
  expect_equal(dim(attr(p, "embedding")), c(150, 2))
  expect_equal(levels(attr(p, "idents")), c("a", "b", "c"))

  txt <- paste(capture.output(BrailleR::VI(p)), collapse = "\n")
  expect_match(txt, "PCA embedding of 150 cells, coloured by label, across 3 groups")
  expect_match(txt, "Modal-cluster grid")
})

test_that("BR_DimPlot on SCE errors without labels or group.by", {
  skip_sce()
  sce <- mock_sce()
  sce$label <- NULL
  expect_error(BR_DimPlot(sce), "no column 'label'")
  expect_s3_class(BR_DimPlot(sce, group.by = "Treatment"), "sc_dimplot")
})

test_that("BR_DimPlot gives identical facts for Seurat and its SCE conversion", {
  skip_sce()
  skip_if_not_installed("Seurat")
  data("pbmc_small", package = "SeuratObject")
  sce <- Seurat::as.SingleCellExperiment(pbmc_small)
  a_seu <- BrailleR::Augment(BR_DimPlot(pbmc_small, reduction = "tsne"))$sc
  a_sce <- BrailleR::Augment(BR_DimPlot(sce, reduction = "TSNE", group.by = "ident"))$sc
  expect_equal(a_sce$per_cluster, a_seu$per_cluster)
  expect_equal(a_sce$modal_grid, a_seu$modal_grid)
})

test_that("BR_ElbowPlot on SCE reads varExplained from runPCA", {
  skip_sce()
  sce <- mock_sce()
  p <- BR_ElbowPlot(sce, ndims = 8)
  expect_s3_class(p, "sc_elbowplot")
  expect_s3_class(p, "ggplot")
  ve <- attr(SingleCellExperiment::reducedDim(sce, "PCA"), "varExplained")
  expect_equal(attr(p, "sdev"), sqrt(ve[1:8]))
  # ndims is capped at the number of PCs computed.
  expect_length(attr(BR_ElbowPlot(sce, ndims = 50), "sdev"), 10)

  txt <- paste(capture.output(BrailleR::VI(p)), collapse = "\n")
  expect_match(txt, "Elbow plot of 8 principal components from reduction PCA")
})

test_that("BR_ElbowPlot on SCE errors clearly without varExplained", {
  skip_sce()
  sce <- mock_sce()
  SingleCellExperiment::reducedDim(sce, "PCA") <- unclass(
    as.matrix(SingleCellExperiment::reducedDim(sce, "PCA")))[, 1:5]
  expect_error(BR_ElbowPlot(sce), "varExplained")
})

test_that("unsupported objects error clearly", {
  expect_error(BR_DimPlot(matrix(1:4, 2)), "supports Seurat and SingleCellExperiment")
  expect_error(BR_ElbowPlot(list()), "supports Seurat and SingleCellExperiment")
})
