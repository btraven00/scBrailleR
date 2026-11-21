### Tests for BR_FeaturePlot and its VI/Augment methods.

test_that(".mean_grid averages per bin with DimPlot orientation", {
  # Four cells in the four corners of a 2x2 grid.
  emb <- cbind(x = c(0, 1, 0, 1), y = c(0, 0, 1, 1))
  v   <- c(1, 2, 3, 4)
  g <- scBrailleR:::.mean_grid(emb, v, grid = c(2, 2))
  # top row = high y: (x=0,y=1)=3, (x=1,y=1)=4
  expect_equal(unname(g), rbind(c(3, 4), c(1, 2)))
  expect_equal(rownames(g), c("2", "1"))
})

test_that(".mean_grid leaves empty bins NA", {
  emb <- cbind(c(0, 1), c(0, 1))
  g <- scBrailleR:::.mean_grid(emb, c(5, 7), grid = c(2, 2))
  expect_equal(sum(is.na(g)), 2)
})

test_that(".feature_facts finds the hotspot and the top group", {
  emb <- cbind(c(0, 0, 1, 1), c(0, 0, 1, 1))
  v   <- c(0, 0, 3, 5)
  ids <- factor(c("a", "a", "b", "b"))
  modal <- scBrailleR:::.modal_grid(emb, ids, grid = c(2, 2))
  F <- scBrailleR:::.feature_facts(emb, v, ids, c(2, 2), modal)
  expect_equal(F$n_expressing, 2)
  expect_equal(F$pct_expressing, 50)
  expect_equal(F$hotspot$mean, 4)
  expect_identical(F$hotspot$group, "b")
  expect_identical(F$by_group$group[1], "b")
  expect_equal(F$by_group$pct_expressing, c(100, 0))
})

test_that("BR_FeaturePlot on Seurat describes each feature", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("BrailleR")
  data("pbmc_small", package = "SeuratObject")
  p <- BR_FeaturePlot(pbmc_small, c("CD3E", "LYZ"), reduction = "tsne")
  expect_s3_class(p, "sc_featureplot")
  expect_equal(dim(attr(p, "expr")), c(ncol(pbmc_small), 2))

  S <- BrailleR::Augment(p)$sc
  expect_named(S$features, c("CD3E", "LYZ"))
  expect_equal(nrow(S$features$CD3E$by_group),
               nlevels(SeuratObject::Idents(pbmc_small)))

  txt <- paste(capture.output(BrailleR::VI(p, Describe = TRUE)), collapse = "\n")
  expect_match(txt, "Feature plot of 2 features \\(CD3E, LYZ\\) on the TSNE embedding")
  expect_match(txt, "CD3E: expressed in \\d+ cells")
  expect_match(txt, "Highest mean expression in bin row \\d, column \\d")
  expect_match(txt, "Top groups by mean:")
  expect_match(txt, "Per-group table")
})

test_that("BR_FeaturePlot gives identical facts for Seurat and its SCE conversion", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("scater")
  skip_if_not_installed("BrailleR")
  data("pbmc_small", package = "SeuratObject")
  sce <- Seurat::as.SingleCellExperiment(pbmc_small)
  a <- BrailleR::Augment(BR_FeaturePlot(pbmc_small, "LYZ", reduction = "tsne"))$sc
  b <- BrailleR::Augment(BR_FeaturePlot(sce, "LYZ", reduction = "TSNE",
                                        group.by = "ident"))$sc
  expect_equal(b$features, a$features)
})

test_that("BR_FeaturePlot on SCE without labels skips group facts", {
  skip_if_not_installed("scater")
  skip_if_not_installed("BrailleR")
  set.seed(1)
  sce <- scater::runPCA(scuttle::logNormCounts(scuttle::mockSCE(ncells = 100, ngenes = 200)),
                        ncomponents = 5)
  p <- BR_FeaturePlot(sce, rownames(sce)[1])
  S <- BrailleR::Augment(p)$sc
  expect_null(S$features[[1]]$by_group)
  expect_null(S$modal_grid)
  txt <- paste(capture.output(BrailleR::VI(p)), collapse = "\n")
  expect_false(grepl("Top groups", txt))
  expect_error(BR_FeaturePlot(sce, "notAGene"), "features not found: notAGene")
})
