### Tests for BR_ElbowPlot and its VI/Augment methods.

# .find_elbow is unexported; reach in for the unit tests.
find_elbow <- scBrailleR:::.find_elbow

test_that(".find_elbow returns NA for too-short input", {
  expect_true(is.na(find_elbow(numeric(0))))
  expect_true(is.na(find_elbow(c(5))))
  expect_true(is.na(find_elbow(c(5, 3))))
})

test_that(".find_elbow finds a sharp knee", {
  # Steep through index 3, then flat from index 4. The bottom of the
  # cliff (index 4) is the furthest point from the p1->pn chord, which
  # is what we want to report as "the elbow".
  y <- c(10, 9, 8, 2, 1.9, 1.8, 1.7, 1.6)
  expect_equal(find_elbow(y), 4)
})

test_that(".find_elbow is scale-invariant", {
  y <- c(10, 9, 8, 2, 1.9, 1.8, 1.7, 1.6)
  expect_equal(find_elbow(y), find_elbow(y * 1000))
  expect_equal(find_elbow(y), find_elbow(y * 1e-4))
})

test_that(".find_elbow on a perfectly linear decay has no strong elbow", {
  # All points lie on the chord from p1 to pn, so distances are ~0.
  # which.max still returns an index; we only assert it returns *something*.
  y <- seq(10, 1, length.out = 10)
  e <- find_elbow(y)
  expect_true(is.integer(e) || is.numeric(e))
  expect_true(e >= 1L && e <= length(y))
})


test_that("BR_ElbowPlot attaches sdev and reduction, and dispatches VI", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("SeuratObject")
  skip_if_not_installed("BrailleR")

  data("pbmc_small", package = "SeuratObject")
  p <- BR_ElbowPlot(pbmc_small, ndims = 10)

  expect_s3_class(p, "sc_elbowplot")
  expect_s3_class(p, "ggplot")
  expect_length(attr(p, "sdev"), 10)
  expect_identical(attr(p, "reduction"), "pca")
  expect_true(all(attr(p, "sdev") >= 0))
  expect_true(all(diff(attr(p, "sdev")) <= 1e-8))  # PC sdev is non-increasing
})

test_that("VI() output mentions PC count, sdev range, elbow PC, variance pct", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("SeuratObject")
  skip_if_not_installed("BrailleR")

  data("pbmc_small", package = "SeuratObject")
  p <- BR_ElbowPlot(pbmc_small, ndims = 10)
  out <- capture.output(BrailleR::VI(p))
  txt <- paste(out, collapse = "\n")

  expect_match(txt, "Elbow plot of 10 principal components")
  expect_match(txt, "reduction pca")
  expect_match(txt, "Standard deviations range from")
  expect_match(txt, "Estimated elbow at PC\\d+")
  expect_match(txt, "% of the variance")
})

test_that("Augment populates the Seurat fact list correctly", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("SeuratObject")
  skip_if_not_installed("BrailleR")

  data("pbmc_small", package = "SeuratObject")
  p <- BR_ElbowPlot(pbmc_small, ndims = 8)
  a <- BrailleR::Augment(p)

  expect_named(a$sc,
               c("reduction", "ndims", "sdev", "drops", "var_pct",
                 "cum_pct", "elbow", "sharpness", "milestones",
                 "table", "var_explained_cum"),
               ignore.order = TRUE)
  expect_equal(a$sc$ndims, 8)
  expect_length(a$sc$sdev, 8)
  expect_true(a$sc$elbow >= 1 && a$sc$elbow <= 8)
  expect_true(a$sc$var_explained_cum > 0 &&
              a$sc$var_explained_cum <= 1)
  # cumulative variance is monotone non-decreasing, ends at ~100%
  expect_true(all(diff(a$sc$cum_pct) >= -1e-8))
  expect_equal(unname(a$sc$cum_pct[8]), 100, tolerance = 1e-6)
  # table has expected columns
  expect_named(a$sc$table,
               c("PC", "sdev", "drop", "var_pct", "cum_pct"))
  expect_equal(nrow(a$sc$table), 8)
  expect_true(is.na(a$sc$table$drop[1]))
})

test_that("VI(Describe = TRUE) prints milestones and the per-PC table", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("SeuratObject")
  skip_if_not_installed("BrailleR")

  data("pbmc_small", package = "SeuratObject")
  p <- BR_ElbowPlot(pbmc_small, ndims = 15)
  out <- capture.output(BrailleR::VI(p, Describe = TRUE))
  txt <- paste(out, collapse = "\n")

  expect_match(txt, "Variance milestones:")
  expect_match(txt, "50% of variance")
  expect_match(txt, "75% of variance")
  expect_match(txt, "90% of variance")
  expect_match(txt, "Per-PC table")
  # header row of the data.frame
  expect_match(txt, "PC\\s+sdev\\s+drop\\s+var_pct\\s+cum_pct")
})

test_that("Sharpness label appears in the short summary when sdev is steep", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("SeuratObject")
  skip_if_not_installed("BrailleR")

  data("pbmc_small", package = "SeuratObject")
  p <- BR_ElbowPlot(pbmc_small, ndims = 15)
  out <- capture.output(BrailleR::VI(p))
  txt <- paste(out, collapse = "\n")
  expect_match(txt, "elbow is (sharp|moderate)\\.|decays gradually")
})

test_that("BR_ElbowPlot errors clearly when Seurat is missing", {
  # We can't actually uninstall Seurat here; instead verify that the guard
  # exists by reading the function body. This is a smoke check, not a deep test.
  src <- deparse(BR_ElbowPlot)
  expect_true(any(grepl("requireNamespace\\(\"Seurat\"", src)))
})


test_that("ndims smaller than 3 still returns a usable plot with NA elbow", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("SeuratObject")
  skip_if_not_installed("BrailleR")

  data("pbmc_small", package = "SeuratObject")
  p <- BR_ElbowPlot(pbmc_small, ndims = 2)
  a <- BrailleR::Augment(p)
  expect_true(is.na(a$sc$elbow))

  out <- capture.output(BrailleR::VI(p))
  expect_true(any(grepl("Too few PCs", out)))
})
