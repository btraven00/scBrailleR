# scBrailleR

Accessible text descriptions of single-cell analysis plots, extending
[BrailleR](https://github.com/ajrgodfrey/BrailleR).

Status: experimental. Currently supports Seurat `ElbowPlot` and `DimPlot`.
`VlnPlot`, `DotPlot`, and `FeaturePlot` are planned.

## Install

scBrailleR is not on CRAN yet. Install it from a local source checkout.

### Requirements

- R >= 4.0
- [`BrailleR`](https://cran.r-project.org/package=BrailleR) (declared as
  `Imports`; installed automatically if you use `remotes::install_local`)
- For actually using the plot wrappers: `Seurat` and `SeuratObject`
  (declared as `Suggests` — install them yourself)

### One-liner from this checkout

```r
# install.packages("remotes")    # if you don't have it
remotes::install_local("~/phd/scBrailleR")
```

`install_local` will pick up `BrailleR` for you. Install the optional
single-cell stack separately:

```r
install.packages(c("Seurat", "SeuratObject"))
```

### From a git clone elsewhere

```sh
git clone https://github.com/btraven00/scBrailleR
```

```r
remotes::install_local("scBrailleR")
```

### Developing locally

If you want to iterate on the package without reinstalling each time:

```r
# from the package root
devtools::load_all(".")
devtools::test()          # run testthat suite
devtools::document()      # regenerate man/*.Rd from roxygen comments
```

### Verifying the install

```r
library(scBrailleR)
library(SeuratObject)
data("pbmc_small")
BrailleR::VI(BR_ElbowPlot(pbmc_small, ndims = 10))
```

You should see a short text summary ending with an "Estimated elbow at
PC..." line.

## Use

```r
library(scBrailleR)
library(SeuratObject)
data("pbmc_small")

# Elbow plot (Seurat::ElbowPlot wrapper)
p <- BR_ElbowPlot(pbmc_small, ndims = 15)
BrailleR::VI(p)                        # short summary
BrailleR::VI(p, Describe = TRUE)       # adds variance milestones + per-PC table

# DimPlot (UMAP / t-SNE / PCA scatter, coloured by cluster)
q <- BR_DimPlot(pbmc_small, reduction = "tsne")
BrailleR::VI(q)                        # per-cluster facts + density and modal grids
BrailleR::VI(q, Describe = TRUE)       # adds per-cluster density grids
```

The same functions accept a `SingleCellExperiment`. The DimPlot is drawn
with `scater::plotReducedDim()` and coloured by `colLabels()` unless you
set `group.by`. The ElbowPlot reads the `varExplained` attribute that
`scater::runPCA()` stores:

```r
library(scater)
sce <- runPCA(scuttle::logNormCounts(scuttle::mockSCE()))
colLabels(sce) <- scran::clusterCells(sce, use.dimred = "PCA")
BrailleR::VI(BR_ElbowPlot(sce))
BrailleR::VI(BR_DimPlot(sce))          # first reducedDim, coloured by label
```

Vignettes: `vignette("seurat", package = "scBrailleR")` and
`vignette("singlecellexperiment", package = "scBrailleR")`.

See `?BR_ElbowPlot` and `?BR_DimPlot` for full argument lists and runnable
examples.

## License

GPL-2, matching upstream BrailleR.

## Dev environment

`pixi run test` builds a conda env (R 4.4, Seurat, SingleCellExperiment,
scater), installs BrailleR from GitHub and `GenomeInfoDbData` from
Bioconductor, then runs the tests.
