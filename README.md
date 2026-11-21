# scBrailleR

Accessible text descriptions of single-cell analysis plots, extending
[BrailleR](https://github.com/ajrgodfrey/BrailleR).

Status: experimental. Describes `ElbowPlot`, `DimPlot` and `FeaturePlot`
for Seurat and SingleCellExperiment objects, and dataset files via
picklerick. `VlnPlot` and `DotPlot` are planned.

## Install

```r
# install.packages("pak")
pak::pak(c("ajrgodfrey/BrailleR", "btraven00/scBrailleR"))
```

BrailleR is installed from GitHub because CRAN does not always have a
build for current R. The object-specific parts are optional; install
what you use:

```r
pak::pak(c("Seurat", "SeuratObject"))                     # Seurat plots
pak::pak(c("SingleCellExperiment", "scater", "scuttle"))  # SCE plots
pak::pak("btraven00/scx/r/picklerick")                    # BR_Inspect(); needs Rust
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

### Describing a dataset file without loading it

`BR_Inspect()` uses [picklerick](https://github.com/btraven00/scx)'s
`inspect()`, which reads only metadata, so it is fast on any size of
`.h5ad`, Zarr, `.h5mu`, H5Seurat, 10x or BPCells file:

```r
BR_Inspect("pbmc3k.h5ad")
#> H5AD file pbmc3k.h5ad: 2,700 cells by 32,738 genes.
#> Main matrix: ... non-zero values (...% of entries). Non-zeros per row ...
#> Cell metadata, 4 columns: ...
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
