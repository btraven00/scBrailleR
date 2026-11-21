### Tests for BR_Inspect (dataset description via picklerick::inspect).

sparse <- function(name = character(0), nnz = integer(0))
  list(name = name, nnz = nnz, nnz_q1 = nnz, nnz_med = nnz, nnz_q3 = nnz, nnz_max = nnz)

fake <- function(...) {
  info <- list(format = "H5AD", n_obs = 3000L, n_vars = 20000L,
               x_stats = list(nnz = 6000000L, nnz_q1 = 1500L, nnz_med = 2000L,
                              nnz_q3 = 2500L, nnz_max = 5000L),
               obs_cols = c("cell_type", "batch"), var_cols = character(0),
               obsm_keys = "X_umap", layers = sparse("counts", 6000000L),
               uns_keys = character(0), obsp = sparse(), varm_keys = character(0))
  scBrailleR:::.as_sc_dataset(utils::modifyList(info, list(...)), "/x/pbmc.h5ad")
}

test_that("VI.sc_dataset describes shape, sparsity and present keys", {
  txt <- paste(capture.output(print(fake())), collapse = "\n")
  expect_match(txt, "H5AD file pbmc.h5ad: 3,000 cells by 20,000 genes.", fixed = TRUE)
  expect_match(txt, "6,000,000 non-zero values (10.0% of entries)", fixed = TRUE)
  expect_match(txt, "median 2000, middle half 1500 to 2500, maximum 5000", fixed = TRUE)
  expect_match(txt, "Cell metadata, 2 columns: cell_type, batch.", fixed = TRUE)
  expect_match(txt, "Cell embeddings (obsm), 1 key: X_umap.", fixed = TRUE)
  expect_match(txt, "Extra layers: counts (6,000,000 non-zero values).", fixed = TRUE)
  expect_match(txt, "Not present: gene metadata, cell graphs (obsp), gene loadings (varm), unstructured metadata (uns).", fixed = TRUE)
})

test_that("long key lists are truncated unless Describe = TRUE", {
  x <- fake(obs_cols = paste0("c", 1:25))
  short <- paste(capture.output(BrailleR::VI(x)), collapse = "\n")
  long  <- paste(capture.output(BrailleR::VI(x, Describe = TRUE)), collapse = "\n")
  expect_match(short, "c20, and 5 more", fixed = TRUE)
  expect_match(long, "c24, c25.", fixed = TRUE)
})

test_that("density does not overflow on large datasets", {
  x <- fake(n_obs = 2000000L, n_vars = 30000L)
  expect_match(paste(capture.output(print(x)), collapse = "\n"), "0.0% of entries")
})

test_that("missing x_stats (MatrixMarket) skips the matrix line", {
  txt <- paste(capture.output(print(fake(x_stats = NULL))), collapse = "\n")
  expect_false(grepl("Main matrix", txt))
})

test_that("BR_Inspect reads real files through picklerick", {
  skip_if_not_installed("picklerick")
  d <- system.file("extdata", package = "picklerick")
  x <- BR_Inspect(file.path(d, "tiny.h5ad"))
  expect_s3_class(x, "sc_dataset")
  txt <- paste(capture.output(print(x)), collapse = "\n")
  expect_match(txt, "H5AD file tiny.h5ad: 6 cells by 5 genes.", fixed = TRUE)
  expect_match(txt, "Cell graphs (obsp): connectivities", fixed = TRUE)

  for (f in c("tiny_10x.h5", "tiny_mtx", "tiny_v4.h5seurat")) {
    expect_output(print(BR_Inspect(file.path(d, f))), "6 cells by 5 genes")
  }
})
