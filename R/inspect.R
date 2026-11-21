### scBrailleR: describe a dataset file without loading it
###
### picklerick::inspect() reads only metadata (shape, column names,
### embedding keys, sparse-matrix index pointers), so this is cheap even
### for millions of cells. We only turn its list into sentences.

#' Describe a single-cell dataset file without loading it
#'
#' Reads the file's metadata with [picklerick::inspect()] (no count
#' matrix is loaded, so this is fast and uses little memory on any size
#' of dataset) and describes it in words: format, number of cells and
#' genes, how sparse the matrix is, and which cell metadata, embeddings,
#' layers, graphs and unstructured keys it holds.
#'
#' Supported formats are those picklerick reads: AnnData `.h5ad` and
#' Zarr, MuData `.h5mu`, H5Seurat, 10x HDF5 and MatrixMarket, Visium and
#' BPCells.
#'
#' @param path Path to a dataset file or directory.
#' @return An object of class `"sc_dataset"` (the [picklerick::inspect()]
#'   list plus the path). Printing it, or calling [BrailleR::VI()] on it,
#'   gives the description; `VI(x, Describe = TRUE)` lists every column
#'   and key instead of the first 20.
#' @examples
#' \donttest{
#' if (requireNamespace("picklerick", quietly = TRUE)) {
#'   f <- system.file("extdata", "tiny.h5ad", package = "picklerick")
#'   BR_Inspect(f)
#' }
#' }
#' @export
BR_Inspect <- function(path) {
  if (!requireNamespace("picklerick", quietly = TRUE)) {
    stop("picklerick is required for BR_Inspect(); ",
         "install it with pak::pak(\"btraven00/scx/r/picklerick\")")
  }
  .as_sc_dataset(picklerick::inspect(path), path)
}

.as_sc_dataset <- function(info, path) {
  info$path <- path
  class(info) <- "sc_dataset"
  info
}

# "a, b, c" or "a, b, ... and 7 more".
.name_list <- function(x, max = 20L) {
  if (length(x) <= max) return(paste(x, collapse = ", "))
  sprintf("%s, and %d more", paste(x[seq_len(max)], collapse = ", "),
          length(x) - max)
}

.n <- function(n, one, many = paste0(one, "s")) {
  sprintf("%s %s", format(n, big.mark = ","), if (n == 1) one else many)
}

#' @exportS3Method BrailleR::VI sc_dataset
VI.sc_dataset <- function(x, Describe = FALSE, ...) {
  max <- if (isTRUE(Describe)) Inf else 20L
  cat(sprintf("%s file %s: %s by %s.\n", x$format, basename(x$path),
              .n(x$n_obs, "cell"), .n(x$n_vars, "gene")))

  s <- x$x_stats
  if (!is.null(s)) {
    # doubles: n_obs * n_vars overflows integer on large datasets.
    dens <- 100 * as.numeric(s$nnz) / (as.numeric(x$n_obs) * as.numeric(x$n_vars))
    cat(sprintf(paste0(
      "Main matrix: %s non-zero values (%.1f%% of entries). ",
      "Non-zeros per row (usually per cell): median %d, middle half %d to %d, maximum %d.\n"),
      format(s$nnz, big.mark = ","), dens,
      s$nnz_med, s$nnz_q1, s$nnz_q3, s$nnz_max))
  }

  absent <- character(0)
  section <- function(keys, label, what) {
    if (length(keys)) {
      cat(sprintf("%s, %s: %s.\n", label, .n(length(keys), what),
                  .name_list(keys, max)))
    } else {
      absent <<- c(absent, tolower(label))
    }
  }
  sparse_section <- function(m, label) {
    if (length(m$name)) {
      cat(sprintf("%s: %s.\n", label, paste(sprintf(
        "%s (%s non-zero values)", m$name, format(m$nnz, big.mark = ",")),
        collapse = "; ")))
    } else {
      absent <<- c(absent, tolower(label))
    }
  }
  section(x$obs_cols,  "Cell metadata", "column")
  section(x$var_cols,  "Gene metadata", "column")
  section(x$obsm_keys, "Cell embeddings (obsm)", "key")
  sparse_section(x$layers, "Extra layers")
  sparse_section(x$obsp, "Cell graphs (obsp)")
  section(x$varm_keys, "Gene loadings (varm)", "key")
  section(x$uns_keys,  "Unstructured metadata (uns)", "key")
  if (length(absent)) cat(sprintf("Not present: %s.\n", paste(absent, collapse = ", ")))
  invisible(x)
}

#' @export
print.sc_dataset <- function(x, ...) VI.sc_dataset(x, ...)
