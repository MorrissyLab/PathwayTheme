## h5ad adapter: aggregate single-cell .h5ad files into per-cluster pseudobulk.
##
## For every `.h5ad` and every cluster within it, take the mean (or sum)
## expression across cells -> one pseudobulk column `<sample>|cl<cluster>`.
##
## Reading is done with hdf5r so the package stays on CRAN-only dependencies.
## Both the modern AnnData layout (`/obs` as a group with `_index`, categorical
## columns as `categories`/`codes` sub-groups) and the legacy compound-dataset
## layout are supported, with dense, CSR and CSC `X`.

.pt_load_h5ad_adapter <- function(config) {
  if (!requireNamespace("hdf5r", quietly = TRUE))
    stop("h5ad input requires the 'hdf5r' package (install.packages(\"hdf5r\"))",
         call. = FALSE)

  paths <- .pt_resolve_h5ad_paths(config)
  if (!length(paths))
    stop("No .h5ad files found for the h5ad adapter", call. = FALSE)

  cluster_col <- config$cluster_col
  sample_col <- config$sample_col
  agg <- tolower(config$aggregate)
  if (!agg %in% c("mean", "sum"))
    stop("input$aggregate must be 'mean' or 'sum', got '", config$aggregate, "'",
         call. = FALSE)

  vecs <- list()
  meta_rows <- list()
  for (fp in paths) {
    f <- hdf5r::H5File$new(fp, mode = "r")
    on.exit(try(f$close_all(), silent = TRUE), add = TRUE)

    obs <- .pt_h5ad_read_frame(f, "obs")
    if (is.null(obs) || !cluster_col %in% names(obs)) {
      f$close_all()
      next
    }
    clusters <- as.character(obs[[cluster_col]])
    sample_id <- if (sample_col %in% names(obs)) {
      vals <- as.character(obs[[sample_col]])
      vals <- vals[!is.na(vals)]
      if (length(vals)) vals[[1L]] else tools::file_path_sans_ext(basename(fp))
    } else {
      tools::file_path_sans_ext(basename(fp))
    }

    genes <- .pt_h5ad_read_index(f, "var")
    n_obs <- length(clusters)

    # length-then-lexicographic cluster order, matching the Python adapter
    uniq <- unique(clusters)
    uniq <- uniq[order(nchar(uniq), uniq, method = "radix")]

    for (cl in uniq) {
      idx <- which(clusters == cl)
      if (!length(idx)) next
      vec <- .pt_h5ad_cluster_profile(f, idx, n_obs, length(genes), agg)
      names(vec) <- genes
      col <- paste0(sample_id, "|cl", cl)
      vecs[[length(vecs) + 1L]] <- list(name = col, values = vec)
      meta_rows[[length(meta_rows) + 1L]] <- data.frame(
        sample_id = as.character(sample_id),
        cluster_id = as.character(cl),
        sample_cluster_id = col,
        modality = "snRNA_cluster",
        source_file = basename(fp),
        n_cells = length(idx),
        stringsAsFactors = FALSE)
    }
    f$close_all()
  }

  if (!length(vecs))
    stop("no .h5ad file contained the cluster column '", cluster_col, "'",
         call. = FALSE)

  expr <- .pt_bind_gene_vectors(vecs)
  meta <- do.call(rbind, meta_rows)
  meta <- meta[!duplicated(meta$sample_cluster_id), , drop = FALSE]
  rownames(meta) <- meta$sample_cluster_id
  meta$sample_cluster_id <- NULL
  meta <- meta[colnames(expr), , drop = FALSE]
  pt_feature_matrix(expr, pt_sample_metadata(meta))
}

.pt_resolve_h5ad_paths <- function(config) {
  if (length(config$h5ad_paths)) return(as.character(config$h5ad_paths))
  if (.pt_is_set(config$h5ad_dir))
    return(sort(list.files(config$h5ad_dir, pattern = "\\.h5ad$",
                           full.names = TRUE), method = "radix"))
  character()
}

## Align a list of named gene vectors on the union of gene names, zero-fill the
## gaps, then average duplicate gene rows -- the R equivalent of
## `pd.concat(axis=1).fillna(0).groupby(level=0).mean()`.
.pt_bind_gene_vectors <- function(vecs) {
  all_genes <- sort(unique(unlist(lapply(vecs, function(v) names(v$values)),
                                  use.names = FALSE)), method = "radix")
  out <- matrix(0, nrow = length(all_genes), ncol = length(vecs),
                dimnames = list(all_genes, vapply(vecs, `[[`, "", "name")))
  for (j in seq_along(vecs)) {
    v <- vecs[[j]]$values
    # duplicate gene names within one file collapse to their mean, as in pandas
    if (anyDuplicated(names(v))) v <- tapply(v, names(v), mean)
    out[names(v), j] <- as.double(v)
  }
  out
}

## ---- HDF5 reading helpers -------------------------------------------------

.pt_h5_attr <- function(obj, name, default = NULL) {
  have <- tryCatch(hdf5r::h5attr_names(obj), error = function(e) character())
  if (!name %in% have) return(default)
  tryCatch(hdf5r::h5attr(obj, name), error = function(e) default)
}

.pt_h5ad_read_index <- function(f, group) {
  g <- f[[group]]
  if (inherits(g, "H5D")) {                      # legacy compound dataset
    df <- g$read()
    idx <- if ("index" %in% names(df)) df[["index"]] else df[[1L]]
    return(as.character(idx))
  }
  key <- .pt_h5_attr(g, "_index", "_index")
  if (!g$exists(key)) key <- if (g$exists("_index")) "_index" else names(g)[[1L]]
  as.character(g[[key]]$read())
}

## Read an AnnData dataframe group (or legacy compound dataset) into a
## data.frame of character/numeric columns.  Categoricals are decoded.
.pt_h5ad_read_frame <- function(f, group) {
  if (!f$exists(group)) return(NULL)
  g <- f[[group]]
  if (inherits(g, "H5D")) {
    df <- as.data.frame(g$read(), stringsAsFactors = FALSE)
    for (nm in names(df)) {
      cat_key <- paste0("uns/", nm, "_categories")
      if (f$exists(cat_key)) {
        cats <- as.character(f[[cat_key]]$read())
        df[[nm]] <- cats[as.integer(df[[nm]]) + 1L]
      }
    }
    return(df)
  }
  keys <- setdiff(names(g), c("__categories"))
  out <- list()
  for (nm in keys) {
    node <- g[[nm]]
    out[[nm]] <- if (inherits(node, "H5Group")) {
      # modern categorical: categories + codes
      if (node$exists("categories") && node$exists("codes")) {
        cats <- as.character(node[["categories"]]$read())
        codes <- as.integer(node[["codes"]]$read())
        ifelse(codes < 0L, NA_character_, cats[codes + 1L])
      } else next
    } else {
      v <- node$read()
      # legacy categorical: codes here, categories under obs/__categories/<nm>
      if (g$exists("__categories") && g[["__categories"]]$exists(nm)) {
        cats <- as.character(g[["__categories"]][[nm]]$read())
        v <- ifelse(as.integer(v) < 0L, NA_character_, cats[as.integer(v) + 1L])
      }
      v
    }
  }
  n <- max(vapply(out, length, 1L), 0L)
  out <- out[vapply(out, length, 1L) == n]
  as.data.frame(out, stringsAsFactors = FALSE, check.names = FALSE)
}

## Mean (or sum) expression profile over the cells in `idx`.
.pt_h5ad_cluster_profile <- function(f, idx, n_obs, n_var, agg) {
  x <- f[["X"]]
  if (inherits(x, "H5D")) {
    # hdf5r returns an h5py (n_obs, n_var) dataset as an R matrix of
    # dim (n_var, n_obs), so cells are already columns
    block <- x$read(args = list(seq_len(n_var), idx))
    block <- matrix(as.double(block), nrow = n_var)
    return(if (agg == "mean") rowMeans(block) else rowSums(block))
  }
  fmt <- .pt_h5_attr(x, "encoding-type", .pt_h5_attr(x, "h5sparse_format", "csr"))
  fmt <- tolower(sub("_matrix$", "", as.character(fmt)[[1L]]))
  data <- as.double(x[["data"]]$read())
  indices <- as.integer(x[["indices"]]$read())
  indptr <- as.integer(x[["indptr"]]$read())

  out <- numeric(n_var)
  if (fmt == "csr") {                     # rows = cells
    for (i in idx) {
      lo <- indptr[[i]] + 1L; hi <- indptr[[i + 1L]]
      if (hi < lo) next
      slice <- lo:hi
      cols <- indices[slice] + 1L
      out[cols] <- out[cols] + data[slice]
    }
  } else {                                 # csc: columns = cells
    keep <- rep(FALSE, n_obs)
    keep[idx] <- TRUE
    for (j in seq_len(n_var)) {
      lo <- indptr[[j]] + 1L; hi <- indptr[[j + 1L]]
      if (hi < lo) next
      slice <- lo:hi
      sel <- keep[indices[slice] + 1L]
      if (any(sel)) out[[j]] <- sum(data[slice][sel])
    }
  }
  if (agg == "mean") out / length(idx) else out
}
