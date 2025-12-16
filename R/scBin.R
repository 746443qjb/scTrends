#' Bin single-cell pseudotime expression profiles
#'
#' @description
#' `scBin()` bins single-cell expression data along a pseudotime axis and
#' computes bin-averaged gene expression profiles. It supports both
#' non-overlapping and overlapping (sliding window) binning, and can
#' optionally filter lowly expressed genes across bins.
#'
#' @param seurat_obj A Seurat object containing the single-cell data.
#' @param pseudotime_col Character scalar. Name of the column in
#'   `seurat_obj@meta.data` that stores pseudotime values.
#' @param assay Character scalar. Assay to use when extracting the
#'   expression matrix. Default is `"RNA"`.
#' @param slot Character scalar. Slot to use when extracting the
#'   expression matrix via `Seurat::GetAssayData()`. Common choices are
#'   `"data"` (normalized) or `"counts"`. Default is `"data"`.
#' @param n_bins Integer. Number of bins along pseudotime. Default is 40.
#' @param overlap_ratio Numeric in [0, 1). If 0, bins are non-overlapping.
#'   Values between 0 and 1 specify the overlap ratio for a sliding window
#'   binning scheme (higher values = more overlap). Default is 0.
#' @param pseudotime_order Either `"increasing"` or `"decreasing"`.
#'   Controls whether cells are ordered from low to high pseudotime or the
#'   reverse. Default is `"increasing"`.
#' @param x_range Numeric vector of length 2 giving the x-axis range (e.g. 0 to 1).
#' @param filter_genes Logical. If `TRUE`, filter genes based on their
#'   expression across bins (see `max_zero_bins_pct` and `min_expression`).
#'   Default is `TRUE`.
#' @param max_zero_bins_pct Numeric between 0 and 1.
#'   Maximum allowed proportion of bins with expression below `min_expression`.
#'   Default is 0.5 (50%).
#' @param min_expression Numeric. Minimum expression threshold used when
#'   computing the proportion of zero/low bins per gene. Default is 0.
#' @param min_cells_per_bin Integer. Minimum number of cells allowed in a
#'   bin. Bins with fewer cells trigger a warning. Default is 10.
#' @param verbose Logical. If `TRUE`, print progress messages. Default is
#'   `TRUE`.
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Extracts pseudotime values from `seurat_obj@meta.data`.
#'   \item Orders cells by pseudotime (increasing or decreasing).
#'   \item Extracts the expression matrix from the specified assay and slot.
#'   \item Bins cells into `n_bins` using either non-overlapping quantile-based
#'         bins or a sliding window with `overlap_ratio`.
#'   \item Computes bin-averaged expression for each gene.
#'   \item Optionally filters genes based on their expression across bins.
#'   \item Creates a pseudotime coordinate vector `x` over `x_range`.
#'   \item Computes summary statistics of pseudotime within each bin.
#' }
#'
#' @return
#' A list of class `"scBin"` with components:
#' \itemize{
#'   \item `bin_means`: A matrix of bin-averaged expression
#'         (genes × bins).
#'   \item `x`: Numeric vector of pseudotime coordinates (length = `n_bins`).
#'   \item `bin_info`: A data frame with per-bin statistics:
#'         `bin_id`, `n_cells`, `pseudotime_median`, `pseudotime_mean`,
#'         `pseudotime_min`, `pseudotime_max`, `x_coord`.
#'   \item `metadata`: A list of parameters and summary values used in the
#'         binning procedure.
#'   \item `cell_to_bin`: A data frame mapping each cell to a bin:
#'         `cell`, `pseudotime`, `bin_id`.
#'   \item `bin_to_cells`: A list where each element is an integer vector of
#'         cell indices belonging to the corresponding bin.
#'   \item `gene_filter_stats`: A data frame with per-gene filtering
#'         information (`gene`, `n_zero_bins`, `zero_bins_pct`, `kept`),
#'         or `NULL` if `filter_genes = FALSE`.
#' }
#'
#' @importFrom Matrix rowMeans
#' @importFrom Seurat DefaultAssay GetAssayData
#' @export
#'
#' @examples
#' \dontrun{
#'   res <- scBin(
#'     seurat_obj = seurat_obj,
#'     pseudotime_col = "pseudotime",
#'     n_bins = 40,
#'     overlap_ratio = 0,
#'     pseudotime_order = "increasing"
#'   )
#'
#'   str(res)
#' }
scBin <- function(seurat_obj,
                  pseudotime_col,
                  assay = "RNA",
                  slot = "data",
                  n_bins = 40,
                  overlap_ratio = 0,
                  pseudotime_order = c("increasing", "decreasing"),
                  x_range = c(0, 1),
                  filter_genes = TRUE,
                  max_zero_bins_pct = 0.5,
                  min_expression = 0,
                  min_cells_per_bin = 10,
                  verbose = TRUE) {

  # ============================================================
  # Argument checks and normalization
  # ============================================================
  pseudotime_order <- match.arg(pseudotime_order)

  if (!is.numeric(overlap_ratio) || overlap_ratio < 0 || overlap_ratio >= 1) {
    stop("`overlap_ratio` must be a numeric value between 0 (no overlap) and < 1.")
  }

  if (verbose) {
    cat("╔══════════════════════════════════════════════════════════╗\n")
    cat("║        scBin: Binning Single-Cell Pseudotime Data       ║\n")
    cat("╚══════════════════════════════════════════════════════════╝\n\n")
  }

  # ============================================================
  # Step 1: Extract pseudotime
  # ============================================================
  if (verbose) cat("Step 1: Extracting pseudotime values...\n")

  if (!pseudotime_col %in% colnames(seurat_obj@meta.data)) {
    stop(
      "Pseudotime column '", pseudotime_col, "' not found in metadata.\n",
      "Available columns: ", paste(colnames(seurat_obj@meta.data), collapse = ", ")
    )
  }

  pt <- seurat_obj@meta.data[[pseudotime_col]]

  if (!is.numeric(pt)) {
    stop(
      "Pseudotime column '", pseudotime_col, "' must be numeric. ",
      "Current type: ", paste(class(pt), collapse = ", ")
    )
  }

  if (any(is.na(pt))) {
    n_na <- sum(is.na(pt))
    warning(
      "Found ", n_na, " cells with NA pseudotime values. ",
      "These cells will be removed."
    )
    valid_cells <- !is.na(pt)
    seurat_obj <- seurat_obj[, valid_cells]
    pt <- pt[valid_cells]
  }

  names(pt) <- colnames(seurat_obj)

  if (verbose) {
    cat("  ✓ Extracted pseudotime for", length(pt), "cells\n")
    cat("  Range: [", round(min(pt), 4), ", ", round(max(pt), 4), "]\n", sep = "")
  }

  # ============================================================
  # Step 2: Order cells by pseudotime
  # ============================================================
  if (verbose) cat("\nStep 2: Ordering cells by pseudotime...\n")

  if (pseudotime_order == "increasing") {
    ord <- order(pt, decreasing = FALSE)
    if (verbose) cat("  ✓ Order: increasing (early → late)\n")
  } else {
    ord <- order(pt, decreasing = TRUE)
    if (verbose) cat("  ✓ Order: decreasing (early → late)\n")
  }

  pt_ord  <- pt[ord]
  obj_ord <- seurat_obj[, ord]

  # ============================================================
  # Step 3: Extract expression matrix
  # ============================================================
  if (verbose) cat("\nStep 3: Extracting expression matrix...\n")

  Seurat::DefaultAssay(obj_ord) <- assay
  expr_mat <- Seurat::GetAssayData(obj_ord, slot = slot, assay = assay)

  if (verbose) {
    cat("  ✓ Matrix dimensions:", nrow(expr_mat), "genes ×", ncol(expr_mat), "cells\n")
    cat("  Assay:", assay, "| Slot:", slot, "\n")
  }

  # ============================================================
  # Step 4: Bin cells (overlapping or non-overlapping)
  # ============================================================
  if (verbose) {
    cat("\nStep 4: Binning cells into", n_bins, "bins")
    if (overlap_ratio > 0) {
      cat(" (overlap ratio:", round(overlap_ratio * 100), "%)\n", sep = "")
    } else {
      cat(" (no overlap)\n")
    }
  }

  n_cells <- length(pt_ord)

  if (overlap_ratio == 0) {
    # ---- Non-overlapping bins (by rank) ----
    cell_rank <- rank(pt_ord, ties.method = "first")
    bin_id <- cut(cell_rank, breaks = n_bins, labels = FALSE)
    bin_id <- as.integer(bin_id)

    bin_to_cells <- lapply(seq_len(n_bins), function(b) which(bin_id == b))

  } else {
    # ---- Overlapping bins: sliding window ----
    window_size <- ceiling(n_cells / n_bins)
    step_size   <- ceiling(window_size * (1 - overlap_ratio))
    if (step_size < 1) step_size <- 1

    bin_to_cells <- list()
    for (i in seq_len(n_bins)) {
      start_idx <- (i - 1) * step_size + 1
      end_idx   <- min(start_idx + window_size - 1, n_cells)

      # Ensure the last bin captures all remaining cells
      if (i == n_bins) {
        end_idx <- n_cells
      }

      bin_to_cells[[i]] <- start_idx:end_idx
    }

    # For compatibility: assign each cell to its first bin
    bin_id <- rep(NA_integer_, n_cells)
    for (b in seq_len(n_bins)) {
      cells_in_bin <- bin_to_cells[[b]]
      unassigned   <- is.na(bin_id[cells_in_bin])
      bin_id[cells_in_bin[unassigned]] <- b
    }

    if (verbose) {
      cat("  Window size:", window_size, "cells\n")
      cat("  Step size:", step_size, "cells\n")
      actual_overlap <- window_size - step_size
      cat(
        "  Actual overlap:", actual_overlap, "cells (",
        round(actual_overlap / window_size * 100, 1), "%)\n", sep = ""
      )
    }
  }

  bin_sizes <- vapply(bin_to_cells, length, integer(1))

  if (verbose) {
    cat("  Bin sizes:\n")
    cat("    Min:", min(bin_sizes), "cells\n")
    cat("    Max:", max(bin_sizes), "cells\n")
    cat("    Mean:", round(mean(bin_sizes), 1), "cells\n")
  }

  if (any(bin_sizes < min_cells_per_bin)) {
    small_bins <- which(bin_sizes < min_cells_per_bin)
    warning(
      "Bins ", paste(small_bins, collapse = ", "),
      " have fewer than ", min_cells_per_bin, " cells. ",
      "Consider reducing `n_bins` or `overlap_ratio`."
    )
  }

  # ============================================================
  # Step 5: Compute bin-averaged expression
  # ============================================================
  if (verbose) cat("\nStep 5: Computing bin-averaged expression...\n")

  bin_means_list <- lapply(seq_len(n_bins), function(b) {
    cells_in_bin <- bin_to_cells[[b]]
    if (length(cells_in_bin) == 0) {
      rep(NA_real_, nrow(expr_mat))
    } else {
      Matrix::rowMeans(expr_mat[, cells_in_bin, drop = FALSE])
    }
  })

  bin_means <- do.call(cbind, bin_means_list)
  rownames(bin_means) <- rownames(expr_mat)
  colnames(bin_means) <- paste0("bin", seq_len(n_bins))

  if (verbose) {
    cat("  ✓ Bin-averaged matrix:", nrow(bin_means), "genes ×", n_bins, "bins\n")
  }

  # ============================================================
  # Step 6: Filter lowly expressed genes
  # ============================================================
  if (filter_genes) {
    if (verbose) cat("\nStep 6: Filtering genes by expression...\n")

    zero_bins <- apply(bin_means, 1, function(gene_expr) {
      sum(gene_expr <= min_expression, na.rm = TRUE)
    })
    zero_bins_pct <- zero_bins / n_bins

    genes_keep  <- zero_bins_pct <= max_zero_bins_pct
    n_filtered  <- sum(!genes_keep)

    if (verbose) {
      cat("  Filtering criteria:\n")
      cat("    Max allowed zero bins:", round(max_zero_bins_pct * 100), "%\n")
      cat("    Min expression threshold:", min_expression, "\n")
      cat(
        "  ✓ Filtered out ", n_filtered, " genes (",
        round(n_filtered / nrow(bin_means) * 100, 1), "%)\n", sep = ""
      )
      cat("  ✓ Retained ", sum(genes_keep), " genes\n", sep = "")
    }

    bin_means <- bin_means[genes_keep, , drop = FALSE]

    gene_filter_stats <- data.frame(
      gene          = rownames(expr_mat),
      n_zero_bins   = zero_bins,
      zero_bins_pct = zero_bins_pct,
      kept          = genes_keep,
      stringsAsFactors = FALSE
    )
  } else {
    if (verbose) cat("\nStep 6: Gene filtering DISABLED\n")
    gene_filter_stats <- NULL
  }

  # ============================================================
  # Step 7: Generate pseudotime coordinates
  # ============================================================
  if (verbose) cat("\nStep 7: Generating pseudotime coordinates...\n")

  x <- seq(x_range[1], x_range[2], length.out = n_bins)

  if (verbose) {
    cat("  ✓ x range: [", x_range[1], ", ", x_range[2], "]\n", sep = "")
    cat(
      "  x values: ",
      paste(round(x[seq_len(min(5, n_bins))], 3), collapse = ", "),
      " ... ", round(x[n_bins], 3), "\n", sep = ""
    )
  }

  # ============================================================
  # Step 8: Compute bin-level pseudotime statistics
  # ============================================================
  bin_pseudotime_median <- vapply(seq_len(n_bins), function(b) {
    cells_in_bin <- bin_to_cells[[b]]
    if (length(cells_in_bin) == 0) {
      NA_real_
    } else {
      median(pt_ord[cells_in_bin])
    }
  }, numeric(1))

  bin_pseudotime_mean <- vapply(seq_len(n_bins), function(b) {
    cells_in_bin <- bin_to_cells[[b]]
    if (length(cells_in_bin) == 0) {
      NA_real_
    } else {
      mean(pt_ord[cells_in_bin])
    }
  }, numeric(1))

  bin_pseudotime_range <- t(vapply(seq_len(n_bins), function(b) {
    cells_in_bin <- bin_to_cells[[b]]
    if (length(cells_in_bin) == 0) {
      c(NA_real_, NA_real_)
    } else {
      c(min(pt_ord[cells_in_bin]), max(pt_ord[cells_in_bin]))
    }
  }, numeric(2)))

  # ============================================================
  # Build result object
  # ============================================================
  if (verbose) {
    cat("\n╔══════════════════════════════════════════════════════════╗\n")
    cat("║                         Summary                          ║\n")
    cat("╚══════════════════════════════════════════════════════════╝\n")
    cat("Total cells: ", ncol(expr_mat), "\n", sep = "")
    cat("Total genes (after filtering): ", nrow(bin_means), "\n", sep = "")
    cat("Number of bins: ", n_bins, "\n", sep = "")
    cat(
      "Binning mode: ",
      ifelse(overlap_ratio == 0,
             "Non-overlapping",
             paste0("Overlapping (", round(overlap_ratio * 100), "%)")),
      "\n", sep = ""
    )
    cat("Pseudotime order: ", pseudotime_order, "\n", sep = "")
    cat("x range: [", x_range[1], ", ", x_range[2], "]\n", sep = "")
    if (filter_genes) {
      cat("Gene filtering: ENABLED\n")
      cat("  Max zero bins: ", round(max_zero_bins_pct * 100), "%\n", sep = "")
      cat(
        "  Genes retained: ", nrow(bin_means), " / ", nrow(expr_mat),
        " (", round(nrow(bin_means) / nrow(expr_mat) * 100, 1), "%)\n",
        sep = ""
      )
    } else {
      cat("Gene filtering: DISABLED\n")
    }
    cat("\n✓ Binning complete!\n\n")
  }

  result <- list(
    bin_means = bin_means,
    x         = x,
    bin_info  = data.frame(
      bin_id           = seq_len(n_bins),
      n_cells          = bin_sizes,
      pseudotime_median = bin_pseudotime_median,
      pseudotime_mean   = bin_pseudotime_mean,
      pseudotime_min    = bin_pseudotime_range[, 1],
      pseudotime_max    = bin_pseudotime_range[, 2],
      x_coord           = x
    ),
    metadata = list(
      pseudotime_col      = pseudotime_col,
      pseudotime_order    = pseudotime_order,
      n_bins              = n_bins,
      overlap_ratio       = overlap_ratio,
      x_range             = x_range,
      filter_genes        = filter_genes,
      max_zero_bins_pct   = max_zero_bins_pct,
      min_expression      = min_expression,
      n_genes_original    = nrow(expr_mat),
      n_genes_filtered    = nrow(bin_means),
      assay               = assay,
      slot                = slot
    ),
    cell_to_bin = data.frame(
      cell      = colnames(expr_mat),
      pseudotime = pt_ord,
      bin_id    = bin_id,
      stringsAsFactors = FALSE
    ),
    bin_to_cells      = bin_to_cells,
    gene_filter_stats = gene_filter_stats
  )

  class(result) <- "scBin"
  return(result)
}
