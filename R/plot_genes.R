#' Plot selected genes along pseudotime
#'
#' @description
#' `plot_genes()` visualizes the pseudotime trajectories of one or more
#' user-specified genes based on a bin-averaged expression matrix
#' (e.g. `scBin()` output). Each gene is plotted as a separate line,
#' with distinct colors and shared x-axis (pseudotime).
#'
#' This function is intentionally style-light: it only defines the core
#' layers (lines and colors) and returns a `ggplot` object that can be
#' further customized with themes, labels, coordinate limits, etc.
#'
#' @param bin_means A numeric matrix of bin-averaged expression values
#'   with genes in rows and pseudotime bins in columns. Row names must
#'   correspond to gene identifiers.
#' @param x Numeric vector of pseudotime coordinates corresponding to
#'   the columns of `bin_means`.
#' @param genes Character vector of gene names to plot. All genes must
#'   be present in the row names of `bin_means`.
#' @param color Optional color specification for genes. Can be:
#'   \itemize{
#'     \item `NULL` (default): a simple default palette is used;
#'     \item a character vector of colors with length \eqn{\ge} number of genes;
#'     \item a named character vector whose names match gene IDs, in which
#'           case colors are matched by gene name.
#'   }
#' @param alpha Numeric transparency level for the lines (0–1).
#' @param size Numeric line width for the gene trajectories.
#'
#' @return
#' A `ggplot` object showing the trajectories of the selected genes
#' along pseudotime. Users are encouraged to add their own themes,
#' labels, and coordinate settings.
#'
#' @examples
#' \dontrun{
#'   # Suppose you have:
#'   #   sb     <- scBin(seurat_obj, pseudotime_col = "pseudotime")
#'   #   genes  <- c("WT1", "NPHS1", "NPHS2")
#'
#'   p <- plot_genes(
#'     bin_means = sb$bin_means,
#'     x         = sb$x,
#'     genes     = genes
#'   ) +
#'     ggplot2::labs(
#'       x = "Pseudotime (0 = Early, 1 = Late)",
#'       y = "Expression level"
#'     ) +
#'     ggplot2::theme_minimal()
#'
#'   print(p)
#' }
#'
#' @export
plot_genes <- function(
    bin_means, x,
    genes,
    color = NULL,
    alpha = 0.3,
    size = 0.5
) {
  # ------------------------------------------------------------
  # 0. Basic checks
  # ------------------------------------------------------------
  if (!is.matrix(bin_means) && !inherits(bin_means, "Matrix")) {
    stop("`bin_means` must be a numeric matrix with genes in rows.")
  }
  if (is.null(rownames(bin_means))) {
    stop("`bin_means` must have rownames corresponding to gene IDs.")
  }
  if (length(x) != ncol(bin_means)) {
    stop("Length of `x` must match the number of columns in `bin_means`.")
  }

  # ------------------------------------------------------------
  # 1. Check and normalize `genes`
  # ------------------------------------------------------------
  if (is.factor(genes)) {
    genes <- as.character(genes)
  }
  if (!is.character(genes)) {
    stop("`genes` should be a character vector of gene names.")
  }
  genes <- unique(genes)

  if (!all(genes %in% rownames(bin_means))) {
    missing_genes <- setdiff(genes, rownames(bin_means))
    stop(
      "The following genes are not found in `bin_means`: ",
      paste(missing_genes, collapse = ", ")
    )
  }

  # ------------------------------------------------------------
  # 2. Extract data for selected genes
  # ------------------------------------------------------------
  gene_data <- bin_means[genes, , drop = FALSE]
  n_genes   <- nrow(gene_data)
  n_bins    <- ncol(gene_data)

  df_plot <- data.frame(
    x         = rep(x, times = n_genes),
    expression = as.vector(t(gene_data)),
    gene      = rep(rownames(gene_data), each = n_bins),
    stringsAsFactors = FALSE
  )

  # ------------------------------------------------------------
  # 3. Color handling
  # ------------------------------------------------------------
  # Default colors (simple palette; users can override via `color`)
  default_colors <- c(
    "red", "blue", "green", "purple", "orange",
    "brown", "pink", "cyan", "yellow", "gray"
  )

  if (is.null(color)) {
    # extend default palette if needed
    if (n_genes > length(default_colors)) {
      color <- rep(default_colors, length.out = n_genes)
    } else {
      color <- default_colors[seq_len(n_genes)]
    }
    names(color) <- rownames(gene_data)
  } else {
    # if named vector, match by gene name
    if (!is.null(names(color))) {
      if (!all(rownames(gene_data) %in% names(color))) {
        missing_cols <- setdiff(rownames(gene_data), names(color))
        stop(
          "Color vector is named, but missing colors for genes: ",
          paste(missing_cols, collapse = ", ")
        )
      }
      color <- color[rownames(gene_data)]
    } else {
      # unnamed color vector: must be at least n_genes
      if (length(color) < n_genes) {
        stop("The number of colors must be at least the number of genes.")
      }
      color <- rep(color, length.out = n_genes)
      names(color) <- rownames(gene_data)
    }
  }

  # ------------------------------------------------------------
  # 4. Build ggplot (minimal styling)
  # ------------------------------------------------------------
  p <- ggplot2::ggplot(
    df_plot,
    ggplot2::aes(x = .data$x, y = .data$expression, color = .data$gene)
  ) +
    ggplot2::geom_line(alpha = alpha, linewidth = size) +
    ggplot2::scale_color_manual(values = color)

  # Note: no theme, labels, or coord limits are imposed here,
  # so users can freely customize the plot:
  #   + ggplot2::theme_minimal()
  #   + ggplot2::labs(...)
  #   + ggplot2::coord_cartesian(...)

  return(p)
}
