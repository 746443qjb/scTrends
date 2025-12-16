#' Plot gene expression trends along pseudotime
#'
#' @description
#' `plot_trends()` visualizes gene expression trajectories along pseudotime
#' based on bin-averaged expression (e.g. from [scBin()]) and trend
#' classification results (e.g. from [scTrends()]). It can:
#' \itemize{
#'   \item filter genes by trend type (e.g. `"Up"`, `"Down"`, `"Complex"`);
#'   \item optionally filter by statistical significance (`p_value` or `q_value`);
#'   \item plot individual gene trajectories;
#'   \item overlay the mean trajectory and an uncertainty band (standard error).
#' }
#'
#' The function is designed to be style-agnostic: it returns a ggplot object
#' that can be further customized with user-defined themes, labels, etc.
#'
#' @param trend_result Either:
#'   \itemize{
#'     \item an object of class `"scTrend"` produced by [scTrends()], or
#'     \item a `data.frame` with at least the columns `gene`, `trend`,
#'           and optionally `p_value`, `q_value`.
#'   }
#' @param bin_means A numeric matrix of bin-averaged expression values with
#'   genes in rows and pseudotime bins in columns (e.g. `scbin_result$bin_means`).
#'   Row names must match the `gene` column in `trend_result`.
#' @param x Numeric vector of pseudotime coordinates corresponding to
#'   columns of `bin_means` (e.g. `scbin_result$x`).
#' @param trend_type Character scalar specifying which trend to plot.
#'   Use `"none"` to ignore trend categories and include all selected genes.
#' @param sig_metric Significance metric to use for gene filtering.
#'   One of `"p_value"`, `"q_value"`, or `"none"` (no significance filtering).
#' @param max_genes Integer. Maximum number of genes to plot after filtering.
#' @param sig_cutoff Numeric threshold for `p_value` or `q_value` when
#'   `sig_metric != "none"`. Default is `0.05`.
#' @param show_individual Logical; if `TRUE`, plot individual gene trajectories.
#' @param show_mean Logical; if `TRUE`, overlay the mean trajectory across
#'   selected genes.
#' @param show_se Logical; if `TRUE` and `show_mean = TRUE`, add a ribbon
#'   representing mean ± `se_scale` × standard error.
#' @param individual_color Color for individual gene trajectories.
#' @param individual_alpha Alpha (transparency) for individual trajectories.
#' @param individual_size Line width for individual trajectories.
#' @param mean_color Color for the mean trajectory.
#' @param mean_size Line width for the mean trajectory.
#' @param se_fill Fill color for the standard error ribbon.
#' @param se_alpha Alpha (transparency) for the standard error ribbon.
#' @param se_scale Numeric multiplier for standard error (e.g. 1.96 ≈ 95% band).
#'
#' @return
#' A `ggplot` object that can be further customized with additional
#' `ggplot2` layers, themes, and labels.
#'
#' @examples
#' \dontrun{
#'   # Suppose you have:
#'   #   sb      <- scBin(seurat_obj, pseudotime_col = "pseudotime")
#'   #   trends  <- scTrends(sb)
#'
#'   p <- plot_trends(
#'     trend_result = trends,          # or trends$results
#'     bin_means    = sb$bin_means,
#'     x            = sb$x,
#'     trend_type   = "Up",
#'     sig_metric   = "q_value",
#'     sig_cutoff   = 0.05,
#'     max_genes    = 50
#'   ) +
#'     ggplot2::labs(
#'       title = "Up-regulated genes",
#'       x = "Pseudotime",
#'       y = "Expression"
#'     ) +
#'     ggplot2::theme_minimal()
#'
#'   print(p)
#' }
#'
#' @export
plot_trends <- function(
    trend_result, bin_means, x,
    trend_type = "none",
    sig_metric = c("p_value", "q_value", "none"),
    max_genes = 30,
    sig_cutoff = 0.05,
    show_individual = TRUE,
    show_mean = TRUE,
    show_se = TRUE,
    individual_color = "gray60",
    individual_alpha = 0.3,
    individual_size = 0.5,
    mean_color = "darkblue",
    mean_size = 1.5,
    se_fill = "steelblue",
    se_alpha = 0.3,
    se_scale = 1.96
) {
  sig_metric <- match.arg(sig_metric)

  # ------------------------------------------------------------
  # 1. Normalize trend_result: accept scTrend or data.frame
  # ------------------------------------------------------------
  if (inherits(trend_result, "scTrend")) {
    tr_df <- trend_result$results
  } else if (is.data.frame(trend_result)) {
    tr_df <- trend_result
  } else {
    stop("`trend_result` must be an 'scTrend' object or a data.frame.")
  }

  required_cols <- c("gene", "trend")
  missing_cols <- setdiff(required_cols, colnames(tr_df))
  if (length(missing_cols) > 0) {
    stop(
      "trend_result is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

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
  # 2. Gene filtering
  # ------------------------------------------------------------
  sig_genes <- tr_df

  # (A) filter by trend_type (if not "none")
  available_trends <- unique(tr_df$trend[!is.na(tr_df$trend)])

  if (trend_type != "none") {
    if (!trend_type %in% available_trends) {
      stop(
        "Trend type '", trend_type, "' not found. Available trends: ",
        paste(available_trends, collapse = ", ")
      )
    }
    sig_genes <- dplyr::filter(sig_genes, .data$trend == trend_type)
  }

  # (B) filter by significance metric (if not "none")
  if (sig_metric != "none") {
    if (!sig_metric %in% colnames(sig_genes)) {
      stop(
        "Selected sig_metric ('", sig_metric,
        "') not found in trend_result."
      )
    }

    if (sig_metric == "p_value") {
      sig_genes <- sig_genes |>
        dplyr::filter(.data$p_value < sig_cutoff) |>
        dplyr::arrange(.data$p_value)
    } else if (sig_metric == "q_value") {
      sig_genes <- sig_genes |>
        dplyr::filter(.data$q_value < sig_cutoff) |>
        dplyr::arrange(.data$q_value)
    }
  }

  # limit number of genes
  sig_genes <- dplyr::slice_head(sig_genes, n = max_genes)

  if (nrow(sig_genes) == 0) {
    stop("No genes selected under current filter settings.")
  }

  # ensure selected genes exist in bin_means
  missing_in_matrix <- setdiff(sig_genes$gene, rownames(bin_means))
  if (length(missing_in_matrix) > 0) {
    warning(
      length(missing_in_matrix), " selected genes were not found in `bin_means` ",
      "and will be removed from plotting."
    )
    sig_genes <- sig_genes[!(sig_genes$gene %in% missing_in_matrix), , drop = FALSE]
  }

  if (nrow(sig_genes) == 0) {
    stop("No genes available for plotting after matching with `bin_means`.")
  }

  # ------------------------------------------------------------
  # 3. Prepare data for plotting
  # ------------------------------------------------------------
  gene_data <- bin_means[sig_genes$gene, , drop = FALSE]

  n_genes <- nrow(gene_data)
  n_bins  <- ncol(gene_data)

  df_plot <- data.frame(
    x         = rep(x, times = n_genes),
    expression = as.vector(t(gene_data)),
    gene      = rep(rownames(gene_data), each = n_bins),
    stringsAsFactors = FALSE
  )

  # mean and standard error across genes
  mean_expr <- colMeans(gene_data, na.rm = TRUE)
  se_expr   <- apply(
    gene_data, 2,
    function(col) stats::sd(col, na.rm = TRUE) / sqrt(sum(!is.na(col)))
  )

  df_mean <- data.frame(
    x        = x,
    mean_expr = mean_expr,
    se       = se_expr
  )

  # ------------------------------------------------------------
  # 4. Build ggplot
  # ------------------------------------------------------------
  p <- ggplot2::ggplot()

  # Individual gene trajectories
  if (show_individual) {
    p <- p +
      ggplot2::geom_line(
        data = df_plot,
        ggplot2::aes(x = .data$x, y = .data$expression, group = .data$gene),
        color = individual_color,
        alpha = individual_alpha,
        linewidth = individual_size
      )
  }

  # standard error ribbon (around mean)
  if (show_mean && show_se) {
    p <- p +
      ggplot2::geom_ribbon(
        data = df_mean,
        ggplot2::aes(
          x = .data$x,
          ymin = .data$mean_expr - se_scale * .data$se,
          ymax = .data$mean_expr + se_scale * .data$se
        ),
        fill  = se_fill,
        alpha = se_alpha
      )
  }

  # mean trajectory
  if (show_mean) {
    p <- p +
      ggplot2::geom_line(
        data = df_mean,
        ggplot2::aes(x = .data$x, y = .data$mean_expr),
        color    = mean_color,
        linewidth = mean_size
      )
  }

  # automatic title
  title <- if (trend_type == "none") {
    paste0("Selected Genes (n = ", n_genes, ")")
  } else {
    paste0("Trend: ", trend_type, " (n = ", n_genes, ")")
  }

  p <- p + ggplot2::ggtitle(title)

  return(p)
}
