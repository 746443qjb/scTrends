#' scTrends: Trend detection for single-cell pseudotime analyses
#'
#' @description
#' The **scTrends** package provides a workflow to:
#' \enumerate{
#'   \item Bin single-cell expression profiles along pseudotime using [scBin()].
#'   \item Automatically tune trend-detection parameters via [scPS()].
#'   \item Classify genes into temporal trend categories using [scTrends()].
#'   \item Visualize gene sets or individual genes along pseudotime with
#'         [plot_trends()] and [plot_genes()].
#' }
#'
#' The typical workflow is:
#' \preformatted{
#'   sb  <- scBin(seurat_obj, pseudotime_col = "pseudotime")
#'   ps  <- scPS(sb)
#'   tr  <- do.call(scTrends, c(list(scbin_result = sb), ps$optimal_params))
#'
#'   # Plot all significant Up genes
#'   p1 <- plot_trend(tr, bin_means = sb$bin_means, x = sb$x,
#'                    trend_type = "Up", sig_metric = "q_value")
#'
#'   # Plot a few marker genes
#'   p2 <- plot_genes(sb$bin_means, sb$x,
#'                    genes = c("WT1", "NPHS1", "NPHS2"))
#' }
#'
#' @keywords internal
"_PACKAGE"
