#' Classify gene expression trends along pseudotime
#'
#' @description
#' `scTrends()` takes the bin-averaged expression matrix produced by
#' [scBin()] and classifies each gene into a temporal trend category
#' along pseudotime using a GAM-based, feature-rich classifier.
#'
#' Each gene is assigned to one of:
#' `"Stable"`, `"Borderline_Stable"`, `"Up"`, `"Down"`,
#' `"Up-Down"`, `"Down-Up"`, `"Complex"`, or `"Insufficient_Data"`.
#'
#' Optionally, permutation-based p-values and FDR-adjusted q-values are
#' computed to assess the significance of the temporal pattern.
#'
#' @param scbin_result An object of class `"scBin"` returned by [scBin()].
#' @param min_cells Minimum number of valid bins required for classification.
#' @param peak_position_margin Fraction of bins excluded at the edges when
#'   checking whether peak/valley lies in the “middle”.
#' @param max_sign_changes_mono Maximum number of sign changes in the
#'   first derivative for a monotone classification.
#' @param gam_k Basis dimension `k` for the GAM smooth term.
#' @param stable_cv_threshold CV threshold for strict “Stable” classification.
#' @param stable_range_threshold Threshold on relative expression range
#'   for “Stable” classification.
#' @param stable_autocorr_threshold Minimum lag-1 autocorrelation supporting
#'   “Stable” classification.
#' @param min_peak_prominence Minimum normalized prominence for peak/valley
#'   patterns.
#' @param min_slope_ratio Minimum normalized slope magnitude around peaks
#'   and valleys.
#' @param total_change_up Minimum normalized net change to label a gene
#'   as `"Up"`.
#' @param total_change_down Maximum normalized net change to label a gene
#'   as `"Down"` (usually negative).
#' @param monotone_consistency Minimum monotonicity score
#'   (fraction of steps in the same direction) for monotone trends.
#' @param min_sign_changes_complex Minimum number of sign changes to support
#'   a `"Complex"` classification.
#' @param complex_variance_threshold Threshold on normalized residual variance
#'   for `"Complex"` trends.
#' @param compute_pvalue Logical; if `TRUE`, perform permutation-based
#'   significance testing per gene.
#' @param n_perm Number of permutations per gene when computing p-values.
#' @param adjust_method Method for multiple testing correction passed to
#'   [stats::p.adjust()], e.g. `"BH"` (default) or `"bonferroni"`.
#'   Use `"none"` to skip adjustment.
#' @param alpha Significance level used to flag significant genes based on
#'   adjusted q-values.
#' @param n_cores Number of cores to use for parallel computation.
#'   Default uses `parallel::detectCores() - 2`.
#' @param use_parallel Logical; if `TRUE`, enable parallel classification.
#' @param return_details Logical; if `TRUE`, store per-gene diagnostics
#'   in the output as well.
#' @param verbose Logical; if `TRUE`, print progress and summary information.
#'
#' @return
#' An object of class `"scTrend"`, which is a list with components:
#' \describe{
#'   \item{results}{Data frame with one row per gene and columns:
#'     `gene`, `trend`, `p_value`, `q_value`, `significant`, plus optional
#'     feature columns (if `return_details = TRUE`).}
#'   \item{summary}{List with global statistics: number of genes, number
#'     of bins, counts per trend category, and number of significant genes.}
#'   \item{parameters}{List of parameter values used by the classifier.}
#'   \item{detailed_results}{(Optional) Raw list of per-gene outputs from
#'     the internal classifier, only if `return_details = TRUE`.}
#' }
#'
#' @importFrom stats p.adjust
#' @importFrom parallel detectCores makeCluster stopCluster clusterSetRNGStream clusterExport clusterEvalQ
#' @importFrom doParallel registerDoParallel
#' @importFrom pbapply pblapply
#' @export
scTrends <- function(scbin_result,
                     # ---------- fixed parameters ----------
                     min_cells = 4,
                     peak_position_margin = 0.15,
                     max_sign_changes_mono = 1,

                     # ---------- core parameters ----------
                     gam_k = 5,
                     stable_cv_threshold = 0.15,
                     stable_range_threshold = 0.12,
                     stable_autocorr_threshold = 0.8,
                     min_peak_prominence = 0.15,
                     min_slope_ratio = 0.08,
                     total_change_up = 0.15,
                     total_change_down = -0.15,
                     monotone_consistency = 0.70,
                     min_sign_changes_complex = 3,
                     complex_variance_threshold = 0.3,

                     # ---------- statistical testing ----------
                     compute_pvalue = TRUE,
                     n_perm = 1000,
                     adjust_method = "BH",
                     alpha = 0.05,

                     # ---------- parallel computation ----------
                     n_cores = parallel::detectCores() - 2,
                     use_parallel = TRUE,

                     # ---------- output control ----------
                     return_details = FALSE,
                     verbose = TRUE) {

  if (!inherits(scbin_result, "scBin")) {
    stop("Input must be an object of class 'scBin' produced by scBin().")
  }

  bin_means <- scbin_result$bin_means
  x         <- scbin_result$x
  n_genes   <- nrow(bin_means)
  n_bins    <- ncol(bin_means)
  gene_names <- rownames(bin_means)

  if (verbose) {
    cat("╔══════════════════════════════════════════════════════════╗\n")
    cat("║           scTrends: Gene Trend Classification            ║\n")
    cat("╚══════════════════════════════════════════════════════════╝\n\n")
    cat("Dataset:", n_genes, "genes ×", n_bins, "bins\n")
    cat(
      "Parameters: Stable(CV<", stable_cv_threshold,
      "), Peak(Prom>", min_peak_prominence,
      "), Mono(Up>", total_change_up,
      ", Down<", total_change_down, ")\n", sep = ""
    )
    cat(
      "P-value: ",
      ifelse(compute_pvalue,
             paste0("ENABLED (", n_perm, " permutations)"),
             "DISABLED"),
      "\n\n", sep = ""
    )
  }

  # Parallel setup
  cl <- NULL
  if (use_parallel && !is.na(n_cores) && n_cores > 1) {
    cl <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)
    if (compute_pvalue) {
      parallel::clusterSetRNGStream(cl, iseed = 123)
    }
    parallel::clusterExport(cl, "classify_single_gene_v4", envir = environment())
    parallel::clusterEvalQ(cl, library(mgcv))
    if (verbose) cat("Parallel cluster initialized (", n_cores, " cores)\n\n", sep = "")
  }

  # Classification
  if (verbose) cat("Classifying genes...\n")

  if (!is.null(cl)) {
    results <- pbapply::pblapply(
      seq_len(n_genes),
      function(i) {
        classify_single_gene_v4(
          gm = bin_means[i, ], x = x,
          min_cells = min_cells,
          peak_position_margin = peak_position_margin,
          max_sign_changes_mono = max_sign_changes_mono,
          gam_k = gam_k,
          stable_cv_threshold = stable_cv_threshold,
          stable_range_threshold = stable_range_threshold,
          stable_autocorr_threshold = stable_autocorr_threshold,
          min_peak_prominence = min_peak_prominence,
          min_slope_ratio = min_slope_ratio,
          total_change_up = total_change_up,
          total_change_down = total_change_down,
          monotone_consistency = monotone_consistency,
          min_sign_changes_complex = min_sign_changes_complex,
          complex_variance_threshold = complex_variance_threshold,
          compute_pvalue = compute_pvalue,
          n_perm = n_perm,
          return_details = return_details,
          verbose = FALSE
        )
      },
      cl = cl
    )
    parallel::stopCluster(cl)
  } else {
    results <- vector("list", n_genes)
    for (i in seq_len(n_genes)) {
      if (verbose) {
        pct <- round(i / n_genes * 100)
        cat(
          "\r[",
          paste0(rep("█", pct %/% 2), collapse = ""),
          paste0(rep("░", 50 - pct %/% 2), collapse = ""),
          "] ", i, "/", n_genes, " (", pct, "%)",
          sep = ""
        )
        flush.console()
      }
      results[[i]] <- classify_single_gene_v4(
        gm = bin_means[i, ], x = x,
        min_cells = min_cells,
        peak_position_margin = peak_position_margin,
        max_sign_changes_mono = max_sign_changes_mono,
        gam_k = gam_k,
        stable_cv_threshold = stable_cv_threshold,
        stable_range_threshold = stable_range_threshold,
        stable_autocorr_threshold = stable_autocorr_threshold,
        min_peak_prominence = min_peak_prominence,
        min_slope_ratio = min_slope_ratio,
        total_change_up = total_change_up,
        total_change_down = total_change_down,
        monotone_consistency = monotone_consistency,
        min_sign_changes_complex = min_sign_changes_complex,
        complex_variance_threshold = complex_variance_threshold,
        compute_pvalue = compute_pvalue,
        n_perm = n_perm,
        return_details = return_details,
        verbose = FALSE
      )
    }
    if (verbose) cat("\n")
  }

  if (verbose) cat("✓ Classification complete!\n\n")

  # Extract main results
  trends <- vapply(
    results,
    function(r) {
      if (is.null(r$trend) || is.na(r$trend)) "Insufficient_Data" else r$trend
    },
    character(1)
  )

  p_values <- vapply(
    results,
    function(r) if (is.null(r$p_value)) NA_real_ else r$p_value,
    numeric(1)
  )

  # FDR correction
  if (compute_pvalue && adjust_method != "none" && any(!is.na(p_values))) {
    valid_idx <- !is.na(p_values)
    q_values <- rep(NA_real_, n_genes)
    q_values[valid_idx] <- stats::p.adjust(p_values[valid_idx], method = adjust_method)
    significant <- rep(NA, n_genes)
    significant[valid_idx] <- q_values[valid_idx] < alpha
  } else {
    q_values   <- rep(NA_real_, n_genes)
    significant <- rep(NA, n_genes)
  }

  # Build results data frame
  result_df <- data.frame(
    gene        = gene_names,
    trend       = trends,
    p_value     = p_values,
    q_value     = q_values,
    significant = significant,
    stringsAsFactors = FALSE,
    row.names   = NULL
  )

  if (return_details) {
    result_df$cv <- vapply(
      results,
      function(r) if (!is.null(r$details$cv)) r$details$cv else NA_real_,
      numeric(1)
    )
    result_df$peak_prominence <- vapply(
      results,
      function(r) if (!is.null(r$details$peak_prominence)) r$details$peak_prominence else NA_real_,
      numeric(1)
    )
    result_df$total_change <- vapply(
      results,
      function(r) if (!is.null(r$details$total_change)) r$details$total_change else NA_real_,
      numeric(1)
    )
    result_df$classification_reason <- vapply(
      results,
      function(r) if (!is.null(r$details$classification_reason)) r$details$classification_reason else NA_character_,
      character(1)
    )
  }

  # Summary counts
  trend_counts <- table(trends)

  if (verbose) {
    cat("╔══════════════════════════════════════════════════════════╗\n")
    cat("║                     Results Summary                      ║\n")
    cat("╚══════════════════════════════════════════════════════════╝\n\n")

    ordered_trends <- c(
      "Stable", "Borderline_Stable",
      "Up", "Down",
      "Up-Down", "Down-Up",
      "Complex", "Insufficient_Data"
    )

    for (t in ordered_trends) {
      if (t %in% names(trend_counts)) {
        n   <- trend_counts[t]
        pct <- round(n / n_genes * 100, 1)
        n_sig <- if (compute_pvalue) {
          sum(trends == t & significant, na.rm = TRUE)
        } else {
          0
        }

        display_name <- if (t == "Borderline_Stable") {
          paste0(t, " *")
        } else {
          t
        }

        cat(sprintf(
          "  %-22s: %5d (%5.1f%%) | %4d sig\n",
          display_name, n, pct, n_sig
        ))
      }
    }

    if ("Borderline_Stable" %in% names(trend_counts)) {
      cat("\n  * Borderline_Stable: genes with weak trends not meeting\n")
      cat("    strict Stable criteria but showing low variability\n")
    }

    cat("\nTotal:", n_genes, "genes")
    if (compute_pvalue) {
      cat(" |", sum(significant, na.rm = TRUE), "significant (",
          round(sum(significant, na.rm = TRUE) / n_genes * 100, 1),
          "%)", sep = "")
    }
    cat("\n\n✓ Analysis complete!\n\n")
  }

  output <- list(
    results = result_df,
    summary = list(
      n_genes        = n_genes,
      n_bins         = n_bins,
      trend_counts   = trend_counts,
      n_significant  = if (compute_pvalue) sum(significant, na.rm = TRUE) else NA,
      stable_types   = list(
        strict_stable      = sum(trends == "Stable", na.rm = TRUE),
        borderline_stable  = sum(trends == "Borderline_Stable", na.rm = TRUE)
      )
    ),
    parameters = list(
      min_cells = min_cells,
      peak_position_margin = peak_position_margin,
      max_sign_changes_mono = max_sign_changes_mono,
      gam_k = gam_k,
      stable_cv_threshold = stable_cv_threshold,
      stable_range_threshold = stable_range_threshold,
      stable_autocorr_threshold = stable_autocorr_threshold,
      min_peak_prominence = min_peak_prominence,
      min_slope_ratio = min_slope_ratio,
      total_change_up = total_change_up,
      total_change_down = total_change_down,
      monotone_consistency = monotone_consistency,
      min_sign_changes_complex = min_sign_changes_complex,
      complex_variance_threshold = complex_variance_threshold,
      compute_pvalue = compute_pvalue,
      n_perm = n_perm,
      adjust_method = adjust_method,
      alpha = alpha
    )
  )

  if (return_details) {
    output$detailed_results <- results
  }

  class(output) <- "scTrends"
  return(output)
}
