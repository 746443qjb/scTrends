#' Classify the temporal trend of a single gene
#'
#' @description
#' Internal helper used by [scTrends()] to classify the temporal trend
#' of a single gene based on bin-averaged expression along pseudotime.
#' It fits a GAM curve, extracts a set of shape and variability features,
#' and assigns the gene to one of several trend classes:
#' `"Stable"`, `"Borderline_Stable"`, `"Up"`, `"Down"`,
#' `"Up-Down"`, `"Down-Up"`, `"Complex"`, or `"Insufficient_Data"`.
#'
#' This function is not intended to be called directly by end users.
#'
#' @param gm Numeric vector of gene expression values across pseudotime bins.
#' @param x Numeric vector of pseudotime coordinates of the same length as `gm`.
#' @param min_cells Minimum number of valid bins required for classification.
#' @param peak_position_margin Fraction of bins excluded at each end when
#'   checking whether a peak/valley lies “in the middle”.
#' @param max_sign_changes_mono Maximum allowed number of sign changes in
#'   the first derivative for a gene to be considered monotone.
#' @param gam_k Basis dimension `k` used in the GAM smoother.
#' @param stable_cv_threshold Coefficient of variation (CV) threshold for
#'   classifying a gene as strictly stable.
#' @param stable_range_threshold Threshold on relative range of expression
#'   for stability classification.
#' @param stable_autocorr_threshold Minimum lag-1 autocorrelation required
#'   to support a “Stable” classification (if available).
#' @param min_peak_prominence Minimum normalized prominence (relative to
#'   global range) required to call a peak/valley.
#' @param min_slope_ratio Minimum normalized slope magnitude required on
#'   each side of a peak/valley.
#' @param total_change_up Minimum normalized net change (end – start) to
#'   classify as `"Up"`.
#' @param total_change_down Maximum normalized net change to classify as
#'   `"Down"` (typically negative).
#' @param monotone_consistency Minimum fraction of steps in the same
#'   direction required for monotone trends.
#' @param min_sign_changes_complex Minimum number of sign changes needed
#'   for a `"Complex"` classification.
#' @param complex_variance_threshold Threshold on normalized residual variance
#'   used to support `"Complex"` classification.
#' @param compute_pvalue Logical; if `TRUE`, compute a permutation-based
#'   p-value for the GAM fit.
#' @param n_perm Number of permutations for p-value estimation.
#' @param return_details Logical; if `TRUE`, return a rich list of
#'   intermediate features and diagnostics in `$details`.
#' @param verbose Logical; if `TRUE`, print messages for debugging.
#'
#' @return
#' A list with elements:
#' \itemize{
#'   \item `trend`: Character scalar with the assigned trend class.
#'   \item `p_value`: Numeric scalar with the permutation p-value
#'         (or `NA` if not computed).
#'   \item `details`: Optional list of diagnostics (only if
#'         `return_details = TRUE`).
#' }
#'
#' @keywords internal
#' @noRd
classify_single_gene_v4 <- function(
    gm, x,
    # ---------- fixed parameters (not recommended to change) ----------
    min_cells = 4,
    peak_position_margin = 0.15,
    max_sign_changes_mono = 1,

    # ---------- core parameters (to be tuned by higher-level methods) ----------
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

    # ---------- output control ----------
    return_details = FALSE,
    verbose = FALSE
) {
  # ------------------------------------------------------------
  # Helper functions
  # ------------------------------------------------------------
  safe_return <- function(trend, pval, details = NULL) {
    result <- list(
      trend = as.character(trend),
      p_value = if (compute_pvalue) as.numeric(pval) else NA_real_
    )
    if (return_details && !is.null(details)) {
      result$details <- details
    }
    return(result)
  }

  calc_slope_robust <- function(y, x_vals) {
    if (length(y) < 2 || length(x_vals) < 2) return(NA_real_)
    tryCatch({
      fit <- stats::lm(y ~ x_vals)
      unname(stats::coef(fit)[2])
    }, error = function(e) NA_real_)
  }

  calc_sign_changes <- function(pred) {
    deriv1 <- diff(pred)
    sign_deriv <- sign(deriv1)
    sum(abs(diff(sign_deriv)) == 2, na.rm = TRUE)
  }

  # ------------------------------------------------------------
  # Phase 1: data cleaning and pre-processing
  # ------------------------------------------------------------
  gm <- as.numeric(gm)
  ok <- !is.na(gm) & !is.na(x)
  if (sum(ok) < min_cells) {
    return(safe_return("Insufficient_Data", NA_real_))
  }
  gm <- gm[ok]
  x2 <- x[ok]
  n_bins <- length(gm)

  # ------------------------------------------------------------
  # Phase 2: Stable classification (multi-criteria)
  # ------------------------------------------------------------
  mean_expr <- mean(gm)
  if (abs(mean_expr) < 1e-6) mean_expr <- 1e-6
  cv <- stats::sd(gm) / abs(mean_expr)

  rng <- max(gm) - min(gm)
  baseline <- stats::median(gm)
  if (baseline < 1e-6) baseline <- 1e-6
  relative_change <- rng / baseline

  if (n_bins >= 3) {
    autocorr <- tryCatch({
      stats::cor(gm[-n_bins], gm[-1], method = "spearman")
    }, error = function(e) NA_real_)
  } else {
    autocorr <- NA_real_
  }

  is_stable <- !is.na(cv) && !is.na(relative_change) &&
    cv < stable_cv_threshold &&
    relative_change < stable_range_threshold

  if (!is.na(autocorr)) {
    is_stable <- is_stable && (autocorr > stable_autocorr_threshold)
  }

  if (is_stable) {
    # P-value under a linear null (optional)
    p_value_stable <- NA_real_
    if (compute_pvalue && n_perm > 0) {
      fit_linear <- stats::lm(gm ~ x2)
      obs_deviance <- sum(stats::residuals(fit_linear)^2)

      null_deviances <- replicate(n_perm, {
        gm_null <- mean(gm) + stats::rnorm(n_bins, 0, stats::sd(gm))
        fit_null <- stats::lm(gm_null ~ x2)
        sum(stats::residuals(fit_null)^2)
      })

      p_value_stable <- mean(null_deviances <= obs_deviance)
    }

    if (verbose) {
      cat("  → Classified as Stable (CV:",
          round(cv, 3), ", RelChange:", round(relative_change, 3), ")\n")
    }

    details <- if (return_details) {
      list(
        cv = cv,
        relative_change = relative_change,
        autocorr = autocorr,
        reason = "Low variability"
      )
    } else NULL

    return(safe_return("Stable", p_value_stable, details))
  }

  # ------------------------------------------------------------
  # Phase 3: GAM fitting
  # ------------------------------------------------------------
  fit <- tryCatch({
    dat <- data.frame(y = gm, x = x2)
    mgcv::gam(y ~ s(x, k = min(gam_k, length(x2) - 2)), data = dat)
  }, error = function(e) NULL, warning = function(w) NULL)

  if (is.null(fit)) {
    details <- if (return_details) {
      list(reason = "GAM fitting failed", n_bins = n_bins)
    } else NULL
    return(safe_return("Insufficient_Data", NA_real_, details))
  }

  pred <- tryCatch({
    stats::predict(fit)
  }, error = function(e) NULL)

  if (is.null(pred) || any(is.na(pred))) {
    return(safe_return("Insufficient_Data", NA_real_))
  }

  # ------------------------------------------------------------
  # Phase 4: permutation-based significance test (one-sided)
  # ------------------------------------------------------------
  p_value <- NA_real_
  if (compute_pvalue && n_perm > 0) {
    tryCatch({
      null_dev  <- fit$null.deviance
      resid_dev <- fit$deviance

      if (!is.na(null_dev) && !is.na(resid_dev) && null_dev > 0) {
        stat_real <- (null_dev - resid_dev) / null_dev

        stat_null <- replicate(n_perm, {
          gm_perm <- sample(gm)
          fit_perm <- tryCatch({
            mgcv::gam(
              gm_perm ~ s(x2, k = min(gam_k, length(x2) - 2)),
              data = data.frame(y = gm_perm, x = x2)
            )
          }, error = function(e) NULL)

          if (is.null(fit_perm)) return(NA_real_)

          null_d  <- fit_perm$null.deviance
          resid_d <- fit_perm$deviance
          if (is.na(null_d) || is.na(resid_d) || null_d == 0) return(NA_real_)
          (null_d - resid_d) / null_d
        })

        stat_null <- stat_null[!is.na(stat_null)]

        if (length(stat_null) > 0) {
          p_value <- mean(stat_null >= stat_real)
        }
      }
    }, error = function(e) {
      p_value <- NA_real_
    })
  }

  # ------------------------------------------------------------
  # Phase 5: feature extraction
  # ------------------------------------------------------------
  tryCatch({
    start_val <- pred[1]
    end_val   <- pred[length(pred)]
    max_val   <- max(pred)
    min_val   <- min(pred)
    max_idx   <- which.max(pred)[1]
    min_idx   <- which.min(pred)[1]

    pred_rng <- max_val - min_val
    if (pred_rng < 1e-6) pred_rng <- 1e-6

    # Peak / valley prominence
    peak_prominence   <- (max_val - max(start_val, end_val)) / pred_rng
    valley_prominence <- (min(start_val, end_val) - min_val) / pred_rng

    # Peak / valley position
    edge_bins <- round(n_bins * peak_position_margin)
    peak_in_middle   <- (max_idx > edge_bins) && (max_idx < (n_bins - edge_bins))
    valley_in_middle <- (min_idx > edge_bins) && (min_idx < (n_bins - edge_bins))

    # Slopes around peak
    left_slope_norm  <- 0
    right_slope_norm <- 0
    if (!is.na(max_idx) && max_idx > 1 && max_idx < length(pred)) {
      left_slope <- calc_slope_robust(pred[1:max_idx], x2[1:max_idx])
      right_slope <- calc_slope_robust(pred[max_idx:length(pred)],
                                       x2[max_idx:length(x2)])
      if (!is.na(left_slope))  left_slope_norm  <- left_slope / pred_rng
      if (!is.na(right_slope)) right_slope_norm <- right_slope / pred_rng
    }

    # Slopes around valley
    valley_left_norm  <- 0
    valley_right_norm <- 0
    if (!is.na(min_idx) && min_idx > 1 && min_idx < length(pred)) {
      valley_left_slope <- calc_slope_robust(pred[1:min_idx], x2[1:min_idx])
      valley_right_slope <- calc_slope_robust(pred[min_idx:length(pred)],
                                              x2[min_idx:length(x2)])
      if (!is.na(valley_left_slope))  valley_left_norm  <- valley_left_slope / pred_rng
      if (!is.na(valley_right_slope)) valley_right_norm <- valley_right_slope / pred_rng
    }

    # Global slope (only sign is used)
    global_slope      <- calc_slope_robust(pred, x2)
    global_slope_sign <- sign(global_slope)

    # Net change
    total_change <- (end_val - start_val) / pred_rng

    # Number of direction changes
    sign_changes <- calc_sign_changes(pred)

    # Monotonicity score
    diffs <- diff(pred)
    up_ratio   <- mean(diffs > 0, na.rm = TRUE)
    down_ratio <- mean(diffs < 0, na.rm = TRUE)
    monotone_score <- max(up_ratio, down_ratio)

    # Residual variance (normalized)
    residuals <- gm - pred
    residual_variance <- stats::var(residuals, na.rm = TRUE) /
      stats::var(gm, na.rm = TRUE)

    # ----------------------------------------------------------
    # Phase 6: trend classification
    # ----------------------------------------------------------
    trend <- "Complex"
    classification_reason <- "Default"

    # Priority 1: Peak-like / valley-like patterns
    if (!is.na(peak_prominence) && !is.na(left_slope_norm) &&
        !is.na(right_slope_norm) &&
        peak_prominence > min_peak_prominence &&
        peak_in_middle &&
        left_slope_norm >  min_slope_ratio &&
        right_slope_norm < -min_slope_ratio) {

      trend <- "Up-Down"
      classification_reason <- paste0(
        "Peak at ", round(max_idx / n_bins * 100), "%, prominence=",
        round(peak_prominence, 2)
      )

    } else if (!is.na(valley_prominence) && !is.na(valley_left_norm) &&
               !is.na(valley_right_norm) &&
               valley_prominence > min_peak_prominence &&
               valley_in_middle &&
               valley_left_norm < -min_slope_ratio &&
               valley_right_norm >  min_slope_ratio) {

      trend <- "Down-Up"
      classification_reason <- paste0(
        "Valley at ", round(min_idx / n_bins * 100), "%, prominence=",
        round(valley_prominence, 2)
      )

    } else if (!is.na(total_change) && !is.na(monotone_score) &&
               !is.na(sign_changes) && !is.na(global_slope_sign)) {

      # Priority 2: Monotone patterns (with dynamic thresholds)
      high_consistency_threshold <- min(0.95, monotone_consistency * 1.15)
      relaxed_up_threshold   <- total_change_up   * 0.8
      relaxed_down_threshold <- total_change_down * 0.8

      # Up
      if ((total_change > total_change_up ||
           (total_change > relaxed_up_threshold &&
            monotone_score > high_consistency_threshold)) &&
          monotone_score > monotone_consistency &&
          sign_changes <= max_sign_changes_mono &&
          global_slope_sign > 0) {

        trend <- "Up"
        classification_reason <- paste0(
          "Monotone increase, total_change=", round(total_change, 2),
          ", consistency=", round(monotone_score, 2)
        )

        # Down
      } else if ((total_change < total_change_down ||
                  (total_change < relaxed_down_threshold &&
                   monotone_score > high_consistency_threshold)) &&
                 monotone_score > monotone_consistency &&
                 sign_changes <= max_sign_changes_mono &&
                 global_slope_sign < 0) {

        trend <- "Down"
        classification_reason <- paste0(
          "Monotone decrease, total_change=", round(total_change, 2),
          ", consistency=", round(monotone_score, 2)
        )

      } else {
        # Priority 3: Complex vs Borderline_Stable
        if ((!is.na(sign_changes) && sign_changes >= min_sign_changes_complex) ||
            (!is.na(residual_variance) &&
             residual_variance > complex_variance_threshold)) {

          trend <- "Complex"
          classification_reason <- paste0(
            "Multiple changes (n=", sign_changes,
            ") or high variance (", round(residual_variance, 2), ")"
          )

        } else {
          # Borderline stability: small net change, low variability
          is_within_mono_range <- (total_change < total_change_up) &&
            (total_change > total_change_down)

          max_mono_threshold <- max(abs(total_change_up), abs(total_change_down))

          if (is_within_mono_range &&
              abs(total_change) < max_mono_threshold * 0.5 &&
              !is.na(cv) && cv < stable_cv_threshold * 1.5) {

            trend <- "Borderline_Stable"
            classification_reason <- paste0(
              "Borderline stable, total_change=", round(total_change, 2),
              ", CV=", round(cv, 2)
            )
          } else {
            trend <- "Complex"
            classification_reason <- "Does not fit a clear pattern"
          }
        }
      }
    } else {
      trend <- "Complex"
      classification_reason <- "Insufficient features"
    }

    # ----------------------------------------------------------
    # Phase 7: build return object
    # ----------------------------------------------------------
    details <- if (return_details) {
      list(
        cv = cv,
        relative_change = relative_change,
        autocorr = autocorr,
        fitted_values = pred,
        residuals = residuals,
        residual_variance = residual_variance,
        peak_prominence = peak_prominence,
        valley_prominence = valley_prominence,
        peak_index = max_idx,
        valley_index = min_idx,
        peak_in_middle = peak_in_middle,
        valley_in_middle = valley_in_middle,
        left_slope_norm = left_slope_norm,
        right_slope_norm = right_slope_norm,
        valley_left_norm = valley_left_norm,
        valley_right_norm = valley_right_norm,
        total_change = total_change,
        sign_changes = sign_changes,
        monotone_score = monotone_score,
        classification_reason = classification_reason
      )
    } else NULL

    if (verbose) {
      cat("  → Classified as", trend, "\n")
      cat("    Reason:", classification_reason, "\n")
    }

    return(safe_return(trend, p_value, details))

  }, error = function(e) {
    if (verbose) cat("  ERROR in classification:", e$message, "\n")
    return(safe_return("Insufficient_Data", NA_real_))
  })
}
