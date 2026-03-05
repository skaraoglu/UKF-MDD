# =============================================================================
# plotting.R
# Unified ggplot2-based plotting functions for the MDD UKF pipeline.
# All base-R plotting variants have been removed; only ggplot versions remain.
#
# Functions:
#   plot_signal_and_smooth()     -- Raw + kernel-smoothed ROI signals
#   plot_ukf_fit()               -- UKF estimated trajectory vs. smoothed data
#   plot_chisq_over_iterations() -- Chi-square convergence trace
#   plot_parameter_distributions() -- Histograms of a, b, k across region pairs
#   plot_parameter_heatmap()     -- Region-pair coupling matrix
#
# Bugs fixed vs. original:
#   - Row indices for observed states computed from N_p, not hardcoded to 3/4
#   - Chi-square loss x-axis labelled "Iteration", not "Time Step"
#   - library() calls removed from inside functions (must be loaded at top)
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(tidyr)
  library(dplyr)
})


# -----------------------------------------------------------------------------
#' plot_signal_and_smooth
#'
#' Overlays raw and kernel-smoothed signals for every ROI in a subject.
#'
#' @param raw_df      Raw data frame (first col = time, rest = ROI signals).
#' @param smoothed_df Smoothed data frame (same structure).
#' @param subject_id  String used in the plot title.
#' @return Invisibly returns a list of ggplot objects (one per ROI).
#' @export
plot_signal_and_smooth <- function(raw_df, smoothed_df,
                                    subject_id = "") {
  region_names <- colnames(smoothed_df)[-1]
  time_col     <- colnames(smoothed_df)[1]
  plots        <- list()

  for (rgn in region_names) {
    df <- data.frame(
      Time     = raw_df[[time_col]],
      Raw      = raw_df[[rgn]],
      Smoothed = smoothed_df[[rgn]]
    )

    p <- ggplot(df, aes(x = Time)) +
      geom_point(aes(y = Raw,      colour = "Raw"),      shape = 3, size = 0.8) +
      geom_line (aes(y = Smoothed, colour = "Smoothed"), linewidth = 0.8) +
      scale_colour_manual(name   = NULL,
                          values = c(Raw = "grey50", Smoothed = "firebrick")) +
      labs(title    = if (nchar(subject_id) > 0)
                        paste0(subject_id, " — ", rgn)
                      else rgn,
           x = "Time (TR)", y = "BOLD signal") +
      theme_minimal(base_size = 11) +
      theme(legend.position = "top")

    plots[[rgn]] <- p
    print(p)
  }

  invisible(plots)
}


# -----------------------------------------------------------------------------
#' plot_ukf_fit
#'
#' Plots the UKF estimated state trajectory against the smoothed input data
#' for a single region pair.
#'
#' @param ukf_out      Return value of UKF_blend or iterative_param_optim.
#' @param smoothed_data  Matrix (T x 3): [time, signal_A, signal_B].
#' @param region_pair  Character vector c("ROI_A", "ROI_B") for axis labels.
#' @param top_title    String title for the first panel.
#' @return Invisibly returns a list of two ggplot objects.
#' @export
plot_ukf_fit <- function(ukf_out, smoothed_data,
                          region_pair = c("Region A", "Region B"),
                          top_title   = "UKF Fit") {
  # FIX: infer N_p from param_est length; observed states start at row N_p+1
  N_p    <- length(as.numeric(ukf_out$param_est))
  T_len  <- ncol(ukf_out$xhat)
  time   <- seq_len(T_len)

  param_label <- paste0(
    paste0("p", seq_len(N_p), "=",
           round(as.numeric(ukf_out$param_est), 3),
           collapse = ",  ")
  )

  make_panel <- function(obs_row, smooth_col, ylabel) {
    df <- data.frame(
      Time     = time,
      UKF      = ukf_out$xhat[obs_row, ],
      Observed = smoothed_data[, smooth_col]
    )
    ggplot(df, aes(x = Time)) +
      geom_point(aes(y = UKF,      colour = "UKF"),      size = 0.7) +
      geom_line (aes(y = Observed, colour = "Observed"), linewidth = 0.8) +
      scale_colour_manual(name   = NULL,
                          values = c(UKF = "steelblue", Observed = "firebrick")) +
      labs(y = ylabel) +
      theme_minimal(base_size = 11) +
      theme(legend.position = "top", axis.title.x = element_blank())
  }

  p1 <- make_panel(N_p + 1, 2, region_pair[1]) +
    ggtitle(top_title)

  p2 <- make_panel(N_p + 2, 3, region_pair[2]) +
    ggtitle(param_label) +
    labs(x = "Time (TR)")

  print(p1)
  print(p2)
  invisible(list(p1 = p1, p2 = p2))
}


# -----------------------------------------------------------------------------
#' plot_chisq_over_iterations
#'
#' Plots chi-square over UKF iterations (from iterative_param_optim output).
#'
#' @param iter_opt  Return value of iterative_param_optim.
#' @param title     Plot title.
#' @return ggplot object.
#' @export
plot_chisq_over_iterations <- function(iter_opt,
                                        title = "Chi-Square Over Iterations") {
  if (is.null(iter_opt$chisq) || length(iter_opt$chisq) == 0)
    stop("iter_opt$chisq is missing or empty.")

  # FIX: x-axis labelled "Iteration", not "Time Step"
  df <- data.frame(
    Iteration   = seq_along(iter_opt$chisq),
    ChiSquare   = iter_opt$chisq
  )

  p <- ggplot(df, aes(x = Iteration, y = ChiSquare)) +
    geom_line(colour = "#2C3E50", linewidth = 0.9) +
    scale_y_log10() +
    labs(title = title, x = "Iteration", y = expression(chi^2)) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5))

  print(p)
  invisible(p)
}


# -----------------------------------------------------------------------------
#' plot_parameter_distributions
#'
#' Histograms of estimated a, b, k and chi-square across all region pairs.
#'
#' @param results_df  Data frame with columns: Region_Pair, a, b, k, chisq.
#' @return Invisibly returns named list of ggplot objects.
#' @export
plot_parameter_distributions <- function(results_df) {
  params <- c("a", "b", "k", "chisq")
  colours <- c(a = "#3498DB", b = "#E74C3C", k = "#2ECC71", chisq = "#9B59B6")
  plots  <- list()

  for (param in params) {
    if (!param %in% colnames(results_df)) {
      message("Column '", param, "' not found in results_df — skipping.")
      next
    }
    vals <- results_df[[param]]
    vals <- vals[is.finite(vals) & vals > 0]   # log10 requires positive values

    df <- data.frame(value = vals)
    p  <- ggplot(df, aes(x = value)) +
      geom_histogram(fill = colours[[param]], colour = "black", alpha = 0.75,
                     bins = 40) +
      scale_x_log10(labels = scales::label_comma()) +
      scale_y_log10() +
      labs(title = paste("Distribution of parameter", param),
           x     = param, y = "Count (log10)") +
      theme_minimal(base_size = 11)

    plots[[param]] <- p
    print(p)
  }
  invisible(plots)
}


# -----------------------------------------------------------------------------
#' plot_parameter_heatmap
#'
#' Visualizes one parameter (e.g. coupling k) as a symmetric region × region
#' heatmap.  Region pairs not present are shown as NA.
#'
#' @param results_df   Data frame with columns Region_Pair, and the param.
#' @param param        Column name to display (default "k").
#' @param title        Plot title.
#' @return ggplot object.
#' @export
plot_parameter_heatmap <- function(results_df,
                                    param = "k",
                                    title = paste("Coupling parameter:", param)) {
  stopifnot(all(c("Region_Pair", param) %in% colnames(results_df)))

  # Split "ROI_A_ROI_B" back into two columns
  parts <- strsplit(results_df$Region_Pair, "_(?=[^_]+$)", perl = TRUE)
  df    <- data.frame(
    Region_A = sapply(parts, `[[`, 1),
    Region_B = sapply(parts, `[[`, 2),
    value    = results_df[[param]],
    stringsAsFactors = FALSE
  )
  # Mirror for symmetric display
  df_mirror <- data.frame(Region_A = df$Region_B,
                           Region_B = df$Region_A,
                           value    = df$value)
  df_full   <- rbind(df, df_mirror)

  p <- ggplot(df_full, aes(x = Region_A, y = Region_B, fill = value)) +
    geom_tile(colour = "white") +
    scale_fill_viridis_c(option = "plasma", na.value = "grey90") +
    labs(title = title, x = NULL, y = NULL, fill = param) +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title  = element_text(hjust = 0.5))

  print(p)
  invisible(p)
}
