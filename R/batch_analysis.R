# =============================================================================
# batch_analysis.R
# Multi-subject, parallelised UKF batch pipeline.
#
# Functions:
#   run_subject_ukf()        -- Run all region pairs for a single subject
#   run_all_subjects_ukf()   -- Parallel outer loop over all subjects
#   load_cached_results()    -- Reload incrementally saved subject results
#   compute_group_stats()    -- Mean / SD / CI of parameters per region pair
#   compute_test_retest()    -- ICC between rest1 and rest2 parameter estimates
#   check_identifiability()  -- Perturbation sensitivity around best parameters
# =============================================================================

source("R/constants.R")
source("R/ukf_engine.R")
source("R/optim.R")
source("R/preprocessing.R")

suppressPackageStartupMessages({
  library(parallel)
  library(doParallel)
  library(foreach)
  library(progress)
  library(dplyr)
})


# -----------------------------------------------------------------------------
#' run_subject_ukf
#'
#' Runs iterative_param_optim for every region pair of one subject.
#' Results are saved incrementally to disk after the subject completes.
#'
#' @param smoothed_df  Smoothed data frame for this subject (from smooth_subject_data).
#' @param subj_id      String identifier used for file naming and messages.
#' @param ode_model    ODE model function.
#' @param N_p          Number of model parameters.
#' @param param_guess  Initial parameter vector (length N_p).
#' @param results_dir  Folder to save per-subject RDS files.
#' @param overwrite    Re-run even if a saved file exists (default FALSE).
#' @param ...          Additional arguments passed to iterative_param_optim.
#' @return Data frame: Region_Pair, a, b, k, chisq, steps, param_norm.
#' @export
run_subject_ukf <- function(smoothed_df, subj_id, ode_model, N_p,
                              param_guess  = rep(1, N_p),
                              results_dir  = "results",
                              overwrite    = FALSE,
                              ...) {

  out_file <- file.path(results_dir, paste0(subj_id, ".rds"))

  if (!overwrite && file.exists(out_file)) {
    message("Skipping ", subj_id, " — cached result found: ", out_file)
    return(readRDS(out_file))
  }

  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

  pairs <- build_region_pairs(smoothed_df)
  N_y   <- 2L

  t_vec   <- smoothed_df[[1]]
  dT      <- t_vec[2] - t_vec[1]
  dt      <- UKF_CONSTANTS$DT_FRACTION * dT
  t_dummy <- 0   # FIX: scalar, not a time vector

  rows <- vector("list", length(pairs))

  pb <- progress_bar$new(
    format = paste0("  ", subj_id, " [:bar] :current/:total pairs | ETA: :eta"),
    total  = length(pairs), clear = FALSE, width = 70
  )

  for (idx in seq_along(pairs)) {
    pair      <- pairs[[idx]]
    pair_name <- paste(pair, collapse = "_")
    ukf_data  <- prepare_ukf_input(smoothed_df, pair)

    result <- tryCatch(
      iterative_param_optim(
        param_guess = param_guess,
        t_dummy     = t_dummy,
        ts_data     = ukf_data,
        ode_model   = ode_model,
        N_p         = N_p,
        N_y         = N_y,
        dt          = dt,
        dT          = dT,
        ...
      ),
      error = function(e) {
        message("\n  Error for pair ", pair_name, " in ", subj_id, ": ", conditionMessage(e))
        NULL
      }
    )

    rows[[idx]] <- if (!is.null(result)) {
      data.frame(
        Subject     = subj_id,
        Region_Pair = pair_name,
        a           = result$param_est[1],
        b           = result$param_est[2],
        k           = result$param_est[3],
        chisq       = result$value,
        steps       = result$steps,
        param_norm  = result$param_norm,
        stringsAsFactors = FALSE
      )
    } else {
      data.frame(Subject = subj_id, Region_Pair = pair_name,
                 a = NA, b = NA, k = NA, chisq = NA, steps = NA,
                 param_norm = NA, stringsAsFactors = FALSE)
    }

    pb$tick()
  }

  subject_df <- do.call(rbind, rows)
  saveRDS(subject_df, out_file)
  message("  Saved: ", out_file)
  subject_df
}


# -----------------------------------------------------------------------------
#' run_all_subjects_ukf
#'
#' Parallelised outer loop that calls run_subject_ukf for every subject in
#' data_list.  Uses all available cores minus one.
#'
#' @param data_list    Named list of raw data frames (from load_subject_data).
#' @param ode_model    ODE model function.
#' @param N_p          Number of model parameters.
#' @param param_guess  Initial parameter vector.
#' @param smooth_bw    Kernel smoothing bandwidth (NULL = auto via dpill).
#' @param results_dir  Folder for per-subject RDS files and final CSV.
#' @param n_cores      Number of parallel workers (default: detectCores()-1).
#' @param overwrite    Overwrite existing per-subject results if TRUE.
#' @param ...          Extra args forwarded to iterative_param_optim.
#' @return Data frame of all subjects and all region pairs.
#' @export
run_all_subjects_ukf <- function(data_list, ode_model, N_p,
                                   param_guess  = rep(1, N_p),
                                   smooth_bw    = NULL,
                                   results_dir  = "results",
                                   n_cores      = max(1L, detectCores() - 1L),
                                   overwrite    = FALSE,
                                   ...) {
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

  subj_ids <- names(data_list)
  message("Running UKF for ", length(subj_ids), " subjects on ",
          n_cores, " core(s).\n")

  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  on.exit(stopCluster(cl), add = TRUE)

  all_results <- foreach(
    subj_id  = subj_ids,
    .packages = c("KernSmooth", "pracma", "Matrix", "MASS", "progress"),
    .export   = c("run_subject_ukf", "smooth_subject_data",
                  "build_region_pairs", "prepare_ukf_input",
                  "iterative_param_optim", "UKF_blend", "UKF_dT",
                  "propagate_model", ".stabilise_pd", "UKF_CONSTANTS"),
    .errorhandling = "pass"
  ) %dopar% {
    source("R/constants.R")     # re-source inside worker
    smoothed <- smooth_subject_data(data_list[[subj_id]], bandwidth = smooth_bw)
    run_subject_ukf(smoothed, subj_id, ode_model, N_p,
                    param_guess = param_guess,
                    results_dir = results_dir,
                    overwrite   = overwrite,
                    ...)
  }

  names(all_results) <- subj_ids

  # Collect and save combined CSV
  combined <- do.call(rbind, Filter(is.data.frame, all_results))
  combined_path <- file.path(results_dir, "all_subjects_results.csv")
  write.csv(combined, combined_path, row.names = FALSE)
  message("\nAll results saved to: ", combined_path)

  combined
}


# -----------------------------------------------------------------------------
#' load_cached_results
#'
#' Re-assembles the full results data frame from per-subject RDS files.
#' Useful for adding new subjects without re-running everything.
#'
#' @param results_dir  Folder containing per-subject .rds files.
#' @return Combined data frame.
#' @export
load_cached_results <- function(results_dir = "results") {
  rds_files <- list.files(results_dir, pattern = "\\.rds$", full.names = TRUE)
  if (length(rds_files) == 0)
    stop("No .rds files found in: ", results_dir)

  do.call(rbind, lapply(rds_files, readRDS))
}


# -----------------------------------------------------------------------------
#' compute_group_stats
#'
#' Summarises estimated parameters across subjects for each region pair.
#'
#' @param results_df  Combined results data frame (all subjects).
#' @param params      Parameter column names to summarise (default: a, b, k).
#' @return Grouped summary data frame with mean, SD, 95 % CI, and n.
#' @export
compute_group_stats <- function(results_df,
                                 params = c("a", "b", "k")) {
  stopifnot(all(params %in% colnames(results_df)),
            "Region_Pair" %in% colnames(results_df))

  results_df %>%
    filter(if_all(all_of(params), is.finite)) %>%
    group_by(Region_Pair) %>%
    summarise(
      n = n(),
      across(all_of(params), list(
        mean = ~mean(.x, na.rm = TRUE),
        sd   = ~sd(.x,   na.rm = TRUE),
        ci_lo = ~mean(.x, na.rm = TRUE) - 1.96 * sd(.x, na.rm = TRUE) / sqrt(sum(!is.na(.x))),
        ci_hi = ~mean(.x, na.rm = TRUE) + 1.96 * sd(.x, na.rm = TRUE) / sqrt(sum(!is.na(.x)))
      )),
      .groups = "drop"
    ) %>%
    arrange(desc(k_mean))
}


# -----------------------------------------------------------------------------
#' compute_test_retest
#'
#' Computes intraclass correlation (ICC) between rest1 and rest2 parameter
#' estimates across subjects for each region pair.
#'
#' @param rest1_results  Results data frame for session 1.
#' @param rest2_results  Results data frame for session 2.
#' @param param          Parameter to assess (default "k").
#' @return Data frame: Region_Pair, ICC, lower CI, upper CI, p-value.
#' @export
compute_test_retest <- function(rest1_results, rest2_results,
                                 param = "k") {
  if (!requireNamespace("irr", quietly = TRUE))
    stop("Package 'irr' is required. Install with: install.packages('irr')")

  shared_pairs <- intersect(rest1_results$Region_Pair,
                             rest2_results$Region_Pair)

  icc_rows <- lapply(shared_pairs, function(pair) {
    v1 <- rest1_results[rest1_results$Region_Pair == pair, param]
    v2 <- rest2_results[rest2_results$Region_Pair == pair, param]

    # Align by subject; drop any missing
    common_subj <- intersect(
      rest1_results$Subject[rest1_results$Region_Pair == pair],
      rest2_results$Subject[rest2_results$Region_Pair == pair]
    )
    v1 <- rest1_results[rest1_results$Region_Pair == pair &
                          rest1_results$Subject %in% common_subj, param]
    v2 <- rest2_results[rest2_results$Region_Pair == pair &
                          rest2_results$Subject %in% common_subj, param]

    if (length(v1) < 3) return(NULL)   # need at least 3 subjects for ICC

    icc_res <- irr::icc(cbind(v1, v2), model = "twoway", type = "agreement")
    data.frame(
      Region_Pair = pair,
      ICC         = icc_res$value,
      lower_CI    = icc_res$lbound,
      upper_CI    = icc_res$ubound,
      p_value     = icc_res$p.value,
      n_subjects  = length(v1),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, Filter(Negate(is.null), icc_rows))
}


# -----------------------------------------------------------------------------
#' check_identifiability
#'
#' Perturbs each parameter independently from its best estimate and measures
#' the chi-square sensitivity.  Near-zero sensitivity flags potential
#' non-identifiability for that parameter.
#'
#' @param best_params  Best-fit parameter vector (length N_p).
#' @param ts_data      Time-series matrix.
#' @param ode_model    ODE model function.
#' @param N_p, N_y, dt, dT, R_scale, Q_scale  UKF settings.
#' @param perturbation  Fractional perturbation size (default 0.1 = 10 %).
#' @param param_names  Names for the parameters (default: p1, p2, ...).
#' @return Data frame: parameter, best_value, perturbed_value,
#'         chisq_best, chisq_perturbed, delta_chisq.
#' @export
check_identifiability <- function(best_params, ts_data, ode_model,
                                    N_p, N_y, dt, dT,
                                    R_scale = 0.01, Q_scale = 0.1,
                                    perturbation = 0.1,
                                    param_names  = paste0("p", seq_len(N_p))) {

  t_dummy  <- 0
  base_run <- UKF_blend(t_dummy, ts_data, ode_model,
                         N_p, N_y, best_params, dt, dT,
                         R_scale, Q_scale,
                         forcePositive = TRUE, seeded = TRUE)
  chisq_best <- base_run$chisq

  rows <- lapply(seq_len(N_p), function(i) {
    p_perturbed    <- best_params
    p_perturbed[i] <- best_params[i] * (1 + perturbation)

    perturbed_run  <- UKF_blend(t_dummy, ts_data, ode_model,
                                 N_p, N_y, p_perturbed, dt, dT,
                                 R_scale, Q_scale,
                                 forcePositive = TRUE, seeded = TRUE)

    data.frame(
      parameter       = param_names[i],
      best_value      = best_params[i],
      perturbed_value = p_perturbed[i],
      chisq_best      = chisq_best,
      chisq_perturbed = perturbed_run$chisq,
      delta_chisq     = perturbed_run$chisq - chisq_best,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}
