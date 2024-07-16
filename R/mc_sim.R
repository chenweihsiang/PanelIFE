#' @title Robust Estimation and Inference in Panels with Interactive Fixed Effects
#'
#' @description This method considers estimation and inference for a regression coefficient in panels with interactive fixed effects (i.e., with a factor structure). As the previously developed estimators and confidence intervals (CIs) might be heavily biased and size-distorted when some of the factors are weak, this method has estimators with improved rates of convergence and bias-aware CIs that are uniformly valid regardless of whether the factors are strong or not.
#'
#' @param Y \eqn{N \times T} matrix of outcomes
#' @param X \eqn{K \times N \times T} tensor of regressors
#' @param R A positive integer, indicates the number of interactive fixed effects in the estimation
#' @param Gamma_LS (Optional) A preliminary LS estimate of the matrix of fixed effects; it will be computed if not provided
#' @param alpha (Optional ) It determines the \code{1 - alpha} coverage of the constructed CI, where default is set to \code{alpha = 0.05}
#'
#' @return A list of results, where
#' \itemize{
#'   \item \code{beta} is the point estimate
#' }
#'
#' @noRd


mc_sim <- function(reps = 500, beta0 = 0, R0 = 1, folder_name = "rlogs",
                   arr_N = c(50, 100), arr_T = c(20, 50), arr_kappa = c(0, 0.1, 0.2, 0.5, 1),
                   ncore = NULL, seed = 123) {
  # ========================================
  # Setting
  # ========================================
  i_case <- 0
  res <- list()
  n_cases <- length(arr_N) * length(arr_T) * length(arr_kappa)
  # ========================================
  # Loop over all settings
  # ========================================
  for(i_N in c(1:length(arr_N))) {
    for(i_T in c(1:length(arr_T))) {
      for(i_kappa in c(1:length(arr_kappa))) {
        ticID_this_design <- Sys.time()
        params$N <- arr_N[i_N]
        params$T <- arr_T[i_T]
        params$kappa <- arr_kappa[i_kappa]
        out <- mc_weak_factors(reps = reps, N = N, T = T, R0 = R0, kappa = kappa, beta0 = beta0,
                               folder_name = folder_name, ncore = ncore)
        res[[paste0(i_N, "_", i_T, "_", i_kappa)]] <- list()
        res[[paste0(i_N, "_", i_T, "_", i_kappa)]]$params <- params
        res[[paste0(i_N, "_", i_T, "_", i_kappa)]]$LS.beta <- out$LS$beta
        res[[paste0(i_N, "_", i_T, "_", i_kappa)]]$LS.se <- out$LS$se
        res[[paste0(i_N, "_", i_T, "_", i_kappa)]]$LS.stats <- out$LS$stats
        res[[paste0(i_N, "_", i_T, "_", i_kappa)]]$H.beta <- out$H$beta
        res[[paste0(i_N, "_", i_T, "_", i_kappa)]]$H.bias <- out$H$wc_bias
        res[[paste0(i_N, "_", i_T, "_", i_kappa)]]$H.se <- out$H$se
        res[[paste0(i_N, "_", i_T, "_", i_kappa)]]$H.stats <- out$H$stats
        i_case <- i_case + 1
        cat(paste0("# Design ", i_case, "/", n_cases, " (N,T,kappa)=(", params$N, ",", params$T, ",", params$kappa, ") ",
                   "took ", round(difftime(Sys.time(), ticID_this_design, units = "secs"), 2), "s\n\n"))
      }
    }
  }
  # ========================================
  # Add least favorable critical value (LFCV)
  # ========================================
  LF_cv_LS <- list()
  LF_cv <- list()
  # ===
  for(i_N in c(1:length(arr_N))) {
    for(i_T in c(1:length(arr_T))) {
      cv_LS_tmp <- array(0, dim = c(dim(res[[1]]$LS.stats$cv), length(arr_kappa)))
      cv_tmp <- array(0, dim = c(dim(res[[1]]$H.stats$cv), length(arr_kappa)))
      for(i_kappa in c(1:length(arr_kappa))) {
        res_tmp <- res[[paste0(i_N, "_", i_T, "_", i_kappa)]]
        cv_LS_tmp[ , , , i_kappa] <- res_tmp$LS.stats$cv
        cv_tmp[ , , , i_kappa] <- res_tmp$H.stats$cv
      }
      LF_cv_LS[[paste0(i_N, "_", i_T)]] <- apply(cv_LS_tmp, c(1:3), max)
      LF_cv[[paste0(i_N, "_", i_T)]] <- apply(cv_tmp, c(1:2), max)
    }
  }
  # ===
  LF_cv_LS_arr <- array(0, dim = c(dim(LF_cv_LS[[1]]), length(LF_cv_LS)))
  LF_cv_arr <- array(0, dim = c(dim(LF_cv[[1]]), length(LF_cv)))
  for(i in c(1:length(LF_cv_LS))) {
    LF_cv_LS_arr[ , , , i] <- LF_cv_LS[[i]]
    LF_cv_arr[ , , i] <- LF_cv[[i]]
  }
  LF_cv_LS_max <- apply(LF_cv_LS_arr, c(1:3), max)
  LF_cv_max <- apply(LF_cv_arr, c(1:2), max)
  # ============================================================
  # Output results
  # ============================================================
  saveRDS(list(res = res,
               LF_cv_LS = LF_cv_LS, LF_cv_LS_max = LF_cv_LS_max,
               LF_cv = LF_cv, LF_cv_max = LF_cv_max),
          file = paste0(params$folder_name, "/MC_merged_w_LFCV_",
                        format(Sys.time(), format = "%Y-%m-%d %H-%M-%S"), ".RDS"))
  # ============================================================
  # Summarize results
  # ============================================================
  res <- readRDS(file = list.files("rlogs/R01_for_paper", full.names = TRUE)[grepl("MC_merged_w_LFCV_", list.files("rlogs/R01_for_paper"))])
  # ===
  res_sum <- data.frame(R = 1, N = 100, T = rep(arr_T, each = length(arr_kappa)), kappa = rep(arr_kappa, length(arr_T)),
                        LS_bias = NA, LS_std = NA, LS_rmse = NA, LS_size = NA, LS_length = NA, LS_length_star = NA,
                        Debiased_bias = NA, Debiased_std = NA, Debiased_rmse = NA, Debiased_size = NA, Debiased_length = NA, Debiased_length_star = NA,
                        Unadjusted_size = NA, Unadjusted_length = NA, Unadjusted_length_star = NA)








}

























