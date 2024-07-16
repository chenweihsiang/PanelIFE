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


mc_output <- function(folder_name = "rlogs", R_0 = 1, N, T) {
  # ============================================================
  # R = 1
  # ============================================================
  arr_N <- c(50, 100, 300)
  arr_T <- c(20, 50, 100, 300)
  arr_kappa <- c(seq(0, 0.25, by = 0.05), 0.50, 1.00)
  # ===
  res <- readRDS(file = list.files("rlogs/R01_for_paper", full.names = TRUE)[grepl("MC_merged_w_LFCV_", list.files("rlogs/R01_for_paper"))])
  # ===
  res_sum <- data.frame(R = 1, N = 100, T = rep(arr_T, each = length(arr_kappa)), kappa = rep(arr_kappa, length(arr_T)),
                        LS_bias = NA, LS_std = NA, LS_rmse = NA, LS_size = NA, LS_length = NA, LS_length_star = NA,
                        Debiased_bias = NA, Debiased_std = NA, Debiased_rmse = NA, Debiased_size = NA, Debiased_length = NA, Debiased_length_star = NA,
                        Unadjusted_size = NA, Unadjusted_length = NA, Unadjusted_length_star = NA)
  # ===
  for(i_N in c(1:length(arr_N))[2]) {
    for(i_T in c(1:length(arr_T))) {
      for(i_kappa in c(1:length(arr_kappa))) {
        res_temp <- res$res[[paste0(i_N, "_", i_T, "_", i_kappa)]]
        R_temp <- res_temp$params$R0
        N_temp <- res_temp$params$N
        T_temp <- res_temp$params$T
        kappa_temp <- res_temp$params$kappa
        res_sum[res_sum$R == R_temp & res_sum$N == N_temp & res_sum$T == T_temp & res_sum$kappa == kappa_temp, -c(1:4)] <-
          c(res_temp$LS.stats$bias[3], res_temp$LS.stats$std[3], res_temp$LS.stats$rmse[3], res_temp$LS.stats$size[2, 3, 2]*100, res_temp$LS.stats$length[2, 3, 2],
            # (alpha = 0.05, bias corrected estimator (using all bias terms), se w/ heteroscedasticity)
            mean(2 * res_temp$LS.se[ , , 2] * res$LF_cv_LS[[paste0(i_N, "_", i_T)]][2, 3, 2]),
            # 2 * (se w/ heteroscedasticity) * (cv with alpha = 0.05, bias corrected estimator (using all bias terms), se w/ heteroscedasticity)
            res_temp$H.stats$bias, res_temp$H.stats$std, res_temp$H.stats$rmse, res_temp$H.stats$size*100, res_temp$H.stats$length,
            mean((res_temp$H.beta[ , 1] + res_temp$H.bias[ , 1] + res_temp$H.se[ , 1] * res$LF_cv[[paste0(i_N, "_", i_T)]][2, 1]) -
                   (res_temp$H.beta[ , 1] - res_temp$H.bias[ , 1] - res_temp$H.se[ , 1] * res$LF_cv[[paste0(i_N, "_", i_T)]][2, 1])),
            res_temp$H.stats$size_unadj * 100,
            mean((res_temp$H.beta[ , 1] + res_temp$H.se[ , 1] * qnorm(1 - 0.05 / 2)) -
                   (res_temp$H.beta[ , 1] - res_temp$H.se[ , 1] * qnorm(1 - 0.05 / 2))),
            mean((res_temp$H.beta[ , 1] + res_temp$H.se[ , 1] * res$LF_cv[[paste0(i_N, "_", i_T)]][2, 1]) -
                   (res_temp$H.beta[ , 1] - res_temp$H.se[ , 1] * res$LF_cv[[paste0(i_N, "_", i_T)]][2, 1])))
        rm(res_temp, R_temp, N_temp, T_temp, kappa_temp)
      }
    }
  }
  rm(i_N, i_T, i_kappa)
  res_sum_1 <- res_sum
  View(sapply(res_sum_1, function(x){round(x, 4)}))


  # ============================================================
  # R = 2
  # ============================================================
  arr_N <- c(100)
  arr_T <- c(20, 50, 100, 300)
  arr_kappa_1 <- c(seq(0, 0.30, by = 0.05), 0.40, 0.50, 1.00)
  arr_kappa_2 <- c(seq(0, 0.30, by = 0.05), 0.40, 0.50, 1.00)
  # ===
  res <- readRDS(file = list.files("rlogs/R02_for_paper", full.names = TRUE)[grepl("MC_merged_w_LFCV_", list.files("rlogs/R02_for_paper"))])
  # ===
  res_sum <- data.frame(R = 2, N = 100, T = 50,
                        kappa_1 = rep(arr_kappa_1, length(c("bias", "std", "rmse", "size", "length", "length_star"))),
                        stats = rep(c("bias", "std", "rmse", "size", "length", "length_star"), each = length(arr_kappa_1)))
  res_sum <- cbind(res_sum, setNames(data.frame(matrix(NA, nrow = nrow(res_sum), ncol = 3*length(arr_kappa_2))),
                                     c(paste0("LS_kappa_2_", arr_kappa_2), paste0("Debiased_kappa_2_", arr_kappa_2), paste0("Unadjusted_kappa_2_", arr_kappa_2))))
  # ===
  for(i_N in c(1:length(arr_N))[1]) {
    for(i_T in c(1:length(arr_T))[2]) {
      for(i_kappa_1 in c(1:length(arr_kappa_1))) {
        for(i_kappa_2 in c(1:i_kappa_1)) {
          res_temp <- res$res[[paste0(i_N, "_", i_T, "_", i_kappa_1, "_", i_kappa_2)]]
          R_temp <- res_temp$params$R0
          N_temp <- res_temp$params$N
          T_temp <- res_temp$params$T
          kappa_1_temp <- res_temp$params$kappa[1]
          kappa_2_temp <- res_temp$params$kappa[2]
          res_sum[res_sum$R == R_temp & res_sum$N == N_temp & res_sum$T == T_temp & res_sum$kappa_1 == kappa_1_temp,
                  paste0("LS_kappa_2_", kappa_2_temp)] <-
            c(res_temp$LS.stats$bias[3], res_temp$LS.stats$std[3], res_temp$LS.stats$rmse[3], res_temp$LS.stats$size[2, 3, 2]*100, res_temp$LS.stats$length[2, 3, 2],
              mean(2 * res_temp$LS.se[ , , 2] * res$LF_cv_LS[[paste0(i_N, "_", i_T)]][2, 3, 2]))
          res_sum[res_sum$R == R_temp & res_sum$N == N_temp & res_sum$T == T_temp & res_sum$kappa_1 == kappa_1_temp,
                  paste0("Debiased_kappa_2_", kappa_2_temp)] <-
            c(res_temp$H.stats$bias, res_temp$H.stats$std, res_temp$H.stats$rmse, res_temp$H.stats$size*100, res_temp$H.stats$length,
              mean((res_temp$H.beta[ , 1] + res_temp$H.bias[ , 1] + res_temp$H.se[ , 1] * res$LF_cv[[paste0(i_N, "_", i_T)]][2, 1]) -
                     (res_temp$H.beta[ , 1] - res_temp$H.bias[ , 1] - res_temp$H.se[ , 1] * res$LF_cv[[paste0(i_N, "_", i_T)]][2, 1])))
          res_sum[res_sum$R == R_temp & res_sum$N == N_temp & res_sum$T == T_temp & res_sum$kappa_1 == kappa_1_temp,
                  paste0("Unadjusted_kappa_2_", kappa_2_temp)] <-
            c(NA, NA, NA, res_temp$H.stats$size_unadj * 100,
              mean((res_temp$H.beta[ , 1] + res_temp$H.se[ , 1] * qnorm(1 - 0.05 / 2)) -
                     (res_temp$H.beta[ , 1] - res_temp$H.se[ , 1] * qnorm(1 - 0.05 / 2))),
              mean((res_temp$H.beta[ , 1] + res_temp$H.se[ , 1] * res$LF_cv[[paste0(i_N, "_", i_T)]][2, 1]) -
                     (res_temp$H.beta[ , 1] - res_temp$H.se[ , 1] * res$LF_cv[[paste0(i_N, "_", i_T)]][2, 1])))
          rm(res_temp, R_temp, N_temp, T_temp, kappa_1_temp, kappa_2_temp)
        }
      }
    }
  }
  rm(i_N, i_T, i_kappa_1, i_kappa_2)
  res_sum_2 <- res_sum
  View(cbind(data.frame(stats = res_sum_2$stats), sapply(res_sum_2[ , -which(colnames(res_sum_2) == "stats")], function(x){round(x, 4)})))


  # ============================================================
  # Output
  # ============================================================
  fwrite(res_sum_1, file = "res_summary_R1.csv")
  # fwrite(sapply(res_sum_1, function(x){round(x, 4)}),
  #        file = "res_summary_R1_round.csv")
  fwrite(res_sum_2, file = "res_summary_R2.csv")
  # fwrite(cbind(data.frame(stats = res_sum_2$stats), sapply(res_sum_2[ , -which(colnames(res_sum_2) == "stats")], function(x){round(x, 4)})),
  #        file = "res_summary_R2_round.csv")

}



