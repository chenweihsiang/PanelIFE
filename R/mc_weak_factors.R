#' @title Robust Estimation and Inference in Panels with Interactive Fixed Effects
#'
#' @description This method considers estimation and inference for a regression coefficient in panels with interactive fixed effects (i.e., with a factor structure). As the previously developed estimators and confidence intervals (CIs) might be heavily biased and size-distorted when some of the factors are weak, this method has estimators with improved rates of convergence and bias-aware CIs that are uniformly valid regardless of whether the factors are strong or not.
#'
#' @param reps Number of total repetitions in the simulation
#' @param N Panel dimensions, N (or individual) dimension
#' @param T Panel dimensions, T (or time) dimension
#' @param R0 True number of factors
#' @param kappa Strength of the factors
#' @param beta0 True regressions parameter
#' @param folder_name File path of saving the result
#' @param ncore Number of threads used in the parallel computation, where default is set to "\code{NULL}", i.e., single thread computation
#'
#' @return A list of results, where
#' \itemize{
#'   \item \code{beta} is the point estimate
#'   \item \code{bias} is the worst-case bias
#'   \item \code{LB} is the lower bound of the \code{1 - alpha} CI
#'   \item \code{UB} is the upper bound of the \code{1 - alpha} CI
#'   \item \code{LB_unadj} is the lower bound of the \code{1 - alpha} CI without bias correction
#'   \item \code{UB_unadj} is the upper bound of the \code{1 - alpha} CI without bias correction
#'   \item \code{A} is the matrix of weights
#' }
#'
#' @noRd

# ============================================================
# mc_weak_factors
# ============================================================
mc_weak_factors <- function(reps, N, T, R0, kappa, beta0, folder_name, ncore = NULL) {
  # ========================================
  # File saving path
  # ========================================
  if ("folder_name" %in% names(params)) {
    folder_name <- params$folder_name
    if (!dir.exists(folder_name)) {
      dir.create(file.path(getwd(), folder_name), recursive = TRUE)
    }
  } else {
    folder_name <- "rlogs"
    if (!dir.exists(folder_name)) {
      dir.create(file.path(getwd(), folder_name), recursive = TRUE)
    }
  }
  # ========================================
  # Report to user
  # ========================================
  cat(paste0("\n",
             "# Running MC-WEAK-FACTORS with the following parameters:", "\n",
             "    rep   = ", reps, "\n",
             "    N     = ", N, "\n",
             "    T     = ", T, "\n",
             "    R0    = ", R0, "\n",
             "    kappa = ", paste0(kappa, collapse = " "), "\n",
             "    beta0 = ", beta0, "\n",
             "# Starting MC runs:", "\n"))
  # ========================================
  # Setting parameter
  # ========================================
  # Additional DGP parameters
  sigma_u <- 1
  sigma_x <- 1
  # ===
  # Bai's LS optimization options
  start_beta <- 0
  repMIN <- 3
  repMAX <- 10
  # ===
  R_arr <- R0  # for now we assume that R0 is known
  N_R <- length(R_arr)
  # ===
  # Initialize matrices with zeros
  beta_LS <- matrix(0, nrow = reps, ncol = N_R)      # LS-estimator using R factors in estimation
  beta_LS_BC1 <- matrix(0, nrow = reps, ncol = N_R)  # bias corrected estimator (only using static biases)
  beta_LS_BC2 <- matrix(0, nrow = reps, ncol = N_R)  # bias corrected estimator (using all bias terms)
  beta_LS_std1 <- matrix(0, nrow = reps, ncol = N_R) # estimated standard error, assuming homoscedasticity
  beta_LS_std2 <- matrix(0, nrow = reps, ncol = N_R) # estimated standard error, assuming heteroscedasticity
  # ===
  beta_h <- matrix(0, nrow = reps, ncol = N_R)
  beta_h_bias <- matrix(0, nrow = reps, ncol = N_R)
  beta_h_se <- matrix(0, nrow = reps, ncol = N_R)
  beta_h_LB <- matrix(0, nrow = reps, ncol = N_R)
  beta_h_UB <- matrix(0, nrow = reps, ncol = N_R)
  beta_h_LB_unadj <- matrix(0, nrow = reps, ncol = N_R)
  beta_h_UB_unadj <- matrix(0, nrow = reps, ncol = N_R)
  beta_h_lind <- matrix(0, nrow = reps, ncol = N_R)
  # ===
  lind_fn <- function(A) {
    max(A^2) / sum(A^2)
  }
  # ========================================
  # Simulation loop
  # ========================================
  if(is.null(ncore)) {
    # Using only one thread
    for(rep in c(1:reps)) {
      ticID_iter <- Sys.time()
      set.seed(1091 + 7717 * rep) # for replicability
      # ===
      # Define DGP:
      lambda0 <- matrix(rnorm(N*R0), N, R0) # factor loading
      f0 <- matrix(rnorm(T*R0), T, R0) # factors
      X <- lambda0 %*% t(f0) + sigma_x * matrix(rnorm(N*T), N, T) # regressors
      Y <- beta0 * X + (lambda0 %*% diag(kappa, length(kappa), length(kappa))) %*% t(f0) + sigma_u * matrix(rnorm(N*T), N, T) # outcome variable
      XX <- array(0, dim = c(1, N, T)) # only one regressor
      XX[1, , ] <- X
      # ===
      beta_LS_tmp <- rep(0, N_R)
      beta_LS_BC1_tmp <- rep(0, N_R)
      beta_LS_BC2_tmp <- rep(0, N_R)
      beta_LS_std1_tmp <- rep(0, N_R)
      beta_LS_std2_tmp <- rep(0, N_R)
      # ===
      beta_h_tmp <- rep(0, N_R)
      beta_h_bias_tmp <- rep(0, N_R)
      beta_h_se_tmp <- rep(0, N_R)
      beta_h_LB_tmp <- rep(0, N_R)
      beta_h_UB_tmp <- rep(0, N_R)
      beta_h_LB_unadj_tmp <- rep(0, N_R)
      beta_h_UB_unadj_tmp <- rep(0, N_R)
      beta_h_lind_tmp <- rep(0, N_R)
      # ===
      for (i_R in c(1:N_R)) {
        R <- R_arr[i_R]
        # Bai's LS estimator
        res <- LS_factor(Y, XX, R, "silent", 10^(-8), "m1", start_beta, repMIN, repMAX, 2, 2)
        if (res$exitflag <= 0) {
          cat(paste0("QMLE OPTIMIZATION PROBLEM, may need to increase repMIN or repMAX or optimset-parameters in function qmle\n"))
        }
        beta_LS_tmp[i_R] <- res$beta                                            # LS-estimator using R factors in estimation
        beta_LS_BC1_tmp[i_R] <- res$beta - res$bcorr2 - res$bcorr3              # bias corrected estimator (only using static biases)
        beta_LS_BC2_tmp[i_R] <- res$beta - res$bcorr1 - res$bcorr2 - res$bcorr3 # bias corrected estimator (using all bias terms)
        beta_LS_std1_tmp[i_R] <- sqrt(res$Vbeta1)                               # estimated standard error, assuming homoscedasticity
        beta_LS_std2_tmp[i_R] <- sqrt(res$Vbeta2)                               # estimated standard error, assuming heteroscedasticity
        # ===
        # Honest weak factors
        res <- NULL
        while(is.null(res) | "try-error" %in% class(res)) {
          # Different operating system may have different BLAS/LAPACK setting, and some of them will lead to rare failure while doing SVD
          # Using try in the while loop will add additional layer of randomness (which still controlled by the seed) to get around this issue
          res <- try(honest_weak_factors(Y, XX, R), silent = TRUE)
        }
        beta_h_tmp[i_R] <- res$beta
        beta_h_bias_tmp[i_R] <- res$bias
        beta_h_se_tmp[i_R] <- res$se
        beta_h_LB_tmp[i_R] <- res$LB
        beta_h_UB_tmp[i_R] <- res$UB
        beta_h_LB_unadj_tmp[i_R] <- res$LB_unadj
        beta_h_UB_unadj_tmp[i_R] <- res$UB_unadj
        beta_h_lind_tmp[i_R] <- lind_fn(res$A)
      }
      # ===
      beta_LS[rep, ] <- beta_LS_tmp
      beta_LS_BC1[rep, ] <- beta_LS_BC1_tmp
      beta_LS_BC2[rep, ] <- beta_LS_BC2_tmp
      beta_LS_std1[rep, ] <- beta_LS_std1_tmp
      beta_LS_std2[rep, ] <- beta_LS_std2_tmp
      # ===
      beta_h[rep, ] <- beta_h_tmp
      beta_h_bias[rep, ] <- beta_h_bias_tmp
      beta_h_se[rep, ] <- beta_h_se_tmp
      beta_h_LB[rep, ] <- beta_h_LB_tmp
      beta_h_UB[rep, ] <- beta_h_UB_tmp
      beta_h_LB_unadj[rep, ] <- beta_h_LB_unadj_tmp
      beta_h_UB_unadj[rep, ] <- beta_h_UB_unadj_tmp
      beta_h_lind[rep, ] <- beta_h_lind_tmp
      # ===
      # Report to user
      if (rep %% 100 == 0 || rep == 1) {
        cat(paste0("    rep ", rep, "/", reps, ", took ", round(difftime(Sys.time(), ticID_iter, units = "secs"), 2), "s\n"))
      }
    }
  } else {
    # Parallel computation
    rep_func <- function(rep) {
      ticID_iter <- Sys.time()
      set.seed(1091 + 7717 * rep) # for replicability
      # ===
      # Define DGP:
      lambda0 <- matrix(rnorm(N*R0), N, R0) # factor loading
      f0 <- matrix(rnorm(T*R0), T, R0) # factors
      X <- lambda0 %*% t(f0) + sigma_x * matrix(rnorm(N*T), N, T) # regressors
      Y <- beta0 * X + (lambda0 %*% diag(kappa, length(kappa), length(kappa))) %*% t(f0) + sigma_u * matrix(rnorm(N*T), N, T) # outcome variable
      XX <- array(0, dim = c(1, N, T)) # only one regressor
      XX[1, , ] <- X
      # ===
      beta_LS_tmp <- rep(0, N_R)
      beta_LS_BC1_tmp <- rep(0, N_R)
      beta_LS_BC2_tmp <- rep(0, N_R)
      beta_LS_std1_tmp <- rep(0, N_R)
      beta_LS_std2_tmp <- rep(0, N_R)
      # ===
      beta_h_tmp <- rep(0, N_R)
      beta_h_bias_tmp <- rep(0, N_R)
      beta_h_se_tmp <- rep(0, N_R)
      beta_h_LB_tmp <- rep(0, N_R)
      beta_h_UB_tmp <- rep(0, N_R)
      beta_h_LB_unadj_tmp <- rep(0, N_R)
      beta_h_UB_unadj_tmp <- rep(0, N_R)
      beta_h_lind_tmp <- rep(0, N_R)
      # ===
      for (i_R in c(1:N_R)) {
        R <- R_arr[i_R]
        # Bai's LS estimator
        res <- LS_factor(Y, XX, R, "silent", 10^(-8), "m1", start_beta, repMIN, repMAX, 2, 2)
        if (res$exitflag <= 0) {
          cat(paste0("QMLE OPTIMIZATION PROBLEM, may need to increase repMIN or repMAX or optimset-parameters in function qmle\n"))
        }
        beta_LS_tmp[i_R] <- res$beta                                            # LS-estimator using R factors in estimation
        beta_LS_BC1_tmp[i_R] <- res$beta - res$bcorr2 - res$bcorr3              # bias corrected estimator (only using static biases)
        beta_LS_BC2_tmp[i_R] <- res$beta - res$bcorr1 - res$bcorr2 - res$bcorr3 # bias corrected estimator (using all bias terms)
        beta_LS_std1_tmp[i_R] <- sqrt(res$Vbeta1)                               # estimated standard error, assuming homoscedasticity
        beta_LS_std2_tmp[i_R] <- sqrt(res$Vbeta2)                               # estimated standard error, assuming heteroscedasticity
        # ===
        # Honest weak factors
        res <- NULL
        while(is.null(res) | "try-error" %in% class(res)) {
          # Different operating system may have different BLAS/LAPACK setting, and some of them will lead to rare failure while doing SVD
          # Using try in the while loop will add additional layer of randomness (which still controlled by the seed) to get around this issue
          res <- try(honest_weak_factors(Y, XX, R), silent = TRUE)
        }
        # ===
        beta_h_tmp[i_R] <- res$beta
        beta_h_bias_tmp[i_R] <- res$bias
        beta_h_se_tmp[i_R] <- res$se
        beta_h_LB_tmp[i_R] <- res$LB
        beta_h_UB_tmp[i_R] <- res$UB
        beta_h_LB_unadj_tmp[i_R] <- res$LB_unadj
        beta_h_UB_unadj_tmp[i_R] <- res$UB_unadj
        beta_h_lind_tmp[i_R] <- lind_fn(res$A)
      }
      return(list(rep = rep,
                  beta_LS_tmp = beta_LS_tmp,
                  beta_LS_BC1_tmp = beta_LS_BC1_tmp, beta_LS_BC2_tmp = beta_LS_BC2_tmp,
                  beta_LS_std1_tmp = beta_LS_std1_tmp, beta_LS_std2_tmp = beta_LS_std2_tmp,
                  beta_h_tmp = beta_h_tmp, beta_h_bias_tmp = beta_h_bias_tmp, beta_h_se_tmp = beta_h_se_tmp,
                  beta_h_LB_tmp = beta_h_LB_tmp, beta_h_UB_tmp = beta_h_UB_tmp,
                  beta_h_LB_unadj_tmp = beta_h_LB_unadj_tmp, beta_h_UB_unadj_tmp = beta_h_UB_unadj_tmp,
                  beta_h_lind_tmp = beta_h_lind_tmp))
    }
    cl <- makeCluster(ncore)
    registerDoParallel(cl) # Setup parallel back-end to use many processors and not to overload
    rep_res <- foreach(rep = c(1:reps), .packages = c("abind", "PanelIFE"), .export = c("LS_factor", "honest_weak_factors")) %dopar% {
      return(rep_func(rep))
    }
    stopCluster(cl) # Stop cluster
    registerDoSEQ()
    # ===
    # Summarize results
    for (rep in c(1:reps)) {
      beta_LS[rep, ] <- rep_res[[rep]]$beta_LS_tmp
      beta_LS_BC1[rep, ] <- rep_res[[rep]]$beta_LS_BC1_tmp
      beta_LS_BC2[rep, ] <- rep_res[[rep]]$beta_LS_BC2_tmp
      beta_LS_std1[rep, ] <- rep_res[[rep]]$beta_LS_std1_tmp
      beta_LS_std2[rep, ] <- rep_res[[rep]]$beta_LS_std2_tmp
      # ===
      beta_h[rep, ] <- rep_res[[rep]]$beta_h_tmp
      beta_h_bias[rep, ] <- rep_res[[rep]]$beta_h_bias_tmp
      beta_h_se[rep, ] <- rep_res[[rep]]$beta_h_se_tmp
      beta_h_LB[rep, ] <- rep_res[[rep]]$beta_h_LB_tmp
      beta_h_UB[rep, ] <- rep_res[[rep]]$beta_h_UB_tmp
      beta_h_LB_unadj[rep, ] <- rep_res[[rep]]$beta_h_LB_unadj_tmp
      beta_h_UB_unadj[rep, ] <- rep_res[[rep]]$beta_h_UB_unadj_tmp
      beta_h_lind[rep, ] <- rep_res[[rep]]$beta_h_lind_tmp
    }
  }
  # ========================================
  # Prepare for outputting results
  # ========================================
  out <- list()
  out$params <- params
  # ===
  out$LS <- list()
  out$LS$beta <- cbind(beta_LS, beta_LS_BC1, beta_LS_BC2)
  out$LS$se <- abind::abind(beta_LS_std1, beta_LS_std2, along = 3)
  out$LS$stats <- gen_stats(out$LS$beta, beta0, out$LS$se)
  # ===
  out$H <- list()
  out$H$beta <- beta_h
  out$H$wc_bias <- beta_h_bias
  out$H$se <- beta_h_se
  out$H$LB <- beta_h_LB
  out$H$UB <- beta_h_UB
  out$H$lind <- beta_h_lind
  out$H$stats <- gen_stats(out$H$beta, beta0, array(out$H$se, dim = c(dim(out$H$se), 1)))
  out$H$lind_stats <- gen_stats(out$H$lind)
  out$H$stats$length <- mean(out$H$UB - out$H$LB)
  out$H$stats$size <- mean((beta0 < out$H$LB) | (beta0 > out$H$UB))
  out$H$LB_unadj <- beta_h_LB_unadj
  out$H$UB_unadj <- beta_h_UB_unadj
  out$H$stats$size_unadj <- mean((beta0 < out$H$LB_unadj) | (beta0 > out$H$UB_unadj))
  # ===
  if ("outputfile" %in% names(params)) {
    outputfile <- params$outputfile
  } else {
    outputfile <- sprintf("%s/mc_N%i_T%i_R0%i_kappa%s_reps%i_%s", folder_name, N, T, R0, paste0(sprintf("%03d", round(kappa*100)), collapse = "_"), reps, format(Sys.time(), format = "%Y-%m-%d %H-%M-%S"))
  }
  # ========================================
  # Output results
  # ========================================
  cat(paste0("# Save result to file: ", outputfile, "\n"))
  saveRDS(out, file = paste0(outputfile, ".RDS"))
  return(out)
}

