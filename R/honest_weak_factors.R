#' @title Robust Estimation and Inference in Panels with Interactive Fixed Effects
#'
#' @description This method considers estimation and inference for a regression coefficient in panels with interactive fixed effects (i.e., with a factor structure). As the previously developed estimators and confidence intervals (CIs) might be heavily biased and size-distorted when some of the factors are weak, this method has estimators with improved rates of convergence and bias-aware CIs that are uniformly valid regardless of whether the factors are strong or not.
#'
#' @param Y \eqn{N \times T} matrix of outcomes
#' @param X \eqn{K \times N \times T} tensor of regressors
#' @param R A positive integer, indicates the number of interactive fixed effects in the estimation
#' @param Gamma_LS (Optional) A preliminary LS estimate of the matrix of fixed effects; it will be computed if not provided
#' @param alpha (Optional) It determines the \code{1 - alpha} coverage of the constructed CI, where default is set to \code{alpha = 0.05}
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
#' @export


honest_weak_factors <- function(Y, X, R, Gamma_LS = NULL, alpha = 0.05) {
  # ========================================
  # Setting parameter
  # ========================================
  N <- dim(Y)[1]
  T <- dim(Y)[2]
  K <- dim(X)[1]
  b <- 2 * R * (sqrt(N) + sqrt(T))
  # ========================================
  # Compute the optimal weights
  # ========================================
  A <- array(0, dim = c(K, N, T))
  mu_mse <- rep(0, K)
  mse <- rep(0, K)
  for (k in 1:K) {
    X_tmp <- X[k, , ]
    svd_res <- svd(X_tmp)
    V <- svd_res$u
    S <- diag(svd_res$d)
    W <- svd_res$v
    if (K == 1) {
      ZZ <- NULL
    } else {
      ZZ <- make_ZZ(array(X[-k, , ], dim = c(dim(X)[1]-1, dim(X)[2:3])))
    }
    # ===
    MSE_crit_fn <- function(mu, X_tmp, ZZ, b, V, S, W) {
      return(compute_mse(X_tmp, ZZ, mu, b, V, S, W)$mse)
    }
    # ===
    # Nonlinear global optimization through "RcppDE" package
    mu_mse <- RcppDE::DEoptim(fn = MSE_crit_fn, lower = 10^(-4), upper = max(diag(S)),
                              RcppDE::DEoptim.control(trace = FALSE, itermax = 50, NP = 10,
                                                      steptol = 20, reltol = 10^(-4)),
                              X_tmp = X_tmp, ZZ = ZZ, b = b, V = V, S = S, W = W)
    mse[k] <- as.numeric(mu_mse$optim$bestmem)
    A[k, , ] <- compute_mse(X = X_tmp, ZZ = ZZ, mu = mse[k], b = b, V = V, S = S, W = W)$A
  }
  # ========================================
  # Compute preliminary estimator
  # ========================================
  if (is.null(Gamma_LS)) {
    start_delta <- matrix(0, nrow = dim(X)[1], ncol = 1)
    LS_factor_res <- LS_factor(Y = Y, X = X, R = R,
                               report = "silent", precision_beta = 10^(-8), method = "m1",
                               start = start_delta, repMIN = 3, repMAX = 10, M1 = 2, M2 = 2)
    delta_LS <- LS_factor_res$beta
    lambda <- LS_factor_res$lambda
    f <- LS_factor_res$f
    Gamma_LS <- lambda %*% t(f)
  }
  Y_tilde_LS <- Y - Gamma_LS
  delta_pre <- rep(0, K)
  Y_diff <- Y
  for (k in 1:K) {
    delta_pre[k] <- sum(A[k, , ] * Y_tilde_LS)
    Y_diff <- Y - X[k, , ] * delta_pre[k]
  }
  # ========================================
  # Compute the final estimate
  # ========================================
  svd_res <- svd(Y_diff)
  svd_res$u <- matrix(svd_res$u[ , 1:R], nrow = nrow(svd_res$u), ncol = R)
  svd_res$v <- matrix(svd_res$v[ , 1:R], nrow = nrow(svd_res$v), ncol = R)
  V <- svd_res$u
  S <- matrix(diag(svd_res$d)[1:R, 1:R], nrow = R, ncol = R)
  W <- svd_res$v
  Gamma_pre <- matrix(V[ , 1:R], ncol = R) %*% matrix(S[1:R, 1:R], nrow = R, ncol = R) %*% t(matrix(W[ , 1:R], ncol = R))
  Y_tilde <- Y - Gamma_pre
  U <- Y_diff - Gamma_pre
  C <- 2 * R * max(svd(U)$d)
  # ===
  beta <- rep(0, K)
  bias <- rep(0, K)
  se <- rep(0, K)
  LB <- rep(0, K)
  UB <- rep(0, K)
  LB_unadj <- rep(0, K)
  UB_unadj <- rep(0, K)
  # ===
  for (k in 1:K) {
    A_k <- A[k, , ]
    beta[k] <- sum(A_k * Y_tilde)
    bias[k] <- C * max(svd(A_k)$d)
    se[k] <- sqrt(sum((A_k * U)^2))
    LB[k] <- beta[k] - bias[k] - se[k] * stats::qnorm(1 - alpha / 2)
    UB[k] <- beta[k] + bias[k] + se[k] * stats::qnorm(1 - alpha / 2)
    LB_unadj[k] <- beta[k] - se[k] * stats::qnorm(1 - alpha / 2)
    UB_unadj[k] <- beta[k] + se[k] * stats::qnorm(1 - alpha / 2)
  }
  # ========================================
  # Output results
  # ========================================
  return(list(beta = beta, bias = bias, se = se, LB = LB, UB = UB, LB_unadj = LB_unadj, UB_unadj = UB_unadj, A = A))
}

