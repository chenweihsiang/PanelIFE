#' @title Robust Estimation and Inference in Panels with Interactive Fixed Effects
#'
#' @description This method considers estimation and inference for a regression coefficient in panels with interactive fixed effects (i.e., with a factor structure). As the previously developed estimators and confidence intervals (CIs) might be heavily biased and size-distorted when some of the factors are weak, this method has estimators with improved rates of convergence and bias-aware CIs that are uniformly valid regardless of whether the factors are strong or not.
#'
#' @details
#' \strong{Disclaimer:} This function is the implementation of Armstrong, Weidner, and Zeleneev (2023, "Robust Estimation and Inference in Panels with Interactive Fixed Effects").
#' This code is offered with no guarantees. Not all features of this code were properly tested. Please let me know if you find any bugs or encounter any problems while using this code. All feedback is appreciated.
#'
#' \strong{Linear panel regression model with weak factors:}
#' \itemize{
#'   \item We consider a linear panel regression model of the form
#'   \deqn{Y_{it} = X_{it} \beta + \sum_{k=1}^{K} Z_{k,it} \delta_{k} + \Gamma_{it} + U_{it}}
#'   where
#'     \itemize{
#'       \item \eqn{Y_{it}}, \eqn{X_{i,t}}, and \eqn{Z_{k,it}} are the observed outcome variable and covariates,
#'       \item \eqn{\Gamma_{it}} is the error component that can be correlated with \eqn{X_{it}} and \eqn{Z_{k,it}},
#'       \item \eqn{U_{it}} is the error component modelled as a mean-zero random shock,
#'       \item and large panel is considered, where both \eqn{N} and \eqn{T} are relatively large.
#'   }
#'   \item Model for \eqn{\Gamma_{it}} is referred to as a factor model with factor loadings \eqn{\lambda_{ir}} and factors \eqn{f_{tr}}, with \eqn{R} being the number of factors.
#'   \item Other than having some strong factors, which requires \eqn{\lambda_{ir}} and \eqn{f_{tr}} to have sufficient variation across \eqn{i} and over \eqn{t}, model here allows weak factors.
#'   \item Proposed method provides the bias-aware confidence intervals that are uniformly valid regardless of whether the factors are strong or not.
#' }
#'
#' \strong{Debiasing approach:}
#' \itemize{
#'   \item Access the preliminary estimate \eqn{\hat{\Gamma}_{pre}} along with a bound \eqn{\hat{C}} on the nuclear norm \eqn{\Vert \Gamma - \hat{\Gamma}_{pre}\Vert_{*}}.
#'   \item Then consider the regression with the augmented outcomes \eqn{\tilde{Y}_{it} := Y_{it} - \hat{\Gamma}_{pre,it}}.
#'   \item Construct the confidence interval with the preliminary estimate \eqn{\hat{\Gamma}_{pre}} and bound \eqn{\hat{C}} on the nuclear norm of its estimation error. Such confidence interval is bias-aware, that is, using the bound \eqn{\hat{C}} will take any remaining bias into account after the previous debias procedure.
#'   \item For detailed implementation, please refer to Armstrong, Weidner, and Zeleneev (2023, "Robust Estimation and Inference in Panels with Interactive Fixed Effects").
#' }
#'
#' @references For a description of the model see Armstrong, Weidner, and Zeleneev (2023, "Robust Estimation and Inference in Panels with Interactive Fixed Effects").
#'
#' @param Y \eqn{N \times T} matrix of outcomes
#' @param X \eqn{K \times N \times T} tensor of regressors
#' @param R A positive integer, indicates the number of interactive fixed effects in the estimation
#' @param Gamma_LS (Optional) A preliminary LS estimate of the matrix of fixed effects; it will be computed if not provided
#' @param alpha (Optional) It determines the \code{1 - alpha} coverage of the constructed confidence interval, where default is set to \code{alpha = 0.05}
#'
#' @return A list of results, where
#' \itemize{
#'   \item \code{beta} is the point estimate
#'   \item \code{bias} is the worst-case bias
#'   \item \code{LB} is the lower bound of the \code{1 - alpha} confidence interval
#'   \item \code{UB} is the upper bound of the \code{1 - alpha} confidence interval
#'   \item \code{LB_unadj} is the lower bound of the \code{1 - alpha} confidence interval without bias correction
#'   \item \code{UB_unadj} is the upper bound of the \code{1 - alpha} confidence interval without bias correction
#'   \item \code{A} is the matrix of weights
#' }
#'
#' @note
#' We assume that all provided input parameters have values and dimensions as described above.
#'
#' @examples
#' dt <- sample_data(N = 100, T = 20, R = 1, kappa = c(0.5))
#' res <- honest_weak_factors(Y = dt$Y, X = dt$X, R = dt$R,
#'                            Gamma_LS = NULL, alpha = 0.05)
#' sum_res <- summary(res)
#' print(sum_res)
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
    LS_factor_res <- ls_factor(Y = Y, X = X, R = R,
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
  return(structure(list(beta = beta, bias = bias, se = se, LB = LB, UB = UB,
                        LB_unadj = LB_unadj, UB_unadj = UB_unadj, A = A),
                   class = "honest_weak_factors"))

}

