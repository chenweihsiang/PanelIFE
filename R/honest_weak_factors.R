#' @title Robust Estimation and Inference in Panels with Interactive Fixed Effects
#'
#' @description This method considers estimation and inference for a regression coefficient in panels with interactive fixed effects (i.e., with a factor structure). As the previously developed estimators and confidence intervals (CIs) might be heavily biased and size-distorted when some of the factors are weak, this method has estimators with improved rates of convergence and bias-aware CIs that are uniformly valid regardless of whether the factors are strong or not.
#'
#' @details
#' \strong{Disclaimer:} This function is the implementation of Armstrong, Weidner, and Zeleneev (2024, "Robust Estimation and Inference in Panels with Interactive Fixed Effects").
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
#'   \item For detailed implementation, please refer to Armstrong, Weidner, and Zeleneev (2024, "Robust Estimation and Inference in Panels with Interactive Fixed Effects").
#' }
#'
#' @references For a description of the model see Armstrong, Weidner, and Zeleneev (2024, "Robust Estimation and Inference in Panels with Interactive Fixed Effects").
#'
#' @param Y \eqn{N \times T} matrix of outcomes
#' @param X \eqn{K \times N \times T} tensor of regressors
#' @param R A positive integer, indicates the number of interactive fixed effects in the estimation. Note that this number does not include the known factors and loadings defined below
#' @param Gamma_LS (Optional) A preliminary LS estimate of the matrix of fixed effects; it will be computed if not provided
#' @param alpha (Optional) It determines the \code{1 - alpha} coverage of the constructed confidence interval, where default is set to \code{alpha = 0.05}
#' @param clustered_se (Optional) Whether or not performing clustered standard error, where default is set to \code{clustered_se = FALSE}
#' @param lambda_known (Optional) \eqn{N \times Rex1} matrix of known factor loadings, e.g., \code{lambda_known = matrix(rep(1, N), nrow = N, ncol = 1)} to control standard time dummies.
#'   Default is set to \code{lambda_known = NA} (or equivalently \code{matrix(NA, nrow = N, ncol = 0)}), i.e., there is no known factor loadings
#' @param f_known (Optional) \eqn{T \times Rex2} matrix of known factor, e.g., \code{f_known = matrix(rep(1, T), nrow = T, ncol = 1)} to control standard individual specific fixed effects.
#'   Default is set to \code{f_known = NA} (or equivalently \code{matrix(NA, nrow = T, ncol = 0)}), i.e., there is no known factor
#' @param itermax (Optional) Maximum iteration allowed while optimizing for the weights A, where default is set to \code{itermax = 50} for faster computation
#' @param reltol (Optional) Relative convergence tolerance for the optimization algorithm to stop if it is unable to reduce the value by \code{reltol * (abs(val) + reltol)} after \code{0.75 * itermax} steps.
#' , where default is set to \code{reltol = 10^(-4)} for faster computation
#'
#' @return A list of results, where
#' \itemize{
#'   \item \code{beta} is the point estimate
#'   \item \code{bias} is the worst-case bias
#'   \item \code{LB} are the lower bounds of the \code{1 - alpha} confidence intervals, from number of weak factors being 0 to R
#'   \item \code{UB} are the upper bounds of the \code{1 - alpha} confidence intervals, from number of weak factors being 0 to R
#'   \item \code{A} is the matrix of weights
#'   \item \code{parameter} is the list of input parameters, including \code{Gamma_LS}, \code{alpha}, and \code{clustered_se}
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


honest_weak_factors <- function(Y, X, R, Gamma_LS = NULL, alpha = 0.05, clustered_se = FALSE,
                                lambda_known = NA, f_known = NA,
                                itermax = 75, reltol = 10^(-6)) {
  # ========================================
  # Setting parameters and project out all known factors and loadings
  # ========================================
  if(!("array" %in% class(lambda_known) | "logical" %in% class(lambda_known))) {
    stop("Format of input 'lambda_known' is incorrect")
  }
  if(!("array" %in% class(f_known) | "logical" %in% class(f_known))) {
    stop("Format of input 'f_known' is incorrect")
  }
  # ===
  if("logical" %in% class(lambda_known)) {
    if(is.na(lambda_known)) {
      lambda_known <- matrix(NA, nrow = nrow(Y), ncol = 0) # Default choice is no known factor loadings
    } else {
      stop("Format of input 'lambda_known' is incorrect")
    }
  }
  if("logical" %in% class(f_known)) {
    if(is.na(f_known)) {
      f_known <- matrix(NA, nrow = ncol(Y), ncol = 0) # Default choice is no known factors
    } else {
      stop("Format of input 'f_known' is incorrect")
    }
  }
  # ===
  if(dim(lambda_known)[1] != nrow(Y)) {
    stop("The dimension of input 'lambda_known' is incorrect")
  }
  if(dim(f_known)[1] != ncol(Y)) {
    stop("The dimension of input 'f_known' is incorrect")
  }
  # ===
  N <- dim(Y)[1]
  T <- dim(Y)[2]
  K <- dim(X)[1]
  b <- 2 * R * (sqrt(N) + sqrt(T))
  # ===
  Rex1 <- dim(lambda_known)[2]
  Rex2 <- dim(f_known)[2]
  # ===
  # Project out (generalized within transformation) all known factor loadings
  if(dim(lambda_known)[2] > 0) {
    Y <- Y - lambda_known %*% pracma::mldivide(lambda_known, Y)
    for(k in 1:K) {
      xx <- X[k, , ]
      X[k, , ] <- xx - lambda_known %*% pracma::mldivide(lambda_known, xx)
    }
  }
  # ===
  # Project out (generalized within transformation) all known factors
  if(dim(f_known)[2] > 0) {
    Y <- t(t(Y) - f_known %*% pracma::mldivide(f_known, t(Y)))
    for(k in 1:K) {
      xx <- X[k, , ]
      X[k, , ] <- t(t(xx) - f_known %*% pracma::mldivide(f_known, t(xx)))
    }
  }
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
    mu_mse <- RcppDE::DEoptim(fn = MSE_crit_fn, lower = min(diag(S)), upper = max(diag(S)),
                              RcppDE::DEoptim.control(trace = FALSE, itermax = itermax, NP = 10 * length(max(diag(S))),
                                                      steptol = round(itermax*0.75, 0), reltol = reltol),
                              X_tmp = X_tmp, ZZ = ZZ, b = b, V = V, S = S, W = W)
    mse[k] <- as.numeric(mu_mse$optim$bestmem)
    A[k, , ] <- compute_mse(X = X_tmp, ZZ = ZZ, mu = mse[k], b = b, V = V, S = S, W = W)$A
  }
  # ========================================
  # Compute preliminary estimator
  # ========================================
  if(is.null(Gamma_LS)) {
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
  for(k in 1:K) {
    delta_pre[k] <- sum(A[k, , ] * Y_tilde_LS)
    Y_diff <- Y_diff - X[k, , ] * delta_pre[k] # Y_diff <- Y - X[k, , ] * delta_pre[k]
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
  # ========================================
  # Loop from number of weak factors Rw = 0 to R
  # ========================================
  beta <- rep(0, K)
  se <- rep(0, K)
  bias <- stats::setNames(cbind(data.frame(0:R), data.frame(matrix(NA, nrow = length(0:R), ncol = K))), c("NumberOfWeakFactors", paste0("Est", 1:K)))
  LB <- stats::setNames(cbind(data.frame(0:R), data.frame(matrix(NA, nrow = length(0:R), ncol = K))), c("NumberOfWeakFactors", paste0("Est", 1:K)))
  UB <- stats::setNames(cbind(data.frame(0:R), data.frame(matrix(NA, nrow = length(0:R), ncol = K))), c("NumberOfWeakFactors", paste0("Est", 1:K)))
  # LB_unadj <- rep(0, K)
  # UB_unadj <- rep(0, K)
  for(Rw in c(0:R)) {
    C <- 2 * Rw * max(svd(U)$d) # Using Rw from this step
    # ===
    for (k in 1:K) {
      A_k <- A[k, , ]
      beta[k] <- sum(A_k * Y_tilde)
      if(clustered_se == TRUE) {
        se[k] <- sqrt(sum((rowSums(A_k * U))^2))
      } else {
        se[k] <- sqrt(sum((A_k * U)^2))
      }
      bias[Rw+1, k+1] <- C * max(svd(A_k)$d)
      LB[Rw+1, k+1] <- beta[k] - bias[Rw+1, k+1] - se[k] * stats::qnorm(1 - alpha / 2)
      UB[Rw+1, k+1] <- beta[k] + bias[Rw+1, k+1] + se[k] * stats::qnorm(1 - alpha / 2)
      # LB_unadj[k] <- beta[k] - se[k] * stats::qnorm(1 - alpha / 2)
      # UB_unadj[k] <- beta[k] + se[k] * stats::qnorm(1 - alpha / 2)
    }
  }
  # ========================================
  # Output results
  # ========================================
  return(structure(list(beta = beta, bias = bias, se = se,
                        LB = LB, UB = UB, A = A,
                        parameter = list(Gamma_LS = Gamma_LS, alpha = alpha, clustered_se = clustered_se)),
                   class = "honest_weak_factors"))

}

