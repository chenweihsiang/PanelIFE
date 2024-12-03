#' @title Least Squares Estimation of Linear Panel Data Models with Interactive Fixed Effects
#'
#' @description This function estimates least squares estimator in a linear panel regression model with factors appearing as interactive fixed effects.
#'
#' @details
#' \strong{Disclaimer:} This function is translated and modified from the Matlab function \code{LS_factor.m} by Martin Weidner, and the documentation details are also mainly from the original function.
#' This code is offered with no guarantees. Not all features of this code were properly tested. Please let me know if you find any bugs or encounter any problems while using this code. All feedback is appreciated.
#'
#' \strong{Different computational methods:}
#' \itemize{
#'   \item Method 1 (recommended default method) iterates the following two steps until convergence:
#'     \itemize{
#'       \item Step 1: forgiven \code{beta} compute update for \code{lambda} and \code{f} as principal components of \code{Y - beta * X} (same as in method 1)
#'       \item Step 2: forgiven \code{lambda} and \code{f} run a pooled OLS regression of \code{Y - lambda * t(f)} on \code{X} to update \code{beta}
#'   }
#'   \item Method 2 (described in Bai, 2009) iterates the following two steps until convergence:
#'     \itemize{
#'       \item Step 1: forgiven \code{beta} compute update for \code{lambda} and \code{f} as principal components of \code{Y - beta * X} (same as in method 1)
#'       \item Step 2: forgiven \code{lambda} and \code{f} run a pooled OLS regression of \code{Y - lambda * t(f)} on \code{X} to update \code{beta}
#'   }
#'   The procedure is repeated multiple times with different starting values.
#' }
#'
#' \strong{Comments:}
#' \itemize{
#'   \item Another method would be to use step 1 as in Method 1 & 2, but to replace step 2 with a regression of \code{Y} on either \code{M_lambda * X} or \code{X * M_f}, i.e. to only project out either lambda or f in the step 2 regression.
#'     Bai (2009) mentions this method and refers to Ahn, Lee, and Schmidt (2001), Kiefer (1980) and Sargan (1964) for this.
#'     We have not tested this alternative method, but we suspect that Method 1 performs better in terms of speed of convergence.
#'   \item This alternative method and the method proposed by Bai (2009), i.e. "method 2" here, have the property of reducing the LS objective function in each step.
#'     This is not true for Method 1 and may be seen as a disadvantage of Method 1.
#'     However, we found this to be a nice feature, because we use this property of Method 1 as a stopping rule:
#'     if the LS objective function does not improve, then we know we are "far away" from a proper minimum, so we stop the iteration and begin the iteration with another randomly chosen starting value. Note that multiple runs with different starting values are required anyways for all methods (because the LS objective function may have multiple local minimal).
#'   \item We recommend method 1, because each iteration step is fast and its rate of convergence in our tests was very good (faster than method 2).
#'     However, we have not much explored the relative sensitivity of the different methods towards the choice of starting value. Note that by choosing the quickest method (method 1) one can try out more different starting values of the procedure in the same amount of time.
#'     Nevertheless, it may well be that method 2 or the alternative method described above perform better in certain situations.
#' }
#'
#' @references For a description of the model and the least squares estimator see e.g. Bai (2009, "Panel data models with interactive fixed effects"), or Moon and Weidner (2017, "Dynamic Linear Panel Regression Models with Interactive Fixed Effects"; 2015, "Linear Regression for Panel with Unknown Number of Factors as Interactive Fixed Effects")
#'
#' @param Y \eqn{N \times T} matrix of outcomes, where we assume a balanced panel, i.e. all elements of \code{Y} are known
#' @param X \eqn{K \times N \times T} tensor of regressors, where we assume a balanced panel, i.e. all elements of \code{X} are known
#' @param R A positive integer, indicates the number of interactive fixed effects in the estimation; this \code{R} does not include the number of known factors and loadings
#' @param lambda_known (Optional) \eqn{N \times Rex1} matrix of known factor loadings, e.g., \code{lambda_known = matrix(rep(1, N), nrow = N, ncol = 1)} to control standard time dummies.
#'   Default is set to \code{lambda_known = matrix(NA, nrow = N, ncol = 0)}, i.e., there is no known factor loadings
#' @param f_known (Optional) \eqn{T \times Rex2} matrix of known factor, e.g., \code{f_known = matrix(rep(1, T), nrow = T, ncol = 1)} to control standard individual specific fixed effects.
#'   Default is set to \code{f_known = matrix(NA, nrow = T, ncol = 0)}, i.e., there is no known factor
#' @param report (Optional) Whether or not to report the progress. \code{"silent"} has the program running silently; \code{"report"} has the program reporting what it is doing
#' @param precision_beta (Optional) Defines stopping criteria for numerical optimization, namely optimization is stopped when difference in beta relative to previous opimtization step is smaller than \code{"precision_beta"} (uniformly over all \eqn{K} components of beta).
#'   Note that the actual precision in beta will typically be lower than precision_beta, depending on the convergence rate of the procedure.
#' @param method (Optional) Optimization method option of choice. Options include \code{"m1"} and \code{"m2"}
#' @param start (Optional) \eqn{K \times 1} vector, first starting value for numerical optimization
#' @param repMIN (Optional) Minimal number, which is a positive integer, of runs of optimization with different starting point
#' @param repMAX (Optional) Maximal number, which is a positive integer, of runs of optimization (in case numerical optimization doesn't terminate properly, we do multiple runs even for \code{repMIN = 1})
#' @param M1 (Optional) A positive integer, bandwidth for bias correction for dynamic bias (bcorr1), \code{M1} is the number of lags of correlation between regressors and errors that is corrected for in dynamic bias correction
#' @param M2 (Optional) A non-negative integer, bandwidth for bias correction for time-serial correlation (\code{bcorr3}), \code{M2 = 0} only corrects for time-series heteroscedasticity, while \code{M2 > 0} corrects for time-correlation in errors up to lag \code{M2}
#' @param DoF_adj (Optional) Whether or not to adjust for degree of freedom, where default is set to \code{DoF_adj = FALSE}
#'
#' @return A list of results, where
#' \itemize{
#'   \item \code{beta} is the parameter estimate
#'   \item \code{exitflag = 1} if iteration algorithm properly converged at optimal beta, and \code{exitflag = -1} if iteration algorithm did not properly converge at optimal beta
#'   \item \code{lambda} is the estimate for factor loading
#'   \item \code{f} is the estimate for factors
#'   \item \code{Vbeta1,2,3} are estimated variance-covariance matrices of beta, assuming
#'     \enumerate{
#'       \item homoscedasticity of errors in both dimensions
#'       \item heteroscedasticity of errors in both dimensions
#'       \item allowing for time-serial correlation up to lag M2 (i.e. if M2 == 0, then Vbeta2 == Vbeta3)
#'     }
#'   \item \code{bcorr1,2,3} are estimates for the three different bias components (needs to be subtracted from beta to correct for the bias), where
#'     \enumerate{
#'       \item is bias due to pre-determined regressors
#'       \item is bias due to cross-sectional heteroscedasticity of errors
#'       \item is bias due to time-serial heteroscedasticity and time-serial correlation of errors
#'     }
#'   \item \code{parameter} is the list of input parameters, including \code{Gamma_LS}, \code{alpha}, and \code{clustered_se}
#' }
#'
#' @note
#' We assume that all provided input parameters have values and dimensions as described above.
#' The program could be improved by checking that this is indeed the case.
#'
#' @examples
#' dt <- sample_data(N = 100, T = 20, R = 1, kappa = c(0.5))
#' res <- ls_factor(Y = dt$Y, X = dt$X, R = dt$R, report = "silent",
#'                  precision_beta = 10^(-8), method = "m1",
#'                  start = c(0), repMIN = 3, repMAX = 10, M1 = 2, M2 = 2)
#' sum_res <- summary(res)
#' print(sum_res)
#'
#' @export


ls_factor <- function(Y, X, R,
                      lambda_known = NA, f_known = NA,
                      report = "report", precision_beta = 10^(-8), method = "m1",
                      start, repMIN, repMAX, M1 = 1, M2 = 0, DoF_adj = FALSE) {
  # ========================================
  # Check if inputs are valid
  # ========================================
  if(!("array" %in% class(Y) & length(dim(Y)) == 2)) {
    stop("Format of input 'Y' is incorrect")
  }
  # ===
  if(!("array" %in% class(X) & length(dim(X)) == 3)) {
    stop("Format of input 'X' is incorrect")
  }
  # ===
  if(!(("numeric" %in% class(R) | "integer" %in% class(R)) & length(R) == 1)) {
      stop("Format of input 'R' is incorrect")
  }
  if("numeric" %in% class(R)) {
    if(!(R %% 1 == 0)) {
      stop("Format of input 'R' is incorrect")
    }
  }
  # ===
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
  if(!("character" %in% class(report) & report %in% c("silent", "report"))) {
    stop("Format of input 'report' is incorrect")
  }
  # ===
  if(!("numeric" %in% class(precision_beta) & length(precision_beta) == 1)) {
    stop("Format of input 'precision_beta' is incorrect")
  }
  # ===
  if(!("character" %in% class(method) & method %in% c("m1", "m2"))) {
    stop("Format of input 'method' is incorrect")
  }
  # ===
  if(!(("numeric" %in% class(start) | "array" %in% class(start)) & length(c(start)) == dim(X)[1])) {
    stop("Format of input 'start' is incorrect")
  }
  # ===
  if(!("numeric" %in% class(repMIN) & length(repMIN) == 1)) {
    stop("Format of input 'repMIN' is incorrect")
  }
  if("numeric" %in% class(repMIN)) {
    if(!(repMIN %% 1 == 0 & repMIN > 0)) {
      stop("Format of input 'repMIN' is incorrect")
    }
  }
  # ===
  if(!("numeric" %in% class(repMAX) & length(repMAX) == 1)) {
    stop("Format of input 'repMAX' is incorrect")
  }
  if("numeric" %in% class(repMAX)) {
    if(!(repMAX %% 1 == 0 & repMAX > 0 & repMAX >= repMIN)) {
      stop("Format of input 'repMAX' is incorrect")
    }
  }
  # ===
  if(!("numeric" %in% class(M1) & length(M1) == 1)) {
    stop("Format of input 'M1' is incorrect")
  }
  if("numeric" %in% class(M1)) {
    if(!(M1 %% 1 == 0 & M1 > 0)) {
      stop("Format of input 'M1' is incorrect")
    }
  }
  # ===
  if(!("numeric" %in% class(M2) & length(M2) == 1)) {
    stop("Format of input 'M2' is incorrect")
  }
  if("numeric" %in% class(M2)) {
    if(!(M2 %% 1 == 0 & M2 > 0)) {
      stop("Format of input 'M2' is incorrect")
    }
  }
  # ========================================
  # Setting parameters
  # ========================================
  K <- dim(X)[1]  # number of regressors
  N <- dim(X)[2]  # cross-sectional dimension
  T <- dim(X)[3]  # time-serial dimension
  # ===
  Rex1 <- dim(lambda_known)[2]
  Rex2 <- dim(f_known)[2]
  # ========================================
  # Set up default parameters
  # ========================================
  # Input parameters that are not provided are given default parameters as follows
  if(missing(start)) {
    start <- rep(0, K)
  }
  if(missing(repMIN)) {
    repMIN <- 30
  }
  if(missing(repMAX)) {
    repMAX <- 10 * repMIN
  }
  # ========================================
  # Project out all known factors and loadings
  # ========================================
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
  # ===
  # If N < T, we permute N and T in order to simplify computation of
  #   eigenvectors and eigenvalues (we only solve the eigenvalue
  #   problems for T*T matrices and we make sure that T <= N)
  trans <- 0
  if(N < T) {
    trans <- 1  # Dummy variable to remember that we exchanged N and T dims
    NN <- N
    N <- T
    T <- NN
    Y <- t(Y)
    X <- aperm(X, c(1, 3, 2))
  }
  # ========================================
  # Numerical optimization to obtain beta
  # ========================================
  beta <- matrix(Inf, nrow = length(start), ncol = 1)   # Best beta found so far
  obj0 <- Inf                                           # Objective function at optimal beta found so far
  count <- 0                                            # How many minimization runs properly converged so far
  exitflag <- -1                                        # No proper solution found, yet
  # ===
  for(i in 1:repMAX) {
    if(count < repMIN) {
      # Choose starting value for optimization
      if(i == 1) {
        # First starting value is given by user (or = 0 by default)
        st <- matrix(start, ncol = 1)
      } else {
        # Choose random starting values up from second run
        st <- matrix(start, ncol = 1) + 10 * matrix(stats::rnorm(length(start) * 1), nrow = length(start), ncol = 1)
      }
      # ===
      # Report to user
      if(report == "report") {
        cat(paste0("    Program ls_factor now starting optimization run number: ", i, " / ", repMAX, "\n",
                   "    number of optimization runs that converged so far: ", count, " / ", repMIN, "\n",
                   "    starting value forcurrent run = ", toString(st), "\n"))
      }
      # ===
      # Run actual optimization over beta
      if(method == "m1") {
        list <- minimize_obj_method1(Y, X, R, st, precision_beta)
        para <- list$beta
        obj <- list$obj
        ef <- list$ef
      } else if(method == "m2") {
        list <- minimize_obj_method2(Y, X, R, st, precision_beta)
        para <- list$beta
        obj <- list$obj
        ef <- list$ef
      }
      # ===
      # Report to user
      if(report == "report") {
        if(ef > 0) {
          cat("    Method ", method, " converged at beta = ", toString(para), "\n")
        } else {
          cat("    Method ", method, " did NOT converge. Stopped at beta = ", toString(para), "\n")
        }
        if(obj < obj0) {
          cat("    Final Objective = ", obj, " ==> NEW BEST MINIMUM FOUND\n\n")
        } else {
          cat("    Final Objective = ", obj, " > Best Objective so Far = ", obj0, "\n\n")
        }
      }
      # ===
      # Update estimator, in case better solution found
      if(obj < obj0) {
        obj0 <- obj
        beta <- para # New "global minimum" found
        if(ef > 0) {
          # optimal beta corresponds to point where iteration algorithm properly converged
          exitflag <- 1
        } else {
          exitflag <- -1
        }
      }
      # ===
      # Update counter of "good" solutions
      if(ef > 0) {
        # If method properly converged, then count how many "good" solutions found
        count <- count + 1
      }
    }
  }
  # end of calculation of beta-estimator
  # ========================================
  # Calculate lambda and f for the optimal beta
  # ========================================
  res1 <- Y
  for(k in 1:K) {
    res1 <- res1 - beta[k] * X[k, , ]
  }
  eigen_decomp <- eigen(t(res1) %*% res1)
  V <- eigen_decomp$vectors
  D <- eigen_decomp$values
  Dsort <- sort(eigen_decomp$values)
  Ind <- order(eigen_decomp$values)
  f <- as.matrix(V[ , Ind[(T-R+1):T]]) # Eigenvectors corresponding to largest eigenvalues
  for(r in 1:R) {
    f[ , r] <- f [ ,r] / norm(f[ , r], type = "2")
    if(mean(f[ , r]) < 0) {
      f[ , r] <- -f[ , r]
    }
  }
  lambda <- res1 %*% f
  res <- res1 - lambda %*% t(f) # Estimate for the residuals
  # ===
  # Undo the interchange of N and T now
  if(trans == 1) {
    save <- lambda
    lambda <- f
    f <- save
    res <- t(res)
    NN <- N
    N <- T
    T <- NN
    Y <- t(Y)
    X <- aperm(X, c(1, 3, 2))
  }
  # ========================================
  # Calculate variance-covariance matrix of beta
  # ========================================
  f_all <- cbind(f, f_known)
  lambda_all <- cbind(lambda, lambda_known)
  Pf_all <- f_all %*% solve(t(f_all) %*% f_all) %*% t(f_all)
  Plambda_all <- lambda_all %*% solve(t(lambda_all) %*% lambda_all) %*% t(lambda_all)
  Mf_all <- diag(T) - Pf_all
  Mlambda_all <- diag(N) - Plambda_all
  W <- matrix(0, K, K)
  Omega <- matrix(0, K, K)
  Omega2 <- matrix(0, K, K)
  for (k1 in 1:K) {
    for (k2 in 1:K) {
      Xk1 <- Mlambda_all %*% X[k1, , ] %*% Mf_all
      Xk2 <- Mlambda_all %*% X[k2, , ] %*% Mf_all
      W[k1, k2] <- (1 / (N * T)) * sum(diag(Mf_all %*% t(X[k1, , ]) %*% Mlambda_all %*% X[k2, , ])) # Hessian
      if(DoF_adj == TRUE) {
        Omega[k1, k2] <- (1 / ((N-R-Rex1) * (T-R-Rex2))) * sum((as.vector(Xk1) * as.vector(Xk2)) * (as.vector(res)^2)) # Variance of Score
        Omega2[k1, k2] <- (1 / ((N-R-Rex1) * (T-R-Rex2))) * sum(trunc(t(res * Xk1) %*% (res * Xk2), M2 + 1, M2 + 1))
        # Omega2[k1, k2] <- (1 / ((N-R-Rex1) * (T-R-Rex2))) * sum(diag(trunc(t(res * Xk1) %*% (res * Xk1), M2 + 1, M2 + 1))) # Old and incorrect code from "LS_factor.m"
      } else {
        Omega[k1, k2] <- (1 / (N * T)) * sum((as.vector(Xk1) * as.vector(Xk2)) * (as.vector(res)^2)) # Variance of Score
        Omega2[k1, k2] <- (1 / (N * T)) * sum(trunc(t(res * Xk1) %*% (res * Xk2), M2 + 1, M2 + 1))
        # Omega2[k1, k2] <- (1 / (N * T)) * sum(diag(trunc(t(res * Xk1) %*% (res * Xk1), M2 + 1, M2 + 1))) # Old and incorrect code from "LS_factor.m"
      }
    }
  }
  if(DoF_adj == TRUE) {
    sigma2 <- sum(diag(t(res) %*% res)) / ((N-R-Rex1) * (T-R-Rex2))
  } else {
    sigma2 <- sum(diag(t(res) %*% res)) / (N * T)
  }
  Vbeta1 <- solve(W) * sigma2 / (N * T)
  Vbeta2 <- solve(W) %*% Omega %*% solve(W) / (N * T)
  Vbeta3 <- solve(W) %*% Omega2 %*% solve(W) / (N * T)
  # ========================================
  # Calculate bias estimators
  # ========================================
  B1 <- c()
  B2 <- c()
  B3 <- c()
  for (k in 1:K) {
    XX <- X[k, , ]
    B1[k] <- 1 / sqrt(N*T) * sum(diag(Pf_all %*% trunc(t(res) %*% XX, 0, M1+1)))
    B2[k] <- 1 / sqrt(N*T) * sum(diag(t(XX) %*% Mlambda_all %*% trunc(res %*% t(res), 1, 1) %*% lambda %*% solve(t(lambda) %*% lambda) %*% solve(t(f) %*% f) %*% t(f)))
    B3[k] <- 1 / sqrt(N*T) * sum(diag(trunc(t(res) %*% res, M2+1, M2+1) %*% Mf_all %*% t(XX) %*% lambda %*% solve(t(lambda) %*% lambda) %*% solve(t(f) %*% f) %*% t(f)))
  }
  bcorr1 <- -solve(W) %*% matrix(B1) / sqrt(N*T)
  bcorr2 <- -solve(W) %*% matrix(B2) / sqrt(N*T)
  bcorr3 <- -solve(W) %*% matrix(B3) / sqrt(N*T)
  # ========================================
  # Output results
  # ========================================
  return(structure(list(beta = beta, exitflag = exitflag,
                        lambda = lambda, f = f,
                        Vbeta1 = Vbeta1, Vbeta2 = Vbeta2, Vbeta3 = Vbeta3,
                        bcorr1 = bcorr1, bcorr2 = bcorr2, bcorr3 = bcorr3,
                        parameter = list(lambda_known = lambda_known, f_known = f_known,
                                         report = report, precision_beta = precision_beta,
                                         method = method, start = start,
                                         repMIN = repMIN, repMAX = repMAX,
                                         M1 = M1, M2 = M2)),
                   class = "ls_factor"))
}

