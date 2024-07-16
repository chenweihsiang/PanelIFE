#' @title Minimizing Objective Function with Method 3 (Method Described in Bai, 2009)
#'
#' @description Method 3 iterate the following two steps until convergence:
#' \itemize{
#'   \item Step 1: for given \code{beta} compute update for \code{lambda} and \code{f} as principal components of \code{Y - beta * X} (same as in method 1)
#'   \item Step 2: for given \code{lambda} and \code{f} run a pooled OLS regression of \code{Y - lambda * t(f)} on \code{X} to update \code{beta}
#' }
#' The procedure is repeated multiple times with different starting values
#'
#' @param Y \eqn{N \times T} matrix
#' @param X \eqn{K \times N \times T} tensor
#' @param R A positive integer, indicates the number of interactive fixed effects in the estimation
#' @param st \eqn{K \times 1} vector of starting value for optimization over regression parameter \code{beta}
#' @param precision_beta Defines stopping criteria, namely optimization is stopped when difference in \code{beta} after one optimization step is smaller than \code{precision_beta} (uniformly over all \eqn{K} components of \code{beta})
#'
#' @return A list of results, where
#' \itemize{
#'   \item \code{beta} is the \eqn{Kx1} of optimal beta that was found
#'   \item \code{obj} is the value of LS objective function at optimal beta
#'   \item \code{ef = 1} if procedure properly terminated (according to \code{precision_beta} criteria), and \code{ef = -1} if procedure failed to converge (namely if objective function did not improve during last step, in principle the procedure could still converge afterwards, but if the objective function does not improve, then it seems more promising to stop the optimization and restart with a different starting value)
#' }
#'
#' @noRd


minimize_obj_method3 <- function(Y, X, R, st, precision_beta) {
  # inputs and outputs as in "minimize_obj_method1"

  N <- nrow(Y)
  T <- ncol(Y)
  K <- dim(X)[1]

  XX <- matrix(NA, nrow = N*T, ncol = K)
  for (k in 1:K) {
    xx <- X[k, , ]
    XX[ , k] <- as.vector(xx) # flatten X, i.e. XX becomes NTxK matrix (needed for step 2 below)
  }

  beta <- st                  # starting value for beta-minimization
  beta_old <- st + Inf
  obj <- 0
  diff_obj <- -Inf
  i <- 0
  SST <- sum(diag(Y %*% t(Y))) / N / T
  obj_save <- rep(Inf, 1000)  # we save the objective functions from the previous 1000 iterations

  while ((max(abs(beta - beta_old)) > precision_beta) && (abs(diff_obj) >= 10^(-7)*SST)) {
    # In case the method does not converge (i.e. if |beta-beta_old| does not
    # become sufficiently small) we need a second stopping criterion,
    # which here we choose relatively conservatively,
    # namely, we stop if the objective function did not improve
    # over the last 1000 iterations by at least 10^-7*trace(Y*Y')/N/T.

    # ---- STEP 1: CALCULATE FACTORS AND FACTOR LOADINGS FOR GIVEN beta BY PRINCIPAL COMPONENTS: ----
    res <- get_residuals(Y, X, beta)
    pc <- principal_components(res, R)
    lambda <- pc$lambda
    f <- pc$f
    res <- res - lambda %*% t(f)             # residuals after subtracting lambda*f'
    obj <- sum(diag(t(res) %*% res)) / N / T # LS objective function
    obj_save[(i %% 1000) + 1] <- obj
    diff_obj <- obj - obj_save[((i+1) %% 1000) + 1]
    # Difference between current objective fct. and objective fct. from
    # 1000 iterations ago. In this method (as opposed to method 1) this
    # difference is negative by construction.

    # ---- STEP 2: CALCULATE OLS ESTIMATOR FROM REGRESSING Y-lambda*f' on X: ----
    YY <- Y - lambda %*% t(f)                # residuals after principal components are subtracted from Y
    YY <- as.vector(YY)                      # flatten Y, i.e. YY now NTx1 vector
    beta_old <- beta                         # store old beta
    beta <- pracma::mldivide(XX, YY)         # calculate OLS estimator

    i <- i + 1 # count the number of iterations
  }

  if (max(abs(beta - beta_old)) <= precision_beta) {
    ef <- 1  # good solution found
  } else {
    ef <- -1 # no good solution found
  }

  obj <- LS_obj(beta, Y, X, R) # calculate objective function for this beta

  return(list(beta = beta, obj = obj, ef = ef))
}

