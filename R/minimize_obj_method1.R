#' @title Minimizing Objective Function with Method 1 (Default Method)
#'
#' @description Method 1 iterate the following two steps until convergence:
#' \itemize{
#'   \item Step 1: for given \code{beta} compute update for \code{lambda} and \code{f} as principal components of \code{Y - beta * X}
#'   \item Step 2: for given \code{lambda} and \code{f} update \code{beta} by running a pooled OLS regression of \code{M_lambda * Y * M_f} (or equivalently of just \code{Y} itself) on \code{M_lambda * X * M_f}
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


minimize_obj_method1 <- function(Y, X, R, st, precision_beta) {
  # INPUT: Y  = NxT
  #        X  = KxNxT
  #        st = Kx1 ... starting value for optimization over regression parameter beta
  #        precision_beta = defines stopping criteria, namely optimization
  #                         is stopped when difference in beta after one
  #                         optimization step is smaller than precision_beta
  #                         (uniformly over all K components of beta)

  N <- dim(Y)[1]
  T <- dim(Y)[2]
  K <- dim(X)[1]

  SST <- sum(diag(Y %*% t(Y))) / N / T

  beta <- st  # starting value for beta-minimization
  beta_old <- st + Inf
  obj <- Inf
  diff_obj <- -Inf

  while ((max(abs(beta - beta_old)) > precision_beta) && (diff_obj <= SST * 10^(-10))) {
    # two stopping criteria for iteration of STEP 1 and STEP 2 below:
    # (1) stop if each component of "beta-beta_old" is smaller than "precision_beta"
    #     ==> this is the good case when we have found a proper local minimum
    # (2) stop if diff_obj > SST*10^-8, i.e. if we made no progress in the objective
    #     function during the last iteration
    #     ==> this is the bad case, where the iteration with that particular
    #     starting value is likely not to converge

    # ---- STEP 1: CALCULATE FACTORS AND FACTOR LOADINGS FOR GIVEN beta BY PRINCIPAL COMPONENTS: ----
    res <- get_residuals(Y, X, beta)
    pc <- principal_components(res, R)
    lambda <- pc$lambda
    f <- pc$f
    rm(pc)
    res <- res - lambda %*% t(f)  # residuals after subtracting lambda*f'
    obj_old <- obj  # save old objective function
    obj <- sum(diag(t(res) %*% res)) / N / T  # LS objective function
    diff_obj <- obj - obj_old  # we hopefully make progress in minimizing objective function,
    # i.e. "diff_obj" should better be negative

    if (diff_obj <= 0) {
      # ---- STEP 2: CALCULATE OLS ESTIMATOR FROM REGRESSING M_lambda*Y*M_f (or just Y) on M_lambda*X*M_f: ----
      YY <- as.vector(Y)                                   # flatten Y, i.e. YY now NTx1 vector
      # alternatively, we could define YY as follows, but it should not matter:
      # YY <- Y - lambda %*% pracma::mldivide(lambda, Y)   # project lambda out of Y
      # YY <- t(t(YY) - f %*% pracma::mldivide(f, t(YY)))  # project f out of Y
      # YY <- as.vector(YY)                                # flatten Y, i.e. YY now NTx1 vector
      XX <- matrix(NA, nrow = N*T, ncol = K)               # initialize XX
      for (k in 1:K) {
        xx <- X[k, , ]
        xx <- xx - lambda %*% pracma::mldivide(lambda, xx) # project lambda out of X
        xx <- t(t(xx) - f %*% pracma::mldivide(f, t(xx)))  # project f out of X
        XX[ , k] <- as.vector(xx)                          # flatten X, i.e. XX becomes NTxK matrix
      }
      beta_old <- beta                                     # store old beta
      beta <- pracma::pinv(t(XX) %*% XX) %*% t(XX) %*% YY  # calculate OLS estimator
      # Here, we use the pseudo-inverse, in case t(XX) %*% XX is not invertible
      # to avoide error messages at some points of the optimization.
      # However, at the optimium we should have that beta = XX\YY = pracma::mldivide(XX, YY).
    }
  }

  if (diff_obj <= 0) {
    ef <- 1  # good solution found
  } else {
    ef <- -1 # no good solution found
  }

  obj <- LS_obj(beta, Y, X, R) # calculate objective function for this beta

  return(list(beta = beta, obj = obj, ef = ef))
}

