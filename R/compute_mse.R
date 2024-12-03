#' @title Compute Weights and MSE
#'
#' @description This function compute the weights \code{A} and the mean squared error defined by definition 2.2 in the paper. The computation follows appendix B: nuclear norm regularized "partialling out" regression
#'
#' @param X \eqn{N \times T} matrix
#' @param ZZ NULL (\code{K == 1} i.e. no covariate) or \eqn{(NT) \times (K-1)} matrix
#' @param mu Scalar, the penalty on the nuclear norm
#' @param b Scalar, the tuning parameter
#' @param V V in the singular value decomposition of SVD(X) = VSW'
#' @param S S in the singular value decomposition of SVD(X) = VSW'
#' @param W W in the singular value decomposition of SVD(X) = VSW'
#'
#' @return A list of results, where
#' \itemize{
#'   \item \code{A} are weights defined by def 2.2 in the paper
#'   \item \code{mse} is mean squared error defined by def 2.2 in the paper
#' }
#'
#' @noRd


compute_mse <- function(X, ZZ, mu, b, V, S, W) {
  if(!is.null(ZZ)) {
    gamma_res <- compute_gamma(X, ZZ, X)
    gamma <- gamma_res$gamma
    X_r_gamma <- gamma_res$X_r_gamma
    tol <- 10^(-4)
    delta_gamma <- 1
    while(delta_gamma > tol) {
      X_r_pi <- compute_X_r_pi(X_r_gamma, mu, X)
      gamma_res <- compute_gamma(X_r_pi, ZZ, X)
      gamma_u <- gamma_res$gamma
      X_r_gamma <- gamma_res$X_r_gamma
      delta_gamma <- max(abs(gamma - gamma_u))
      gamma <- gamma_u
    }
    Omega <- X_r_pi - t(array(ZZ %*% gamma, dim = c(dim(X)[2], dim(X)[1])))
  } else {
    Omega <- V %*% diag(pmin(mu, diag(S))) %*% t(W)
  }
  # Following definition 2.2 in the paper
  A <- Omega / sum(X * Omega) # Normalize A to have <A,X>_F = 1, also see appendix B.2
  mse <- b^2 * max(svd(A)$d)^2 + sum(A^2)
  return(list(mse = mse, A = A))
}

