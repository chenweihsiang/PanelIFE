#' @title Computing X after Partialling Out and with the Nuclear Norm Regularization
#'
#' @description This function compute the \eqn{X} after partialling out the potential correlation with the estimation error and perform the nuclear norm regularization.
#'
#' @param X_r_gamma \eqn{N \times T} matrix, \eqn{X} after partialling out the potential correlation with the estimation error \code{gamma}
#' @param mu Scalar, the penalty on the nuclear norm
#' @param X \eqn{N \times T} matrix
#'
#' @return \eqn{N \times T} matrix, \eqn{X} after partialling out the potential correlation with the estimation error and perform the nuclear norm regularization
#'
#' @noRd


compute_X_r_pi <- function(X_r_gamma, mu, X) {
  svd_res <- svd(X_r_gamma) # SVD: X_r_gamma = V S W'
  V <- svd_res$u
  S <- diag(svd_res$d)
  W <- svd_res$v
  X_r_pi <- X - V %*% diag(pmax(diag(S) - mu, 0)) %*% t(W)
  return(X_r_pi)
}

