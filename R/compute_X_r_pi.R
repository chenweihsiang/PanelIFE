#' @title compute_X_r_pi
#'
#' @description compute_X_r_pi
#'
#' @param X_r_gamma X_r_gamma
#' @param mu mu
#' @param X X
#'
#' @return \code{X_r_pi}
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

