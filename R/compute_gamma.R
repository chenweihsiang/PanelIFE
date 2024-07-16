#' @title compute_gamma
#'
#' @description compute_gamma
#'
#' @param X_r_pi X_r_pi
#' @param ZZ ZZ
#' @param X X
#'
#' @return A list of \code{gamma} and \code{X_r_gamma}
#'
#' @noRd


compute_gamma <- function(X_r_pi, ZZ, X) {
  XX_r_pi <- as.vector(t(X_r_pi))
  gamma <- pracma::mldivide(ZZ, XX_r_pi)
  X_r_gamma <- X - t(array(ZZ %*% gamma, dim = c(dim(X)[2],dim(X)[1])))
  return(list(gamma = gamma, X_r_gamma = X_r_gamma))
}

