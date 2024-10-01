#' @title Computing Gamma
#'
#' @description This function compute the \code{gamma}, which is the error component that can be correlated with \eqn{X} and \eqn{Z}. Also, it partials out the potential correlation of \eqn{X} with the estimation error \code{gamma}.
#'
#' @param X_r_pi \eqn{N \times T} matrix, \eqn{X} after partialling out the potential correlation with the estimation error and perform the nuclear norm regularization
#' @param ZZ \eqn{(NT) \times (K-1)} matrix
#' @param X \eqn{N \times T} matrix
#'
#' @return A list of results, where
#' \itemize{
#'   \item \code{gamma} is the error component that can be correlated with \eqn{X} and \eqn{Z}
#'   \item \code{X_r_gamma} is \eqn{X} after partialling out the potential correlation with the estimation error \code{Gamma}
#' }
#'
#' @noRd


compute_gamma <- function(X_r_pi, ZZ, X) {
  XX_r_pi <- as.vector(t(X_r_pi))
  gamma <- pracma::mldivide(ZZ, XX_r_pi)
  X_r_gamma <- X - t(array(ZZ %*% gamma, dim = c(dim(X)[2],dim(X)[1])))
  return(list(gamma = gamma, X_r_gamma = X_r_gamma))
}

