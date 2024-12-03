#' @title Calculate Least Squares Objective
#'
#' @description This function calculates least squares objective, i.e. SSR/N/T of \code{Y - beta * X} after subtracting \code{R} largest principal components.
#'
#' @details Within \code{ls_factor} it is guaranteed that T <= N, so below we diagonalize a \eqn{T \times T} matrix (not an \eqn{N \times N}) matrix. When using this function outside \code{ls_factor}, one should check whether \eqn{T < N} or \eqn{N > T} and switch dimensions accordingly, if necessary.
#'
#' @param beta \eqn{K \times 1} vector
#' @param Y \eqn{N \times T} matrix
#' @param X \eqn{K \times N \times T} tensor
#' @param R A positive integer, indicates the number of interactive fixed effects in the estimation
#'
#' @return Scalar objective function
#'
#' @noRd


ls_obj <- function(beta, Y, X, R) {
  res <- get_residuals(Y, X, beta)
  ev <- sort(eigen(t(res) %*% res)$values)
  obj <- sum(ev[1:(dim(Y)[2] - R)]) / dim(Y)[1] / dim(Y)[2]
  return(obj)
}

