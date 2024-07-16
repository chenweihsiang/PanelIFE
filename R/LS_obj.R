#' @title Calculate LS Objective
#'
#' @description This function calculates LS objective, i.e. SSR/N/T of \code{Y - beta * X} after subtracting \code{R} largest principal components. Within "\code{LS_factor}" it is guaranteed that \eqn{T \leq N}, so below we diagonalize a \eqn{TxT} matrix (not an \eqn{NxN}) matrix. When using this function outside "\code{LS_factor}" one should check whether \eqn{T<N} or \eqn{N>T} and switch dimensions accordingly, if necessary.
#'
#' @param beta \eqn{K \times 1} vector
#' @param Y \eqn{N \times T} matrix
#' @param X \eqn{K \times N \times T} tensor
#' @param R A positive integer, indicates the number of interactive fixed effects in the estimation
#'
#' @return Scalar objective function
#'
#' @noRd


LS_obj <- function(beta, Y, X, R) {
  # Calculate LS objective, i.e. SSR/N/T of Y-beta*X after subtracting R
  # largest principal components.
  # INPUT: Kx1 beta, NxT Y, KxNxT X, integer R>=0
  # OUTPUT: scalar objective function

  # COMMENT: within "LS_factor" it is guaranteed that T<=N, so below we
  # diagonalize a TxT matrix (not an NxN) matrix. When using this function
  # outside "LS_factor" one should check whether T<N or N>T and switch
  # dimensions accordingly, if neccessary.

  res <- get_residuals(Y, X, beta)
  ev <- sort(eigen(t(res) %*% res)$values)
  obj <- sum(ev[1:(dim(Y)[2] - R)]) / dim(Y)[1] / dim(Y)[2]
  return(obj)
}

