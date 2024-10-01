#' @title Extract R leading Principal Components
#'
#' @description This function extracts the "\code{R}" leading principal components out of the \eqn{N \times T} matrix of residuals \code{res}.
#'
#' @details Within \code{ls_factor} it is guaranteed that T <= N, so below we diagonalize a \eqn{T \times T} matrix (not an \eqn{N \times N}) matrix. When using this function outside \code{ls_factor}, one should check whether \eqn{T < N} or \eqn{N > T} and switch dimensions accordingly, if necessary.
#'
#' @param res \eqn{N \times T} matrix of residuals
#' @param R A positive integer, indicates the number of interactive fixed effects in the estimation
#'
#' @return A list of results, where
#' \itemize{
#'   \item \code{lambda} is the \eqn{N \times R} matrix of factor loadings
#'   \item \code{f} is the \eqn{T \times R} matrix of factors
#' }
#'
#' @noRd


principal_components <- function(res, R) {
  # Within "ls_factor" it is guaranteed that T <= N,
  #   so below we diagonalize a T*T matrix (not an N*N) matrix.
  # When using this function outside "ls_factor",
  #   one should check whether T < N or N > T
  #   and switch dimensions accordingly, if necessary.
  # ===
  T <- ncol(res)
  eigen_decomp <- eigen(t(res) %*% res)
  Dsort <- sort(eigen_decomp$values)
  Ind <- order(eigen_decomp$values)
  f <- as.matrix(eigen_decomp$vectors[ , Ind[(T-R+1):T]])
  for(r in 1:R) {
    f[ , r] <- f[ , r] / norm(f[ , r], type = "2")
    if(mean(f[ , r]) < 0) {
      f[ , r] <- -f[ , r]
    }
  }
  lambda <- res %*% f
  return(list(lambda = lambda, f = f))
}

