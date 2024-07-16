#' @title Extract R leading Principal Components
#'
#' @description This function extracts the "\code{R}" leading principal components out of the \eqn{NxT} matrix of residuals "\code{res}".
#'
#' @param res \eqn{N \times T} matrix of residuals
#' @param R A positive integer, indicates the number of interactive fixed effects in the estimation
#'
#' @return \eqn{N \times R} matrix \code{lambda} and \eqn{T \times R} matrix \code{f}, which are factor loadings and factors
#'
#' @noRd


principal_components <- function(res, R) {
  # Extract the "R" leading principal components out of the NxT matrix "res".
  # Output: NxR matrix lambda, TxR matrix f (=factor loadings and factors)

  # %COMMENT: within "LS_factor" it is guaranteed that T<=N, so below we
  # diagonalize a TxT matrix (not an NxN) matrix. When using this function
  # outside "LS_factor" one should check whether T<N or N>T and switch
  # dimensions accordingly, if neccessary.

  T <- ncol(res)
  eigen_decomp <- eigen(t(res) %*% res)
  Dsort <- sort(eigen_decomp$values)
  Ind <- order(eigen_decomp$values)
  f <- as.matrix(eigen_decomp$vectors[ , Ind[(T-R+1):T]])
  for (r in 1:R) {
    f[ , r] <- f[ , r] / norm(f[ , r], type = "2")
    if (mean(f[ , r]) < 0) {
      f[ , r] <- -f[ , r]
    }
  }
  lambda <- res %*% f
  return(list(lambda = lambda, f = f))
}

