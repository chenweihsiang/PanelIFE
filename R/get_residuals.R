#' @title Calculate Residuals
#'
#' @description This function calculates residuals \code{Y - beta * X}
#'
#' @param Y \eqn{N \times T} matrix
#' @param X \eqn{K \times N \times T} tensor
#' @param beta \eqn{K \times 1} vector
#'
#' @return \eqn{N \times T} matrix of residuals
#'
#' @noRd


get_residuals <- function(Y, X, beta) {
  res <- Y
  for (k in 1:dim(X)[1]) {
    res <- res - beta[k, 1] * X[k, , ]
  }
  return(res)
}


