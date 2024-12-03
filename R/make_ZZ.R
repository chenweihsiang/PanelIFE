#' @title Stacking Tensor to Matrix
#'
#' @description This function converts a \eqn{K \times N \times T} tensor to \eqn{NT \times K} matrix by stacking them.
#'
#' @param Z \eqn{K \times N \times T} tensor
#'
#' @return \eqn{NT \times K} matrix
#'
#' @examples
#' # For example, a 2*3*4 tensor Z[1, , ] = [ 1  7 13 19
#' #                                          3  9 15 21
#' #                                          5 11 17 23]
#' #                             Z[2, , ] = [ 2  8 14 20
#' #                                          4 10 16 22
#' #                                          6 12 18 24]
#' # will then become a 12*2 matrix ZZ = [ 1  7 13 19  3  9 15 21  5 11 17 23
#' #                                       2  8 14 20  4 10 16 22  6 12 18 24]'
#' # Z <- array(c(1:(2*3*4)), dim = c(2, 3, 4))
#' # ZZ <- make_ZZ(Z)
#'
#' @noRd


make_ZZ <- function(Z) {
  N <- dim(Z)[2]
  T <- dim(Z)[3]
  ZZ <- matrix(0, N*T, dim(Z)[1])
  for (k in 1:dim(Z)[1]) {
    ZZ[ , k] <- as.vector(t(Z[k, , ]))
  }
  return(ZZ)
}

