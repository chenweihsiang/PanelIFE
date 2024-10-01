#' @title Truncate Symmetric Matrix
#'
#' @description This function truncates a symmetric matrix \code{A} to the \code{(M1+M2-1)} first diagonals
#'
#' @param A Symmetric matrix
#' @param M1 Integer
#' @param M2 Integer
#'
#' @return The \code{(M1+M2-1)} first diagonals of a symmetric matrix \code{A}
#'
#' @noRd


trunc <- function(A, M1, M2) {
  # Truncates a symmetric matrix A to the (M1+M2-1) first diagonals
  NN <- nrow(A)
  AT <- matrix(0, nrow = NN, ncol = NN)
  for(i in 1:NN) {
    for(j in max(i-M1+1, 1):min(i+M2-1, NN)) {
      if(max(i-M1+1, 1) <= min(i+M2-1, NN)) {
        AT[i, j] <- A[i, j]
      }
    }
  }
  return(AT)
}


