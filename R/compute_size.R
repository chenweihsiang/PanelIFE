#' @title Compute Size
#'
#' @description This function computes size with given parameter, estimates, standard errors, and confidence level
#'
#' @param beta Estimates
#' @param beta0 True parameter
#' @param se Standard errors
#' @param alpha It determines the confidence level \code{1 - alpha}
#'
#' @return A list of results, where
#' \itemize{
#'   \item \code{rp}
#'   \item \code{CI_length}
#'   \item \code{size}
#'   \item \code{cv}
#'   \item \code{CI_oracle_length}
#' }
#'
#' @noRd


compute_size <- function(beta, beta0, se, alpha) {
  n_e <- dim(beta)[2]
  if(n_e > dim(se)[2]) {
    se <- matrix(se, nrow = nrow(se), ncol = n_e)
  }
  # ===
  t_stat <- abs((beta - beta0) / se)
  # ===
  rp <- matrix(0, nrow = length(alpha), ncol = n_e)
  CI_length <- matrix(0, nrow = length(alpha), ncol = n_e)
  size <- matrix(0, nrow = length(alpha), ncol = n_e)
  CI_oracle_length <- matrix(0, nrow = length(alpha), ncol = n_e)
  cv <- matrix(0, nrow = length(alpha), ncol = n_e)
  # ===
  for(i_e in c(1:ncol(beta))) {
    rp[ , i_e] <- colMeans(sapply(alpha, function(x) {ifelse(matrix(t_stat[ , i_e], ncol = 1) >= stats::qnorm(1 - x / 2), 1, 0)}))
    CI_length[ , i_e] <- colMeans(2 * matrix(se[ , i_e], ncol = 1) %*% matrix(stats::qnorm(1 - alpha / 2), nrow = 1))
    size[ , i_e] <- colMeans(beta0 < matrix(beta[ , i_e], nrow = dim(beta)[1], ncol = length(alpha)) -
                               matrix(se[ , i_e], ncol = 1) %*% matrix(stats::qnorm(1 - alpha / 2), nrow = 1) |
                               beta0 > matrix(beta[ , i_e], nrow = dim(beta)[1], ncol = length(alpha)) +
                               matrix(se[ , i_e], ncol = 1) %*% matrix(stats::qnorm(1 - alpha / 2), nrow = 1))
    cv[ , i_e] <- as.numeric(stats::quantile(t_stat[ , i_e], 1 - alpha, type = 5))
    CI_oracle_length[ , i_e] <- colMeans(2 * matrix(se[ , i_e], ncol = 1) %*% t(cv[ , i_e]))
  }
  # ===
  return(list(rp = rp, CI_length = CI_length, size = size, cv = cv, CI_oracle_length = CI_oracle_length))
}

