#' @title Generate Statistics
#'
#' @description This function generates statistics with given parameter, estimates, standard errors, and confidence level
#'
#' @param beta \code{reps} times \code{n_est} (number of estimators) matrix of estimates
#' @param beta0 The true value of beta (scalar for now)
#' @param se \code{reps} times \code{n_est} (number of estimators) times \code{n_se} (number of std errs); could also be reps times 1 times n_se
#'
#' @return A list of results, where
#' \itemize{
#'   \item \code{rp}
#'   \item \code{length}
#'   \item \code{size}
#'   \item \code{oracle_length}
#'   \item \code{cv}
#' }
#'
#' @noRd


gen_stats <- function(beta, beta0, se) {
  stats <- list()
  stats$mean <- apply(as.matrix(beta), 2, mean)
  stats$med <- apply(as.matrix(beta), 2, stats::median)
  stats$std <- apply(as.matrix(beta), 2, function(x) {sqrt(sum((x - mean(x))^2) / length(x))})
  stats$iqr <- apply(as.matrix(beta), 2, function(x) {stats::IQR(x, type = 5)}) / 1.3490 # To be consistent with Matlab code, use type = 5
  # ===
  if(!missing(beta0)) {
    stats$bias <- stats$mean - beta0
    stats$med_bias <- stats$med - beta0
    stats$rmse <- sqrt(stats$bias^2 + stats$std^2)
  }
  # ===
  if(!missing(se)) {
    alpha <- c(0.10, 0.05 , 0.01, 0.01, 0.01)
    rp_arr <- array(0, dim = c(length(alpha), dim(beta)[2], dim(se)[3]))
    CI_length_arr <- array(0, dim = c(length(alpha), dim(beta)[2], dim(se)[3]))
    size_arr <- array(0, dim = c(length(alpha), dim(beta)[2], dim(se)[3]))
    CI_oracle_length_arr <- array(0, dim = c(length(alpha), dim(beta)[2], dim(se)[3]))
    cv_arr <- array(0, dim = c(length(alpha), dim(beta)[2], dim(se)[3]))
    # ===
    for(i_se in c(1:(dim(se)[3]))) {
      res <- compute_size(beta, beta0, matrix(se[ , , i_se], nrow = dim(se)[1], ncol = dim(se)[2]), alpha)
      rp_arr[ , , i_se] <- res$rp
      CI_length_arr[ , , i_se] <- res$CI_length
      size_arr[ , , i_se] <- res$size
      CI_oracle_length_arr[ , , i_se] <- res$CI_oracle_length
      cv_arr[ , , i_se] <- res$cv
    }
    # ===
    stats$rp <- rp_arr
    stats$length <- CI_length_arr
    stats$size <- size_arr
    stats$oracle_length <- CI_oracle_length_arr
    stats$cv <- cv_arr
  }
  # ===
  return(stats)
}

