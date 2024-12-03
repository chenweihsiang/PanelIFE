#' @title Summarizing Least Squares Estimation of Linear Panel Data Models with Interactive Fixed Effects
#'
#' @description \code{summary} method for object of class \code{ls_factor} and returned object is of class {summary.ls_factor}
#'
#' @details
#' Summarizing the fitted \code{ls_factor} object and computing
#' \itemize{
#'   \item Estimates under different bias correction schemes.
#'   \item Variance-covariance matrix under different assumptions.
#'   \item Confidence intervals and their lengths under different settings.
#' }
#'
#' @param object The \code{ls_factor} fitted result of class \code{ls_factor}.
#' @param ... Confidence level \code{alpha} can be set for computing confidence intervals. Default is set to \eqn{0.05}.
#' @param x The summarized result of class \code{summary.ls_factor}.
#' @param digits Number of digits to use when printing. Default is set to \eqn{4}.
#' @param labels Option of whether or not to print labels in the summary. Default is set to \code{TRUE}. Only parameters estimates will be printed when it's set to \code{FALSE}.
#'
#' @return A list of results, where
#' \itemize{
#'   \item \code{beta} is the table of parameters estimates under different bias correction schemes.
#'   \item \code{CI} is the table of confidence intervals and their lengths under different settings.
#'   \item \code{var_cov} are the variance-covariance matrices under different assumptions.
#' }
#'
#' @method summary ls_factor
#' @export


summary.ls_factor <- function(object, ...) {
  if (!inherits(object, "ls_factor")) {
    stop("Invalid input")
  }
  if(!("alpha" %in% ls())) {
    # default confidence level
    alpha <- 0.05
  }
  # ===
  sum_res <- list(beta = cbind(stats::setNames(data.frame(t(cbind(object$beta,
                                                           object$beta - object$bcorr2 - object$bcorr3,
                                                           object$beta - object$bcorr1 - object$bcorr2 - object$bcorr3))),
                                        paste0("beta", c(1:dim(object$beta)[1]))),
                               data.frame(bias_correction = c("na", "static_biases", "all_biases"))),
                  var_cov = list(homoscedasticity = object$Vbeta1,
                                 heteroscedasticity = object$Vbeta2,
                                 serialcorrelation = object$Vbeta3),
                  CI = data.frame())
  sum_res$CI <- data.frame(est = rep(paste0("beta", c(1:dim(object$beta)[1])), each = 6),
                           alpha = alpha, CI_LB = NA, CI_UB = NA, CI_length = NA,
                           bias_correction = rep(c("na",
                                                   "static_biases",
                                                   "all_biases",
                                                   "na",
                                                   "static_biases",
                                                   "all_biases"),
                                                 dim(object$beta)[1]),
                           var_assumption = rep(c("homoscedasticity",
                                                  "homoscedasticity",
                                                  "homoscedasticity",
                                                  "heteroscedasticity",
                                                  "heteroscedasticity",
                                                  "heteroscedasticity"),
                                                dim(object$beta)[1]))
  for(i in c(1:dim(object$beta)[1])) {
    sum_res$CI$CI_LB[sum_res$CI$est == paste0("beta", i)] <-
      c(sum_res$beta[[paste0("beta", i)]] - sqrt(sum_res$var_cov$homoscedasticity[i, i]) * stats::qnorm(1 - alpha / 2),
        sum_res$beta[[paste0("beta", i)]] - sqrt(sum_res$var_cov$heteroscedasticity[i, i]) * stats::qnorm(1 - alpha / 2))
    sum_res$CI$CI_UB[sum_res$CI$est == paste0("beta", i)] <-
      c(sum_res$beta[[paste0("beta", i)]] + sqrt(sum_res$var_cov$homoscedasticity[i, i]) * stats::qnorm(1 - alpha / 2),
        sum_res$beta[[paste0("beta", i)]] + sqrt(sum_res$var_cov$heteroscedasticity[i, i]) * stats::qnorm(1 - alpha / 2))
    sum_res$CI$CI_length[sum_res$CI$est == paste0("beta", i)] <-
      sum_res$CI$CI_UB[sum_res$CI$est == paste0("beta", i)] - sum_res$CI$CI_LB[sum_res$CI$est == paste0("beta", i)]
  }
  # ===
  return(structure(sum_res, class = "summary.ls_factor"))
}


#' @rdname summary.ls_factor
#' @method print summary.ls_factor
#' @export


print.summary.ls_factor <- function(x, digits = 4, labels = TRUE, ...) {
  if (labels == TRUE) {
    est <- rbind(paste0("(beta ", c(1:(dim(x$beta)[2]-1)), ")"),
                 data.frame(matrix("NA", nrow = 3, ncol = ncol(x$beta)-1)))
    for(k in 1:(ncol(x$beta)-1)) {
      est[2:4, k] <- sprintf(paste0("%.", digits, "f"), x$beta[ , k])
    }
    est <- stats::setNames(cbind(data.frame(c("",
                                       " Without Correction:",
                                       " Static Biases Correction:",
                                       " All Biases Correction:")),
                          format(est, width = max(sapply(unlist(est), nchar)), justify = "centre")), NULL)
    # ===
    var1 <- rbind(c(paste0("(beta ", c(1:(dim(x$beta)[2]-1)), ")")),
                  data.frame(matrix("NA", nrow = nrow(x$var_cov$homoscedasticity), ncol = ncol(x$var_cov$homoscedasticity))))
    var2 <- data.frame(var1)
    for(i in 1:(nrow(var1)-1)) {
      for(j in 1:ncol(var1)) {
        var1[i+1, j] <- sprintf(paste0("%.", digits, "f"), x$var_cov$homoscedasticity[i, j])
        var2[i+1, j] <- sprintf(paste0("%.", digits, "f"), x$var_cov$heteroscedasticity[i, j])
      }
    }
    var1 <- stats::setNames(cbind(c(paste0(rep(" ", 26), collapse = ""), paste0("(beta ", c(1:(dim(x$beta)[2]-1)), ")")),
                           format(var1, width = max(sapply(unlist(var1), nchar)), justify = "centre")), NULL)
    var2 <- stats::setNames(cbind(c(paste0(rep(" ", 26), collapse = ""), paste0("(beta ", c(1:(dim(x$beta)[2]-1)), ")")),
                           format(var2, width = max(sapply(unlist(var2), nchar)), justify = "centre")), NULL)
    # ===
    cat("\n")
    cat("Parameters Estimates:")
    print(est, row.names = FALSE)
    cat("\nVariance-Covariance Matrix:\n")
    cat("  (1) Homoscedasticity:")
    print(var1, row.names = FALSE)
    cat("  (2) Heteroscedasticity:")
    print(var2, row.names = FALSE)
  } else {
    est <- rbind(paste0("(beta ", c(1:(dim(x$beta)[2]-1)), ")"),
                 data.frame(matrix("NA", nrow = 3, ncol = ncol(x$beta)-1)))
    for(k in 1:(ncol(x$beta)-1)) {
      est[2:4, k] <- sprintf(paste0("%.", digits, "f"), x$beta[ , k])
    }
    est <- stats::setNames(cbind(data.frame(c("",
                                       " Without Correction:",
                                       " Static Biases Correction:",
                                       " All Biases Correction:")),
                          format(est, width = max(sapply(unlist(est), nchar)), justify = "centre")), NULL)
    # ===
    cat("\n")
    cat("Parameters Estimates:")
    print(est, row.names = FALSE)
  }
}

