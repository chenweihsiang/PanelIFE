#' @title Summarizing Robust Estimation and Inference in Panels with Interactive Fixed Effects
#'
#' @description \code{summary} method for object of class \code{honest_weak_factors} and returned object is of class {summary.honest_weak_factors}
#'
#' @details
#' Summarizing the fitted \code{honest_weak_factors} object and computing
#' \itemize{
#'   \item Debiased parameter estimate.
#'   \item Standard error.
#'   \item Worst-case bias of the estimator.
#'   \item With and without bias-aware confidence intervals.
#' }
#'
#' @param object The \code{honest_weak_factors} fitted result of class \code{honest_weak_factors}.
#' @param ... Confidence level \code{alpha} can be set for computing confidence intervals. Default is set to \eqn{0.05}.
#' @param x The summarized result of class \code{summary.ls_factor}.
#' @param digits Number of digits to use when printing. Default is set to \eqn{4}.
#' @param labels Option of whether or not to print labels in the summary. Default is set to \code{TRUE}. Only parameters estimates will be printed when it's set to \code{FALSE}.
#'
#' @return A list of results, where
#' \itemize{
#'   \item \code{beta} is the debiased parameter estimate.
#'   \item \code{se} is the standard errors of the parameters estimator.
#'   \item \code{bias} is the worst-case bias of the estimator.
#'   \item \code{CI} is the bias-aware confidence interval.
#'   \item \code{CI_unadj} is the confidence interval without bias adjustment.
#'   \item \code{A} is the choice of weight matrix.
#' }
#'
#' @method summary honest_weak_factors
#' @export


summary.honest_weak_factors <- function(object, ...) {
  if (!inherits(object, "honest_weak_factors")) {
    stop("Invalid input")
  }
  if(!("alpha" %in% ls())) {
    # default confidence level
    alpha <- 0.05
  }
  # ===
  sum_res <- list(beta = object$beta,
                  se = object$se,
                  bias = object$bias,
                  CI = c(object$LB, object$UB),
                  CI_unadj = c(object$LB_unadj, object$UB_unadj),
                  A = object$A)
  # ===
  class(sum_res) <- "summary.honest_weak_factors"
  return(sum_res)
}


#' @rdname summary.honest_weak_factors
#' @method print summary.honest_weak_factors
#' @export


print.summary.honest_weak_factors <- function(x, digits = 4, labels = TRUE, ...) {
  if (labels == TRUE) {
    cat("\n")
    cat("Parameter Estimate:  ", round(x$beta, digits), "\n", sep = "")
    cat("Standard Error:      ", round(x$se, digits), "\n", sep = "")
    cat("Worst-Case Bias:     ", round(x$bias, digits), "\n", sep = "")
    cat("Confidence Interval: ", "\n", sep = "")
    cat("  (1) Bias-Aware:         (", round(x$CI[1], digits), " , ", round(x$CI[2], digits), ")", "\n", sep = "")
    cat("  (2) Without Adjustment: (", round(x$CI_unadj[1], digits), " , ", round(x$CI_unadj[2], digits), ")", "\n", sep = "")
  } else {
    cat("\n")
    cat("Parameter Estimate:  ", round(x$beta, digits), "\n", sep = "")
  }
}

