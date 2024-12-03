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
                  CI_LB = object$LB,
                  CI_UB = object$UB,
                  A = object$A,
                  parameter = object$parameter)
  # ===
  return(structure(sum_res, class = "summary.honest_weak_factors"))
}


#' @rdname summary.honest_weak_factors
#' @method print summary.honest_weak_factors
#' @export


print.summary.honest_weak_factors <- function(x, digits = 4, labels = TRUE, ...) {
  if (labels == TRUE) {
    cat("\n")
    cat("Parameter Estimate:  ", round(x$beta, digits), "\n", sep = "")
    cat("Standard Error:      ", round(x$se, digits), "\n", sep = "")
    # ===
    df_bias <- stats::setNames(rbind(c(" # Weak Factors", paste0("Bias (", c(1:length(x$beta)), ")")),
                              data.frame(matrix("NA", nrow = nrow(x$bias), ncol = ncol(x$bias)))),
                        NULL)
    df_bias[2:nrow(df_bias), 1] <- paste0(x$bias[, 1])
    for(k in 1:length(x$beta)) {
      df_bias[2:nrow(df_bias), k+1] <- sprintf(paste0("%.", digits, "f"), x$bias[ , k+1])
    }
    # ===
    df_CI <- stats::setNames(rbind(c("", c(rbind(rep(paste0((1-x$parameter$alpha)*100, "% ", "CI"), length(x$beta)),
                                          paste0("(", 1:length(x$beta), ")")))),
                            c(" # Weak Factors", rbind(rep("Lower Bound", length(x$beta)), rep("Upper Bound", length(x$beta)))),
                            data.frame(matrix("NA", nrow = nrow(x$bias), ncol = ncol(x$bias)*2-1))),
                      NULL)
    df_CI[3:nrow(df_CI), 1] <- paste0(x$bias[, 1])
    for(k in 1:length(x$beta)) {
      df_CI[3:nrow(df_CI), k*2] <- sprintf(paste0("%.", digits, "f"), x$CI_LB[ , k+1])
      df_CI[3:nrow(df_CI), k*2+1] <- sprintf(paste0("%.", digits, "f"), x$CI_UB[ , k+1])
    }
    CI_text <- paste0((1-x$parameter$alpha)*100, "% CI")
    max_width <- max(sapply(unlist(df_CI), nchar))
    df_CI[1, ] <- c("",
                    rbind(paste0(paste0(rep(" ", max_width - nchar(CI_text)), collapse = ""), rep(CI_text, length(x$beta))),
                          paste0("(", 1:length(x$beta), ")", paste0(rep(" ", max_width - nchar(length(x$beta)) - 2), collapse = ""))))
    # ===
    cat("\n")
    cat("Worst-Case Bias:     ")
    print(format(df_bias, width = max(sapply(unlist(df_bias), nchar)), justify = "centre"), row.names = FALSE)
    cat("\n")
    cat("Confidence Interval: ")
    print(format(df_CI, width = max(sapply(unlist(df_CI), nchar)), justify = "centre"), row.names = FALSE)
    cat("\n")
    cat("(# weak factors = 0 indicates CI without adjustment; > 0 indicates (robust) bias-aware CI)\n")
  } else {
    cat("\n")
    cat("Parameter Estimate:  ", round(x$beta, digits), "\n", sep = "")
  }
}

