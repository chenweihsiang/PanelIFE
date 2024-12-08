#' @title Empirical Dataset (Divorce Rate in the USA from 1956 to 1988)
#'
#' @description Data of the divorce rate in the USA from 1956 to 1988, originally from Friedberg (1998, "Did Unilateral Divorce Raise Divorce Rates? Evidence from Panel Data") and Wolfers (2006, "Did Unilateral Divorce Laws Raise Divorce Rates? A Reconciliation and New Results").
#' Here, this empirical dataset follows Kim and Oka (2014, "Divorce Law Reforms and Divorce Rates in the USA: An Interactive Fixed‐Effects Approach") and use their data to construct a balanced panel with \code{N = 48} states and \code{T = 33} years.
#'
#' @format A \code{data.frame} with 1,650 rows and 88 variables:
#' \describe{
#'   \item{st}{Two-letter state code}
#'   \item{id_st}{Numeric ID for the states}
#'   \item{year}{Year of observation}
#'   \item{div_rate_rev01}{Divorces per 1000 people, 1956-1998}
#'   \item{div_rate_rev02}{Divorces per 1000 people, 1956-1998}
#'   \item{stpop}{State population}
#'   \item{unilateral}{Dummy for unilateral law in Friedberg (1998)}
#'   \item{dyn_uni2}{Dynamic treatment effects (years 1-2): dummy for year of unilateral law in Friedberg (1998)}
#'   \item{dyn_uni3}{Dynamic treatment effects (years 3-4): dummy for year of unilateral law in Friedberg (1998)}
#'   \item{dyn_uni4}{Dynamic treatment effects (years 5-6): dummy for year of unilateral law in Friedberg (1998)}
#'   \item{dyn_uni5}{Dynamic treatment effects (years 7-8): dummy for year of unilateral law in Friedberg (1998)}
#'   \item{dyn_uni6}{Dynamic treatment effects (years 9-10): dummy for year of unilateral law in Friedberg (1998)}
#'   \item{dyn_uni7}{Dynamic treatment effects (years 11-12): dummy for year of unilateral law in Friedberg (1998)}
#'   \item{dyn_uni8}{Dynamic treatment effects (years 13-14): dummy for year of unilateral law in Friedberg (1998)}
#'   \item{dyn_uni9}{Dynamic treatment effects (years 15+): dummy for year of unilateral law in Friedberg (1998)}
#'   \item{dyn_gruber_2}{Dynamic treatment effects (years 1-2): dummy for year of divorce law reform according to Gruber (2004)}
#'   \item{dyn_gruber_3}{Dynamic treatment effects (years 3-4): dummy for year of divorce law reform according to Gruber (2004)}
#'   \item{dyn_gruber_4}{Dynamic treatment effects (years 5-6): dummy for year of divorce law reform according to Gruber (2004)}
#'   \item{dyn_gruber_5}{Dynamic treatment effects (years 7-8): dummy for year of divorce law reform according to Gruber (2004)}
#'   \item{dyn_gruber_6}{Dynamic treatment effects (years 9-10): dummy for year of divorce law reform according to Gruber (2004)}
#'   \item{dyn_gruber_7}{Dynamic treatment effects (years 11-12): dummy for year of divorce law reform according to Gruber (2004)}
#'   \item{dyn_gruber_8}{Dynamic treatment effects (years 13-14): dummy for year of divorce law reform according to Gruber (2004)}
#'   \item{dyn_gruber_9}{Dynamic treatment effects (years 15+): dummy for year of divorce law reform according to Gruber (2004)}
#'   \item{dyn_johnson_2}{Dynamic treatment effects (years 1-2): dummy for year of divorce law reform according to Johnson and Mazingo (2000)}
#'   \item{dyn_johnson_3}{Dynamic treatment effects (years 3-4): dummy for year of divorce law reform according to Johnson and Mazingo (2000)}
#'   \item{dyn_johnson_4}{Dynamic treatment effects (years 5-6): dummy for year of divorce law reform according to Johnson and Mazingo (2000)}
#'   \item{dyn_johnson_5}{Dynamic treatment effects (years 7-8): dummy for year of divorce law reform according to Johnson and Mazingo (2000)}
#'   \item{dyn_johnson_6}{Dynamic treatment effects (years 9-10): dummy for year of divorce law reform according to Johnson and Mazingo (2000)}
#'   \item{dyn_johnson_7}{Dynamic treatment effects (years 11-12): dummy for year of divorce law reform according to Johnson and Mazingo (2000)}
#'   \item{dyn_johnson_8}{Dynamic treatment effects (years 13-14): dummy for year of divorce law reform according to Johnson and Mazingo (2000)}
#'   \item{dyn_johnson_9}{Dynamic treatment effects (years 15+): dummy for year of divorce law reform according to Johnson and Mazingo (2000)}
#'   \item{dyn_mechoulan_2}{Dynamic treatment effects (years 1-2): dummy for year of divorce law reform according to Mechoulan (2001)}
#'   \item{dyn_mechoulan_3}{Dynamic treatment effects (years 3-4): dummy for year of divorce law reform according to Mechoulan (2001)}
#'   \item{dyn_mechoulan_4}{Dynamic treatment effects (years 5-6): dummy for year of divorce law reform according to Mechoulan (2001)}
#'   \item{dyn_mechoulan_5}{Dynamic treatment effects (years 7-8): dummy for year of divorce law reform according to Mechoulan (2001)}
#'   \item{dyn_mechoulan_6}{Dynamic treatment effects (years 9-10): dummy for year of divorce law reform according to Mechoulan (2001)}
#'   \item{dyn_mechoulan_7}{Dynamic treatment effects (years 11-12): dummy for year of divorce law reform according to Mechoulan (2001)}
#'   \item{dyn_mechoulan_8}{Dynamic treatment effects (years 13-14): dummy for year of divorce law reform according to Mechoulan (2001)}
#'   \item{dyn_mechoulan_9}{Dynamic treatment effects (years 15+): dummy for year of divorce law reform according to Mechoulan (2001)}
#'   \item{dyn_ellmanlohr1_2}{Dynamic treatment effects (years 1-2): dummy for year of divorce law reform according to Ellman and Lohr (1998) (definition a)}
#'   \item{dyn_ellmanlohr1_3}{Dynamic treatment effects (years 3-4): dummy for year of divorce law reform according to Ellman and Lohr (1998) (definition a)}
#'   \item{dyn_ellmanlohr1_4}{Dynamic treatment effects (years 5-6): dummy for year of divorce law reform according to Ellman and Lohr (1998) (definition a)}
#'   \item{dyn_ellmanlohr1_5}{Dynamic treatment effects (years 7-8): dummy for year of divorce law reform according to Ellman and Lohr (1998) (definition a)}
#'   \item{dyn_ellmanlohr1_6}{Dynamic treatment effects (years 9-10): dummy for year of divorce law reform according to Ellman and Lohr (1998) (definition a)}
#'   \item{dyn_ellmanlohr1_7}{Dynamic treatment effects (years 11-12): dummy for year of divorce law reform according to Ellman and Lohr (1998) (definition a)}
#'   \item{dyn_ellmanlohr1_8}{Dynamic treatment effects (years 13-14): dummy for year of divorce law reform according to Ellman and Lohr (1998) (definition a)}
#'   \item{dyn_ellmanlohr1_9}{Dynamic treatment effects (years 15+): dummy for year of divorce law reform according to Ellman and Lohr (1998) (definition a)}
#'   \item{dyn_ellmanlohr2_2}{Dynamic treatment effects (years 1-2): dummy for year of divorce law reform according to Ellman and Lohr (1998) (definition b)}
#'   \item{dyn_ellmanlohr2_3}{Dynamic treatment effects (years 3-4): dummy for year of divorce law reform according to Ellman and Lohr (1998) (definition b)}
#'   \item{dyn_ellmanlohr2_4}{Dynamic treatment effects (years 5-6): dummy for year of divorce law reform according to Ellman and Lohr (1998) (definition b)}
#'   \item{dyn_ellmanlohr2_5}{Dynamic treatment effects (years 7-8): dummy for year of divorce law reform according to Ellman and Lohr (1998) (definition b)}
#'   \item{dyn_ellmanlohr2_6}{Dynamic treatment effects (years 9-10): dummy for year of divorce law reform according to Ellman and Lohr (1998) (definition b)}
#'   \item{dyn_ellmanlohr2_7}{Dynamic treatment effects (years 11-12): dummy for year of divorce law reform according to Ellman and Lohr (1998) (definition b)}
#'   \item{dyn_ellmanlohr2_8}{Dynamic treatment effects (years 13-14): dummy for year of divorce law reform according to Ellman and Lohr (1998) (definition b)}
#'   \item{dyn_ellmanlohr2_9}{Dynamic treatment effects (years 15+): dummy for year of divorce law reform according to Ellman and Lohr (1998) (definition b)}
#'   \item{dyn_brinigbuckley_2}{Dynamic treatment effects (years 1-2): dummy for year of divorce law reform according to Brinig and Buckley (1998)}
#'   \item{dyn_brinigbuckley_3}{Dynamic treatment effects (years 3-4): dummy for year of divorce law reform according to Brinig and Buckley (1998)}
#'   \item{dyn_brinigbuckley_4}{Dynamic treatment effects (years 5-6): dummy for year of divorce law reform according to Brinig and Buckley (1998)}
#'   \item{dyn_brinigbuckley_5}{Dynamic treatment effects (years 7-8): dummy for year of divorce law reform according to Brinig and Buckley (1998)}
#'   \item{dyn_brinigbuckley_6}{Dynamic treatment effects (years 9-10): dummy for year of divorce law reform according to Brinig and Buckley (1998)}
#'   \item{dyn_brinigbuckley_7}{Dynamic treatment effects (years 11-12): dummy for year of divorce law reform according to Brinig and Buckley (1998)}
#'   \item{dyn_brinigbuckley_8}{Dynamic treatment effects (years 13-14): dummy for year of divorce law reform according to Brinig and Buckley (1998)}
#'   \item{dyn_brinigbuckley_9}{Dynamic treatment effects (years 15+): dummy for year of divorce law reform according to Brinig and Buckley (1998)}
#'   \item{dyn_nakonezny_2}{Dynamic treatment effects (years 1-2): dummy for year of divorce law reform according to Nakonezny, Shull and Rodgers (1995)}
#'   \item{dyn_nakonezny_3}{Dynamic treatment effects (years 3-4): dummy for year of divorce law reform according to Nakonezny, Shull and Rodgers (1995)}
#'   \item{dyn_nakonezny_4}{Dynamic treatment effects (years 5-6): dummy for year of divorce law reform according to Nakonezny, Shull and Rodgers (1995)}
#'   \item{dyn_nakonezny_5}{Dynamic treatment effects (years 7-8): dummy for year of divorce law reform according to Nakonezny, Shull and Rodgers (1995)}
#'   \item{dyn_nakonezny_6}{Dynamic treatment effects (years 9-10): dummy for year of divorce law reform according to Nakonezny, Shull and Rodgers (1995)}
#'   \item{dyn_nakonezny_7}{Dynamic treatment effects (years 11-12): dummy for year of divorce law reform according to Nakonezny, Shull and Rodgers (1995)}
#'   \item{dyn_nakonezny_8}{Dynamic treatment effects (years 13-14): dummy for year of divorce law reform according to Nakonezny, Shull and Rodgers (1995)}
#'   \item{dyn_nakonezny_9}{Dynamic treatment effects (years 15+): dummy for year of divorce law reform according to Nakonezny, Shull and Rodgers (1995)}
#'   \item{divx1}{Friedberg (1998)'s dummies for coding breaks; see appendix of the paper}
#'   \item{divx2}{Friedberg (1998)'s dummies for coding breaks; see appendix of the paper}
#'   \item{divx3}{Friedberg (1998)'s dummies for coding breaks; see appendix of the paper}
#'   \item{divx4}{Friedberg (1998)'s dummies for coding breaks; see appendix of the paper}
#'   \item{divx5}{Friedberg (1998)'s dummies for coding breaks; see appendix of the paper}
#'   \item{divx6}{Friedberg (1998)'s dummies for coding breaks; see appendix of the paper}
#'   \item{divx7}{Friedberg (1998)'s dummies for coding breaks; see appendix of the paper}
#'   \item{divx8}{Friedberg (1998)'s dummies for coding breaks; see appendix of the paper}
#'   \item{divx9}{Friedberg (1998)'s dummies for coding breaks; see appendix of the paper}
#'   \item{divx10}{Friedberg (1998)'s dummies for coding breaks; see appendix of the paper}
#'   \item{divx11}{Friedberg (1998)'s dummies for coding breaks; see appendix of the paper}
#'   \item{divx12}{Friedberg (1998)'s dummies for coding breaks; see appendix of the paper}
#'   \item{divx13}{Friedberg (1998)'s dummies for coding breaks; see appendix of the paper}
#'   \item{divx14}{Friedberg (1998)'s dummies for coding breaks; see appendix of the paper}
#'   \item{divx15}{Friedberg (1998)'s dummies for coding breaks; see appendix of the paper}
#'   \item{divx16}{Friedberg (1998)'s dummies for coding breaks; see appendix of the paper}
#'   \item{divx17}{Friedberg (1998)'s dummies for coding breaks; see appendix of the paper}
#' }
#'
#' @references For the detail of the data and the contruction of the balanced panel, see Friedberg (1998, "Did Unilateral Divorce Raise Divorce Rates? Evidence from Panel Data"), Wolfers (2006, "Did Unilateral Divorce Laws Raise Divorce Rates? A Reconciliation and New Results"), and Kim and Oka (2014, "Divorce Law Reforms and Divorce Rates in the USA: An Interactive Fixed‐Effects Approach").
#'
#' @note
#' \strong{Examples} section provides replication code for Table 5 and 9 in Armstrong, Weidner, and Zeleneev (2024, "Robust Estimation and Inference in Panels with Interactive Fixed Effects").
#'
#' @examples
#' # Following replication requires some time for computation
#'
#' # # Replication of Table 5 in Armstrong, Weidner, and Zeleneev (2024)
#' # ## Data cleaning
#' # df <- PanelIFE::empirical_data[ , -1]       # Load data without two-letter state code
#' # df <- df[df$id_st != 16 & df$id_st != 33, ] # Drop IN (id == 16) and NM (id == 33)
#' # XX <- df[ , grep("dyn_uni", colnames(df))]
#' # # ===
#' # T <- max(df$year) - min(df$year) + 1
#' # N <- nrow(df) / T
#' # K <- ncol(XX)
#' # # ===
#' # Y <- matrix(df$div_rate_rev02, nrow = N, ncol = T, byrow = TRUE)
#' # X <- array(0, dim = c(K, N, T))
#' # for(k in 1:K) {
#' #   X[k, , ] <- matrix(XX[ , k], nrow = N, ncol = T, byrow = TRUE)
#' # }
#' # # ===
#' # D <- array(0, dim = c(1, N, T))             # Using one regressor only
#' # D[1, , ] <- apply(X, c(2, 3), sum)
#' # X <- D
#' # K <- 1
#' # ## Setting estimation parameters
#' # set.seed(1)
#' # beta0 <- rep(0, K)
#' # alpha <- 0.05                               # 1 - confidence level
#' # Rmax <- 6                                   # Maximum number of the factors
#' # lambda_known <- matrix(rep(1, N))           # Standard time effects
#' # f_known <- cbind(rep(1, T), (1:T), (1:T)^2) # Individual effects + time trend + time trend^2
#' # ## Least squares estimation of linear panel data model with interactive fixed effects (LS factors)
#' # LS_summary <- setNames(data.frame(matrix(NA, nrow = Rmax, ncol = 4)),
#' #                        c("R", "beta", "CI_LB", "CI_UB"))
#' # for(R in 1:Rmax) {
#' #   res <- ls_factor(Y = Y, X = X, R = R,
#' #                    lambda_known = lambda_known, f_known = f_known,
#' #                    report = "silent", precision_beta = 10^(-8), method = "m1",
#' #                    start = beta0, repMIN = 30, repMAX = 300, M1 = 2, M2 = 2)
#' #   LS_summary$R[R] <- R
#' #   LS_summary$beta[R] <- res$beta - res$bcorr2 - res$bcorr3
#' #   LS_summary$CI_LB[R] <- LS_summary$beta[R] - sqrt(diag(res$Vbeta2)) * qnorm(1 - alpha/2)
#' #   LS_summary$CI_UB[R] <- LS_summary$beta[R] + sqrt(diag(res$Vbeta2)) * qnorm(1 - alpha/2)
#' # }
#' # ## Robust estimation and inference in panels with interactive fixed effects (honest weak factors)
#' # Debiased_summary <- setNames(data.frame(matrix(NA, nrow = Rmax * 3, ncol = 5)),
#' #                              c("R", "Rw", "beta", "CI_LB", "CI_UB"))
#' # for(R in 1:Rmax) {
#' # res <- honest_weak_factors(Y = Y, X = X, R = R,
#' #                            Gamma_LS = NULL, alpha = 0.05, clustered_se = FALSE,
#' #                            lambda_known = lambda_known, f_known = f_known,
#' #                            itermax = 75, reltol = 10^(-6))
#' #   Debiased_summary$R[(3*R-2):(3*R)] <- R
#' #   Debiased_summary$Rw[(3*R-2):(3*R)] <- c(0, 1, R)
#' #   Debiased_summary$beta[(3*R-2):(3*R)] <- res$beta
#' #   Debiased_summary$CI_LB[(3*R-2):(3*R)] <- res$LB[c(1:2, R+1), 2]
#' #   Debiased_summary$CI_UB[(3*R-2):(3*R)] <- res$UB[c(1:2, R+1), 2]
#' # }
#' # ## Print results
#' # print(round(LS_summary, 3), row.names = FALSE)
#' # print(round(Debiased_summary, 3), row.names = FALSE)
#'
#' # # Replication of Table 10 in Armstrong, Weidner, and Zeleneev (2024)
#' # ## Data cleaning
#' # df <- PanelIFE::empirical_data[ , -1]       # Load data without two-letter state code
#' # df <- df[df$id_st != 16 & df$id_st != 33, ] # Drop IN (id == 16) and NM (id == 33)
#' # XX <- df[ , grep("dyn_uni", colnames(df))]
#' # # ===
#' # T <- max(df$year) - min(df$year) + 1
#' # N <- nrow(df) / T
#' # K <- ncol(XX)
#' # # ===
#' # Y <- matrix(df$div_rate_rev02, nrow = N, ncol = T, byrow = TRUE)
#' # X <- array(0, dim = c(K, N, T))
#' # for(k in 1:K) {
#' #   X[k, , ] <- matrix(XX[ , k], nrow = N, ncol = T, byrow = TRUE)
#' # }
#' # # ===
#' # D <- array(0, dim = c(4, N, T))             # Using four regressors
#' # D[1, , ] <- apply(X[1:2, , ], c(2, 3), sum)
#' # D[2, , ] <- apply(X[3:4, , ], c(2, 3), sum)
#' # D[3, , ] <- apply(X[5:6, , ], c(2, 3), sum)
#' # D[4, , ] <- apply(X[7:8, , ], c(2, 3), sum)
#' # X <- D
#' # K <- 4
#' # ## Setting estimation parameters
#' # set.seed(1)
#' # beta0 <- rep(0, K)
#' # alpha <- 0.05                               # 1 - confidence level
#' # Rmax <- 6                                   # Maximum number of the factors
#' # lambda_known <- matrix(rep(1, N))           # Standard time effects
#' # f_known <- cbind(rep(1, T), (1:T), (1:T)^2) # Individual effects + time trend + time trend^2
#' # ## Least squares estimation of linear panel data model with interactive fixed effects (LS factors)
#' # LS_summary <- setNames(data.frame(matrix(NA, nrow = Rmax * K, ncol = 5)),
#' #                        c("R", "k", "beta", "CI_LB", "CI_UB"))
#' # for(R in 1:Rmax) {
#' #   res <- ls_factor(Y = Y, X = X, R = R,
#' #                    lambda_known = lambda_known, f_known = f_known,
#' #                    report = "silent", precision_beta = 10^(-8), method = "m1",
#' #                    start = beta0, repMIN = 30, repMAX = 300, M1 = 2, M2 = 2)
#' #   LS_summary$k[(K*R-(K-1)):(K*R)] <- 1:K
#' #   LS_summary$R[(K*R-(K-1)):(K*R)] <- R
#' #   LS_summary$beta[(K*R-(K-1)):(K*R)] <- as.numeric(res$beta - res$bcorr2 - res$bcorr3)
#' #   LS_summary$CI_LB[(K*R-(K-1)):(K*R)] <-
#' #     LS_summary$beta[(K*R-(K-1)):(K*R)] - sqrt(diag(res$Vbeta2)) * qnorm(1 - alpha/2)
#' #   LS_summary$CI_UB[(K*R-(K-1)):(K*R)] <-
#' #     LS_summary$beta[(K*R-(K-1)):(K*R)] + sqrt(diag(res$Vbeta2)) * qnorm(1 - alpha/2)
#' # }
#' # ## Robust estimation and inference in panels with interactive fixed effects (honest weak factors)
#' # Debiased_summary <- setNames(data.frame(matrix(NA, nrow = Rmax * K * 3, ncol = 6)),
#' #                              c("R", "k", "Rw", "beta", "CI_LB", "CI_UB"))
#' # for(R in 1:Rmax) {
#' #   res <- honest_weak_factors(Y = Y, X = X, R = R,
#' #                              Gamma_LS = NULL, alpha = 0.05, clustered_se = FALSE,
#' #                              lambda_known = lambda_known, f_known = f_known,
#' #                              itermax = 50, reltol = 10^(-4))
#' #   Debiased_summary$k[(3*K*R-(3*K-1)):(3*K*R)] <- rep(1:K, each = 3)
#' #   Debiased_summary$R[(3*K*R-(3*K-1)):(3*K*R)] <- R
#' #   Debiased_summary$Rw[(3*K*R-(3*K-1)):(3*K*R)] <- rep(c(0, 1, R), K)
#' #   Debiased_summary$beta[(3*K*R-(3*K-1)):(3*K*R)] <- rep(res$beta, each = 3)
#' #   Debiased_summary$CI_LB[(3*K*R-(3*K-1)):(3*K*R)] <- c(as.matrix(res$LB[c(1:2, R+1), -1]))
#' #   Debiased_summary$CI_UB[(3*K*R-(3*K-1)):(3*K*R)] <- c(as.matrix(res$UB[c(1:2, R+1), -1]))
#' # }
#' # ## Print results
#' # print(round(LS_summary, 3), row.names = FALSE)
#' # print(round(Debiased_summary, 3), row.names = FALSE)
#'


"empirical_data"

