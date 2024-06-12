#' External Normal Data for Propensity Score Matching
#'
#' This is a simulated dataset used to illustrate bayesian dynamic borrowing in
#' the normal case.
#'
#' @format ## `ex_norm_df`
#' A data frame with 150 rows and 6 columns:
#' \describe{
#'   \item{subjid}{Unique subject ID}
#'   \item{cov1}{Covariate 1, which is normally distributed around 50 with an sd of 10}
#'   \item{cov2}{Covariate 2, which is binary with about a 20% response rate}
#'   \item{cov3}{Covariate 3, which is binary with about a 60% response rate}
#'   \item{cov4}{Covariate 4, which is binary with about a 30% response rate}
#'   \item{y}{Response, which is normally distributed around 45 with an sd of 5}
#' }
"ex_norm_df"

#' Internal Normal Data for Propensity Score Matching
#'
#' This is a simulated dataset used to illustrate bayesian dynamic borrowing in
#' the normal case.
#'
#' @format ## `int_norm_df`
#' A data frame with 120 rows and 7 columns:
#' \describe{
#'   \item{subjid}{Unique subject ID}
#'   \item{cov1}{Covariate 1, which is normally distributed around 55 with an sd of 8}
#'   \item{cov2}{Covariate 2, which is binary with about a 30% response rate}
#'   \item{cov3}{Covariate 3, which is binary with about a 50% response rate}
#'   \item{cov4}{Covariate 4, which is binary with about a 30% response rate}
#'   \item{trt}{Treatment indicator, where 0 = control and 1 = active treatment}
#'   \item{y}{Response, which is normally distributed around 50 with an sd of 3}
#' }
"int_norm_df"

#' External Binary Data for Propensity Score Matching
#'
#' This is a simulated dataset used to illustrate bayesian dynamic borrowing in
#' the binary case.
#'
#' @format ## `ex_binary_df`
#' A data frame with 150 rows and 6 columns:
#' \describe{
#'   \item{subjid}{Unique subject ID}
#'   \item{cov1}{Covariate 1, which is normally distributed around 65 with an sd of 15}
#'   \item{cov2}{Covariate 2, which is binary with about a 30% response rate}
#'   \item{cov3}{Covariate 3, which is binary with about a 40% response rate}
#'   \item{cov4}{Covariate 4, which is binary with about a 50% response rate}
#'   \item{y}{Response, which is binary}
#' }
"ex_binary_df"

#' Internal Binary Data for Propensity Score Matching
#'
#' This is a simulated dataset used to illustrate bayesian dynamic borrowing in
#' the binary case.
#'
#' @format ## `int_binary_df`
#' A data frame with 150 rows and 6 columns:
#' \describe{
#'   \item{subjid}{Unique subject ID}
#'   \item{cov1}{Covariate 1, which is normally distributed around 62 with an sd of 5}
#'   \item{cov2}{Covariate 2, which is binary with about a 60% response rate}
#'   \item{cov3}{Covariate 3, which is binary with about a 40% response rate}
#'   \item{cov4}{Covariate 4, which is binary with about a 30% response rate}
#'   \item{trt}{Treatment indicator, where 0 = control and 1 = active treatment}
#'   \item{y}{Response, which is binary}
#' }
"int_binary_df"

