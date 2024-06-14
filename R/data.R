#' External Normal Control Data for Propensity Score Matching
#'
#' This is a simulated dataset used to illustrate Bayesian dynamic borrowing in
#' the case when borrowing from an external control arm with a normal endpoint.
#'
#' @format ## `ex_norm_df`
#' A data frame with 150 rows and 6 columns:
#' \describe{
#'   \item{subjid}{Unique subject ID}
#'   \item{cov1}{Covariate 1, which is normally distributed around 50 with an sd of 10}
#'   \item{cov2}{Covariate 2, which is binary (0 vs. 1) with about 20% of participants having level 1}
#'   \item{cov3}{Covariate 3, which is binary (0 vs. 1) with about 60% of participants having level 1}
#'   \item{cov4}{Covariate 4, which is binary (0 vs. 1) with about 30% of participants having level 1}
#'   \item{y}{Response, which is normally distributed around 45 with an sd of 5}
#' }
"ex_norm_df"

#' Internal Normal Data for Propensity Score Matching
#'
#' This is a simulated dataset used to illustrate Bayesian dynamic borrowing in
#' the case when borrowing from an external control arm with a normal endpoint.
#'
#' @format ## `int_norm_df`
#' A data frame with 120 rows and 7 columns:
#' \describe{
#'   \item{subjid}{Unique subject ID}
#'   \item{cov1}{Covariate 1, which is normally distributed around 55 with an sd of 8}
#'   \item{cov2}{Covariate 2, which is binary (0 vs. 1) with about 30% of participants having level 1}
#'   \item{cov3}{Covariate 3, which is binary (0 vs. 1) with about 50% of participants having level 1}
#'   \item{cov4}{Covariate 4, which is binary (0 vs. 1) with about 30% of participants having level 1}
#'   \item{trt}{Treatment indicator, where 0 = control and 1 = active treatment}
#'   \item{y}{Response, which is normally distributed around 50 with an sd of 3}
#' }
"int_norm_df"

#' External Binary Control Data for Propensity Score Matching
#'
#' This is a simulated dataset used to illustrate Bayesian dynamic borrowing in
#' the case when borrowing from an external control arm with a binary endpoint.
#'
#' @format ## `ex_binary_df`
#' A data frame with 150 rows and 6 columns:
#' \describe{
#'   \item{subjid}{Unique subject ID}
#'   \item{cov1}{Covariate 1, which is normally distributed around 65 with an sd of 15}
#'   \item{cov2}{Covariate 2, which is binary (0 vs. 1) with about 30% of participants having level 1}
#'   \item{cov3}{Covariate 3, which is binary (0 vs. 1) with about 40% of participants having level 1}
#'   \item{cov4}{Covariate 4, which is binary (0 vs. 1) with about 50% of participants having level 1}
#'   \item{y}{Response, which is binary (0 vs. 1) with a response rate of 58%}
#' }
"ex_binary_df"

#' Internal Binary Data for Propensity Score Matching
#'
#' This is a simulated dataset used to illustrate Bayesian dynamic borrowing in
#' the case when borrowing from an external control arm with a binary endpoint.
#'
#' @format ## `int_binary_df`
#' A data frame with 150 rows and 6 columns:
#' \describe{
#'   \item{subjid}{Unique subject ID}
#'   \item{cov1}{Covariate 1, which is normally distributed around 62 with an sd of 5}
#'   \item{cov2}{Covariate 2, which is binary (0 vs. 1) with about 60% of participants having level 1}
#'   \item{cov3}{Covariate 3, which is binary (0 vs. 1) with about 40% of participants having level 1}
#'   \item{cov4}{Covariate 4, which is binary (0 vs. 1) with about 30% of participants having level 1}
#'   \item{trt}{Treatment indicator, where 0 = control and 1 = active treatment}
#'   \item{y}{Response, which is binary (0 vs. 1)}
#' }
"int_binary_df"

