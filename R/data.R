#' External Normal Control Data for Propensity Score Balancing
#'
#' This is a simulated dataset used to illustrate Bayesian dynamic borrowing in
#' the case when borrowing from an external control arm with a normal endpoint,
#' where the baseline covariate distributions of the internal and external data
#' are balanced via inverse probability weighting.
#'
#' @format ## `ex_norm_df`
#' A data frame with 150 rows and 6 columns:
#' \describe{
#'   \item{subjid}{Unique subject ID}
#'   \item{cov1}{Covariate 1, which is normally distributed around 50 with a SD of 10}
#'   \item{cov2}{Covariate 2, which is binary (0 vs. 1) with about 20% of participants having level 1}
#'   \item{cov3}{Covariate 3, which is binary (0 vs. 1) with about 60% of participants having level 1}
#'   \item{cov4}{Covariate 4, which is binary (0 vs. 1) with about 30% of participants having level 1}
#'   \item{y}{Response, which is normally distributed with a SD of 0.15}
#' }
"ex_norm_df"

#' Internal Normal Data for Propensity Score Balancing
#'
#' This is a simulated dataset used to illustrate Bayesian dynamic borrowing in
#' the case when borrowing from an external control arm with a normal endpoint,
#' where the baseline covariate distributions of the internal and external data
#' are balanced via inverse probability weighting.
#'
#' @format ## `int_norm_df`
#' A data frame with 120 rows and 7 columns:
#' \describe{
#'   \item{subjid}{Unique subject ID}
#'   \item{cov1}{Covariate 1, which is normally distributed around 55 with a SD of 8}
#'   \item{cov2}{Covariate 2, which is binary (0 vs. 1) with about 30% of participants having level 1}
#'   \item{cov3}{Covariate 3, which is binary (0 vs. 1) with about 50% of participants having level 1}
#'   \item{cov4}{Covariate 4, which is binary (0 vs. 1) with about 30% of participants having level 1}
#'   \item{trt}{Treatment indicator, where 0 = control and 1 = active treatment}
#'   \item{y}{Response, which is normally distributed with a SD of 0.15}
#' }
"int_norm_df"

#' External Binary Control Data for Propensity Score Balancing
#'
#' This is a simulated dataset used to illustrate Bayesian dynamic borrowing in
#' the case when borrowing from an external control arm with a binary endpoint,
#' where the baseline covariate distributions of the internal and external data
#' are balanced via inverse probability weighting.
#'
#' @format ## `ex_binary_df`
#' A data frame with 150 rows and 6 columns:
#' \describe{
#'   \item{subjid}{Unique subject ID}
#'   \item{cov1}{Covariate 1, which is normally distributed around 65 with a SD of 10}
#'   \item{cov2}{Covariate 2, which is binary (0 vs. 1) with about 30% of participants having level 1}
#'   \item{cov3}{Covariate 3, which is binary (0 vs. 1) with about 40% of participants having level 1}
#'   \item{cov4}{Covariate 4, which is binary (0 vs. 1) with about 50% of participants having level 1}
#'   \item{y}{Response, which is binary (0 vs. 1)}
#' }
"ex_binary_df"

#' Internal Binary Data for Propensity Score Balancing
#'
#' This is a simulated dataset used to illustrate Bayesian dynamic borrowing in
#' the case when borrowing from an external control arm with a binary endpoint,
#' where the baseline covariate distributions of the internal and external data
#' are balanced via inverse probability weighting.
#'
#' @format ## `int_binary_df`
#' A data frame with 160 rows and 7 columns:
#' \describe{
#'   \item{subjid}{Unique subject ID}
#'   \item{cov1}{Covariate 1, which is normally distributed around 62 with an sd of 8}
#'   \item{cov2}{Covariate 2, which is binary (0 vs. 1) with about 40% of participants having level 1}
#'   \item{cov3}{Covariate 3, which is binary (0 vs. 1) with about 40% of participants having level 1}
#'   \item{cov4}{Covariate 4, which is binary (0 vs. 1) with about 60% of participants having level 1}
#'   \item{trt}{Treatment indicator, where 0 = control and 1 = active treatment}
#'   \item{y}{Response, which is binary (0 vs. 1)}
#' }
"int_binary_df"

#' External Time-to-Event Control Data for Propensity Score Balancing
#'
#' This is a simulated dataset used to illustrate Bayesian dynamic borrowing in
#' the case when borrowing from an external control arm with a time-to-event endpoint,
#' where the baseline covariate distributions of the internal and external data
#' are balanced via inverse probability weighting.
#'
#' @format ## `ex_tte_df`
#' A data frame with 150 rows and 9 columns:
#' \describe{
#'   \item{subjid}{Unique subject ID}
#'   \item{y}{Response (observed event/censored time)}
#'   \item{enr_time}{Enrollment time}
#'   \item{total_time}{Time from study start}
#'   \item{event}{Response, which is binary (0 vs. 1)}
#'   \item{cov1}{Covariate 1, which is normally distributed around 65 with a SD of 10}
#'   \item{cov2}{Covariate 2, which is binary (0 vs. 1) with about 30% of participants having level 1}
#'   \item{cov3}{Covariate 3, which is binary (0 vs. 1) with about 40% of participants having level 1}
#'   \item{cov4}{Covariate 4, which is binary (0 vs. 1) with about 50% of participants having level 1}
#' }
"ex_tte_df"

#' Internal Time to Event Control Data for Propensity Score Balancing
#'
#' This is a simulated dataset used to illustrate Bayesian dynamic borrowing in
#' the case when borrowing from an external control arm with a time to event endpoint,
#' where the baseline covariate distributions of the internal and external data
#' are balanced via inverse probability weighting.
#'
#' @format ## `int_tte_df`
#' A data frame with 150 rows and 6 columns:
#' \describe{
#'   \item{subjid}{Unique subject ID}
#'   \item{y}{Response, which is binary (0 vs. 1)}
#'   \item{enr_time}{Enrollment time}
#'   \item{total_time}{Time from study start}
#'   \item{event}{Response, which is binary (0 vs. 1)}
#'   \item{trt}{Treatment indicator, where 0 = control and 1 = active treatment}
#'   \item{cov1}{Covariate 1, which is normally distributed around 65 with a SD of 10}
#'   \item{cov2}{Covariate 2, which is binary (0 vs. 1) with about 30% of participants having level 1}
#'   \item{cov3}{Covariate 3, which is binary (0 vs. 1) with about 40% of participants having level 1}
#'   \item{cov4}{Covariate 4, which is binary (0 vs. 1) with about 50% of participants having level 1}
#' }
"int_tte_df"
