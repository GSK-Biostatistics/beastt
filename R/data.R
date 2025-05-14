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
#'   \item{y}{Response (observed time at which the participant either had an event or was censored)}
#'   \item{enr_time}{Enrollment time}
#'   \item{total_time}{Time from study start}
#'   \item{event}{Event indicator (1: event; 0: censored)}
#'   \item{cov1}{Covariate 1, which is normally distributed around 65 with a SD of 10}
#'   \item{cov2}{Covariate 2, which is binary (0 vs. 1) with about 30% of participants having level 1}
#'   \item{cov3}{Covariate 3, which is binary (0 vs. 1) with about 40% of participants having level 1}
#'   \item{cov4}{Covariate 4, which is binary (0 vs. 1) with about 50% of participants having level 1}
#' }
"ex_tte_df"

#' Internal Time-to-Event Control Data for Propensity Score Balancing
#'
#' This is a simulated dataset used to illustrate Bayesian dynamic borrowing in
#' the case when borrowing from an external control arm with a time-to-event endpoint,
#' where the baseline covariate distributions of the internal and external data
#' are balanced via inverse probability weighting.
#'
#' @format ## `int_tte_df`
#' A data frame with 160 rows and 10 columns:
#' \describe{
#'   \item{subjid}{Unique subject ID}
#'   \item{y}{Response (observed time at which the participant either had an event or was censored)}
#'   \item{enr_time}{Enrollment time}
#'   \item{total_time}{Time from study start}
#'   \item{event}{Event indicator (1: event; 0: censored)}
#'   \item{trt}{Treatment indicator, where 0 = control and 1 = active treatment}
#'   \item{cov1}{Covariate 1, which is normally distributed around 62 with a SD of 8}
#'   \item{cov2}{Covariate 2, which is binary (0 vs. 1) with about 40% of participants having level 1}
#'   \item{cov3}{Covariate 3, which is binary (0 vs. 1) with about 40% of participants having level 1}
#'   \item{cov4}{Covariate 4, which is binary (0 vs. 1) with about 60% of participants having level 1}
#' }
"int_tte_df"

#' Binary Simulation Data
#'
#' This is an example of output from a simulation study that investigates the
#' operating characteristics of inverse probability weighted Bayesian dynamic
#' borrowing for the case with a binary outcome. This output was generated
#' based on the binary simulation template. For this simulation study, only the
#' degree of covariate imbalance (as indicated by `population`) and the
#' marginal treatment effect were varied.
#'
#' @format ## `binary_sim_df` A data frame with 255 rows and 6 columns:
#' \describe{
#'   \item{population}{Populations defined by different covariate imbalances}
#'   \item{marg_trt_eff}{Marginal treatment effect}
#'   \item{true_control_RR}{True control response rate on the marginal scale}
#'   \item{reject_H0_yes}{Probability of rejecting the null hypothesis in the case with borrowing}
#'   \item{no_borrowing_reject_H0_yes}{Probability of rejecting the null hypothesis without borrowing}
#'   \item{pwr_prior}{Vector of power priors (or some other informative
#'   prior distribution for the control marginal parameter of interest
#'   based on the external data) as distributional objects}
#' }
"binary_sim_df"


#' Time to Event Simulation Data
#'
#' This is an example of output from a simulation study that investigates the
#' operating characteristics of inverse probability weighted Bayesian dynamic
#' borrowing for the case with a time-to-event outcome. This output was generated
#' based on the time-to-event simulation template. For this simulation study, only the
#' degree of covariate imbalance (as indicated by `population`) and the
#' marginal treatment effect were varied.
#'
#' @format ## `tte_sim_df` A data frame with 18 rows and 7 columns:
#' \describe{
#'   \item{population}{Populations defined by different covariate imbalances}
#'   \item{marg_trt_eff}{Marginal treatment effect}
#'   \item{true_control_surv_prop}{True control survival probability on the marginal scale}
#'   \item{reject_H0_yes}{Probability of rejecting the null hypothesis in the case with borrowing}
#'   \item{no_borrowing_reject_H0_yes}{Probability of rejecting the null hypothesis without borrowing}
#'   \item{pwr_prior}{Vector of IPW power priors as distributional objects}
#'   \item{mix_prior}{Vector of mixture priors, the robustified IPW power priors, as distributional objects}
#' }
"tte_sim_df"
