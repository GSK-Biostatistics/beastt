#' Simulate participant accrual times
#'
#' @param n Number of participants
#' @param accrual_periods Vector of right endpoints defining the time periods of
#'   accrual, e.g., c(6,8) defines 0<=x<6, 6<=x<8.
#' @param accrual_props Vector indicating the proportion of participants that are
#'   expected to be enrolled during each of the defined accrual periods. Must sum
#'   to 1.
#' @return Vector of accrual times corresponding to when each participant enters
#'   the study
#' @export
#' @importFrom stats runif
#'
#' @examples
#' at <- sim_accrual(n = 100000, accrual_periods = c(6, 8), accrual_props = c(.5, .5))
#' hist(at, breaks = 100, main = "Histogram of Enrollment Times", xlab = "Enrollment Time")
sim_accrual <- function(n, accrual_periods, accrual_props){
  if(length(accrual_periods) != length(accrual_props)){
    cli_abort("{.arg accrual_periods} and {.arg accrual_props} should have equal lengths")
  }

  # Calculate accrual probabilities for each period
  accrual_periods <- c(0, accrual_periods)

  # Randomly sample period during which each patient is enrolled
  periods <- sample(1:(length(accrual_periods) - 1), n, replace = TRUE, prob = accrual_props)

  # Randomly sample accrual time within specified interval for each patient
  t <- runif(n, min = accrual_periods[periods], accrual_periods[periods + 1])
  return(t)

}

#' Simulate event times for each individual from a piecewise constant hazard model
#'
#' @param n Number of individuals
#' @param hazard_periods Vector of break points between time periods with separate constant
#'        hazards, e.g., c(6,8) defines [0,6), [6,8), [8, infinity). Leave as NULL if
#'        defining only one hazard period.
#' @param hazard_values Vector of constant hazard values associated with the time intervals
#'
#' @return Vector of simulated times from the time-to-event distribution
#' @export
#' @importFrom stats rexp
#'
#' @examples
#' tte_dat <- sim_pw_const_haz(n = 100000, hazard_periods = c(6, 8), hazard_values = c(0.1, 0.1, 0.1))
#' hist(tte_dat, breaks = 100, main = "Event Time Distribution", xlab = "Event Time")
sim_pw_const_haz <- function(n, hazard_periods = NULL, hazard_values){
  if(length(hazard_periods)+1 != length(hazard_values)){
    cli_abort("{.arg hazard_values} should have length equal to one more than the length of {.arg hazard_periods}")
  } else if(any(c(hazard_periods, hazard_values) < 0)){
    cli_abort("{.arg hazard_values} and {.arg hazard_periods} can only have positive values")
  }

  # Define hazard periods, including 0 and (essentially) infinity
  hazard_periods <- c(0, hazard_periods, 1e50)

  # Simulate event times for all individuals using the hazard value for the first hazard period
  t <- rexp(n, rate = hazard_values[1])

  # Identify which individuals have initial simulated event times beyond the endpoint of the first hazard
  # period, and then sample an updated time using the hazard value for the second hazard period. Repeat
  # this process for all intervals to obtain the event times for all individuals.
  for(interval in 2:length(hazard_values)){
    r <- which(t > hazard_periods[interval])
    t[r] <- hazard_periods[interval] + rexp(length(r), rate = hazard_values[interval])
  }
  return(t)

}

#' Simulate event times for each participant from a Weibull proportional hazards
#' regression model
#'
#' @param weibull_ph_mod `survreg` object corresponding to a Weibull
#'   proportional hazards model fit using the external data
#' @param samp_df Data frame of covariates corresponding to the sample arm
#'   (control or treated) for which event times should be simulated. The column
#'   names should correspond to the covariate names in the `survreg` object.
#' @param cond_drift Optional value of the conditional drift by which the
#'   intercept in the Weibull proportional hazards regression model should be
#'   increased/decreased to incorporate the impact of unmeasurable sources of
#'   drift. Default is 0.
#' @param cond_trt_effect Optional value of the conditional treatment effect by
#'   which the intercept in the Weibull proportional hazards regression model
#'   should be increased/decreased if simulating event data for a treated arm.
#'   Default is 0.
#'
#' @return Vector of simulated event times from a Weibull proportional hazards
#'   regression model
#' @export
#' @importFrom stats predict rweibull
#'
#' @examples
#' library(dplyr)
#' library(survival)
#' # Model "true" regression coefficients and shape parameter using the external data
#' weibull_ph_mod <- survreg(Surv(y, event) ~ cov1 + cov2 + cov3 + cov4, data = ex_tte_df,
#'                           dist = "weibull")
#'
#' # Sample covariates for internal control arm via bootstrap from external data
#' samp_int_ctrl <- bootstrap_cov(ex_tte_df, n = 100) |>
#'   select(c(cov1, cov2, cov3, cov4))     # keep only covariate columns
#' tte_dat <- sim_weib_ph(weibull_ph_mod, samp_df = samp_int_ctrl)
sim_weib_ph <- function(weibull_ph_mod, samp_df, cond_drift = 0,
                                 cond_trt_effect = 0){
  if(length(cond_drift) > 1 | !is.numeric(cond_drift)){
    cli_abort("{.arg cond_drift} must be a single number")
  } else if(length(cond_trt_effect) > 1 | !is.numeric(cond_trt_effect)){
    cli_abort("{.arg cond_trt_effect} must be a single number")
  } else if(!inherits(weibull_ph_mod, "survreg") || weibull_ph_mod["dist"] != "weibull") {
    cli_abort("{.arg weibull_ph_mod} must be a `survreg` object with a Wiebull distribution")
  }
  # Calculate Weibull shape parameter
  shape <- 1 / weibull_ph_mod$scale

  # Calculate Weibull scale parameter for each individual using linear predictors
  lp_vec <- predict(weibull_ph_mod, newdata = samp_df, type = "lp")
  sigma <- exp(lp_vec - cond_drift - cond_trt_effect)

  # Sample event times
  t <- rweibull(n = nrow(samp_df), shape = shape, scale = sigma)
  t
}



#' Calculate Conditional Drift and Treatment Effect for Time-to-Event Outcome Models
#'
#' In order to properly generate time-to-event (TTE) outcome data for the
#' internal trial as part of a simulation study that investigates inverse
#' probability weighting, we need to translate the desired marginal drift and
#' treatment effect to the corresponding conditional drift and treatment effect
#' that can then be added into a TTE outcome model (e.g., Weibull proportional
#' hazards regression model) used to simulate response data.
#'
#' @param population A very large data frame (e.g., number of rows \eqn{\ge}
#'   100,000) where the columns correspond to the covariates defined in the
#'   `survreg` object for the Weibull proportional hazards model. This data
#'   frame should be constructed to represent the population of the internal
#'   trial according to the assumed covariate distributions (possibly imbalanced
#'   from the external data).
#' @param weibull_ph_mod `survreg` object corresponding to a Weibull
#'   proportional hazards model fit using the external data
#' @param marg_drift Vector of marginal drift values
#' @param marg_trt_eff Vector of marginal treatment effect values
#' @param analysis_time A single time point when survival probabilities will be
#'   calculated
#'
#' @details In simulation studies that investigate the properties of inverse
#'   probability weighted Bayesian dynamic borrowing, scenarios should be
#'   considered in which the underlying survival probabilities at some
#'   prespecified time \eqn{t} (`analysis_time`) for the internal and external
#'   control populations differ by varying amounts due to unmeasured confounding
#'   (i.e., drift, where positive values indicate a higher survival probability
#'   for the internal population). While values of drift and treatment effect
#'   (i.e., difference between the survival probabilities at time \eqn{t} for
#'   the treated and control populations) can be defined on the marginal scale
#'   for simulation studies, we must first convert these values to the
#'   conditional scale and then include these terms, along with covariates, in a
#'   Weibull proportional hazards (PH) regression outcome model when generating
#'   time-to-event (TTE) data for the internal arms. Doing so allows us to
#'   assume a relationship between the covariates and the response variable
#'   while properly accounting for drift and treatment effect.
#'
#'   To identify the conditional drift and treatment effect that correspond to
#'   specified values of marginal drift and treatment effect, we first bootstrap
#'   covariate vectors from the external data (e.g., \eqn{N \ge 100,000}) to
#'   construct a "population" that represents both the internal trial
#'   (possibly incorporating intentional covariate imbalance) and the external
#'   trial \emph{after} standardizing it to match the covariate distributions
#'   of the internal trial (allowing us to control for measured confounding
#'   from potential imbalance in the covariate distributions). Measured
#'   confounding can be incorporated into the data generation by bootstrapping
#'   a very large data frame (`population`) in which the distribution of at
#'   least one covariate is intentionally varied from that of the external data;
#'   additional \emph{unmeasured} drift can be incorporated through the
#'   translation of specified marginal values (`marg_drift`) to conditional
#'   values.
#'
#'   Let \eqn{\Delta} and \eqn{\delta} denote the marginal and conditional drift,
#'   respectively. For a specified value of \eqn{\Delta}, we can identify the
#'   corresponding \eqn{\delta} as the value that, when added as an additional
#'   term in the Weibull PH model survival function (i.e., additive change in
#'   the intercept) for each individual in the population, increases/decreases
#'   the population-averaged conditional probabilities of survival at time
#'   \eqn{t} by an amount approximately equal to \eqn{\Delta}. That is, the
#'   optimal \eqn{\delta} minimizes
#'
#'   \deqn{\left| \left( \frac{1}{N} \sum_{i=1}^N \exp \left( -\left\{ \exp
#'   \left( \boldsymbol{x}_i^\prime \boldsymbol{\beta}_{EC} + \delta \right)
#'   \times t \right\}^{\alpha_{EC}} \right) - \frac{1}{N} \sum_{i=1}^N \exp
#'   \left( -\left\{ \exp \left( \boldsymbol{x}_i^\prime \boldsymbol{\beta}_{EC}
#'   \right) \times t \right\}^{\alpha_{EC}} \right) \right) - \Delta \right|,}
#'
#'   where \eqn{\alpha_{EC}} is the Weibull shape parameter,
#'   \eqn{\boldsymbol{\beta}_{EC}} is a vector of regression coefficients, and
#'   \eqn{\boldsymbol{x}_i} is a vector of covariates (including an intercept
#'   term) from the bootstrapped population of size \eqn{N}. We note that
#'   \eqn{\alpha_{EC} = 1/\sigma_{EC}} and \eqn{\boldsymbol{\beta}_{EC} =
#'   -\boldsymbol{\xi}_{EC}} are calculated as functions of the scale parameter
#'   (\eqn{\sigma_{EC}}) and coefficients (\eqn{\boldsymbol{\xi}_{EC}})
#'   estimated by the `survreg` object that was fit to the external data, and we
#'   assume here that these estimates are the "true" shape and covariate effects
#'   when generating response data. In the formula above, the first and second
#'   terms correspond to the population-averaged conditional survival functions
#'   (i.e., the marginal survival probabilities) at time \eqn{t} for the
#'   internal control population with drift and the external control population
#'   (with covariate distributions standardized to match the internal trial),
#'   respectively.
#'
#'   If we now denote the marginal and conditional treatment effect by
#'   \eqn{\Gamma} and \eqn{\gamma}, respectively, we can use a similar process
#'   to identify the optimal \eqn{\gamma} that approximately corresponds to the
#'   specified value of \eqn{\Gamma}, which is done by minimizing the following:
#'
#'   \deqn{\left| \left( \frac{1}{N} \sum_{i=1}^N \exp \left( -\left\{ \exp
#'   \left( \boldsymbol{x}_i^\prime \boldsymbol{\beta}_{EC} + \delta + \gamma
#'   \right) \times t \right\}^{\alpha_{EC}} \right) - \frac{1}{N} \sum_{i=1}^N
#'   \exp \left( -\left\{ \exp \left( \boldsymbol{x}_i^\prime
#'   \boldsymbol{\beta}_{EC} + \delta \right) \times t \right\}^{\alpha_{EC}}
#'   \right) \right) - \Gamma \right|,}
#'
#'   where the first term is the average of the conditional survival functions
#'   (i.e., the marginal survival probabilities) at time \eqn{t} for the
#'   internal treated population.
#'
#'   See [here](https://github.com/GSK-Biostatistics/beastt/blob/e2b41fe90f639924d10c0d94ceff04a74d0ce617/inst/templates/tte-template.R)
#'   for a survival simulation example.
#'
#' @returns tibble of all combinations of the marginal drift and treatment
#'   effect. For each row the conditional drift and treatment effect has been
#'   calculated as well as the true marginal survival probabilities at time `t`
#'   for the control and treatment populations.
#' @export
#'
#' @examples
#' library(dplyr)
#' library(survival)
#' # Model "true" regression coefficients using the external data
#' weibull_ph_mod <- survreg(Surv(y, event) ~ cov1 + cov2 + cov3 + cov4, data = ex_tte_df,
#'                           dist = "weibull")
#'
#' # Bootstrap internal control "population" with imbalance w.r.t. covariate 2
#' pop_int_ctrl <- bootstrap_cov(ex_tte_df, n = 100000, imbal_var = cov2,
#'                               imbal_prop = 0.25, ref_val = 0) |>
#'   select(c(cov1, cov2, cov3, cov4))     # keep only covariate columns
#'
#' # Convert the marginal drift and treatment effects to conditional
#' calc_cond_weibull(population = pop_int_ctrl, weibull_ph_mod,
#'                   marg_drift = c(-.1, 0, .1), marg_trt_eff = c(0, .10),
#'                   analysis_time = 12)
#'
#' @importFrom dplyr left_join ungroup
#' @importFrom purrr map2_dbl
calc_cond_weibull <- function(population, weibull_ph_mod, marg_drift, marg_trt_eff,
                              analysis_time){
  if(!inherits(weibull_ph_mod, "survreg")){
    cli_abort("{.arg weibull_ph_mod} must be a survreg object")
  } else if(weibull_ph_mod$dist != "weibull"){
    cli_abort("{.arg weibull_ph_mod} must use a weibull distribution")
  }

  if(length(analysis_time) > 1 | !is.numeric(analysis_time)){
    cli_abort("{.arg analysis_time} must be a single number")
  }

  if(!inherits(population, "data.frame")){
    cli_abort("{.arg population} must be a tibble or dataframe. If you are using lists, check you haven't converted the dataframe into a list of vectors")
  }

  cov_vec <-   names(weibull_ph_mod$coefficients) |>
    discard(\(x) x == "(Intercept)")
  beta_coefs <- weibull_ph_mod$coefficients * -1   # transform coefficients to match desired parameterization
  shape <- 1 / weibull_ph_mod$scale                # transform Weibull shape to match desired parameterization

  if(all(cov_vec %in% colnames(population))){

    # Construct design matrix (with intercept) for the large sample ("population")
    # corresponding to the internal control (IC) population
    X_IC = as.matrix(cbind(int = 1, select(population,cov_vec)))

    # Marginal survival probability (SP) at the analysis time for the external control
    # (EC) population AFTER standardizing it to match the covariate distributions
    # of the IC population (possibly imbalanced)
    EC_SP_star <- mean( exp( -(exp(X_IC %*% beta_coefs) * analysis_time)^shape ) )

    # Remove scenarios for which the EC surv probability + marginal drift + marginal trt effect
    # is outside the range (0,1)
    scenarios <- crossing(
      marg_drift,marg_trt_eff
    ) |>
      filter(marg_drift + marg_trt_eff + EC_SP_star > 0 &
               marg_drift + marg_trt_eff + EC_SP_star < 1)

    # Identify the optimal conditional drift value ("delta") that corresponds to
    # the defined marginal drift value ("Delta"). If we calculate the marginal
    # survival probability at the analysis time for the IC population by averaging
    # over each individual's conditional probability of survival (conditional on
    # their covariates, the shape parameter, and delta), the optimal value of the
    # conditional drift is the delta that results in a difference between the IC
    # survival probability and the EC survival probability (after standardization)
    # that is approximately equal to the marginal drift Delta.
    delta_df <- scenarios |>
      pull(marg_drift) |>
      unique() |>
      map(function(Delta){
        if(Delta == 0){
          delta_val <- 0    # set delta = 0 if Delta = 0
        } else {
          optim_delta <- function(x){
            IC_SP <- mean( exp( -(exp(X_IC %*% beta_coefs + x) *
                                    analysis_time)^shape ) )   # IC surv prob
            abs( Delta - ( IC_SP - EC_SP_star ) )
          }
          delta_val <- optimize( f = optim_delta, lower = -5, upper = 5 )$minimum
        }

        # Calculate the "true" IC survival probability at the analysis time for
        # each defined value of drift
        IC_SP_true_val <- mean( exp( -(exp(X_IC %*% beta_coefs + delta_val) *
                                         analysis_time)^shape ) )
        c("marg_drift" = Delta, "conditional_drift" = delta_val,
          "true_control_surv_prob" = IC_SP_true_val)
      }) |>
      bind_rows()

    # Identify the optimal conditional trt effect value ("gamma") that corresponds
    # to the defined marginal trt effect value ("Gamma"). If we calculate the
    # marginal survival probability at the analysis time for the internal treated
    # (IT) population by averaging over each individual's conditional probability
    # of survival (conditional on their covariates, the shape parameter, delta,
    # and gamma), the optimal value of the conditional trt effect is the gamma
    # that results in a difference between the IT survival probability and the
    # IC survival probability that is approximately equal to the marginal trt
    # effect (Gamma).

    cond_df <- scenarios |>
      left_join(delta_df, by = "marg_drift") |>
      mutate(conditional_trt_eff =
               map2_dbl(marg_trt_eff, .data$conditional_drift,
                        function(Gamma, delta){
                          if(Gamma == 0){
                            gamma_val <- 0    # set gamma = 0 if Gamma = 0
                          } else {
                            optim_gamma <- function(x){
                              IT_SP <- mean( exp( -(exp(X_IC %*% beta_coefs + delta + x) *
                                                      analysis_time)^shape ) )   # IT surv prob
                              abs( Gamma - ( IT_SP - mean( exp( -(exp(X_IC %*% beta_coefs + delta) *
                                                                    analysis_time)^shape ) ) ) )
                            }
                            gamma_val <- optimize( f = optim_gamma, lower = -5, upper = 5 )$minimum
                          }
                          gamma_val
                        })) |>
      rowwise() |>
      # Calculate the "true" IT survival probability at the analysis time for each defined value
      # of drift and treatment effect
      mutate(true_trt_surv_prob = mean( exp( -(exp(X_IC %*% beta_coefs + .data$conditional_drift +
                                                     .data$conditional_trt_eff) * analysis_time)^shape ) ))|>
      ungroup()
  } else {
    cli_abort("Not all covariates in {.arg weibull_ph_mod} are in the population")
  }
  cond_df

}

#' Calculate the analysis time based on a target number of events
#'
#' @param study_time Vector of study (accrual + observed) times
#' @param observed_time Vector of observed times
#' @param event_indicator Vector of boolean values (True/False or 1/0) indicating if the observed time value is an event or censoring
#' @param target_events Number of target events, if only using target follow-up time leave NULL
#' @param target_follow_up Target follow-up for each subject, if only using target events leave NULL
#'
#' @returns Time of analysis
#' @export
#'
#' @examples
#' library(dplyr)
#' # Determining analysis time by reaching a number of events
#' ex_tte_df |> mutate(
#'   analysis_time = calc_study_duration(study_time = total_time, observed_time = y,
#'                                      event_indicator = event, target_events = 30)
#' )
#' # Determining analysis time by minimum follow-up time
#' ex_tte_df |> mutate(
#'   analysis_time = calc_study_duration(study_time = total_time, observed_time = y,
#'                                      event_indicator = event, target_follow_up = 12)
#' )
#' # Or use both and whichever happens first
#' ex_tte_df |> mutate(
#'   analysis_time = calc_study_duration(study_time = total_time, observed_time = y,
#'                                      event_indicator = event,
#'                                      target_events = 30, target_follow_up = 12)
#' )
calc_study_duration <- function(study_time, observed_time, event_indicator,
                               target_events = NULL, target_follow_up = NULL){
  if(is.null(target_events) & is.null(target_follow_up)){
    cli_abort("Either {.arg target_events} or {.arg target_follow_up} need to be not NULL")
  }
  analy_time <- max(study_time) # Maximum time
  event_indicator <- as.logical(event_indicator)
  if(!is.null(target_events)){
    #The study time where the number of events equals the target (order and then just get the ith event)
    event_st <- study_time[event_indicator]
    event_time <- event_st[order(event_st)][target_events]
    analy_time <- if_else(event_time < analy_time, event_time, analy_time,
                          missing = analy_time)
  }
  if(!is.null(target_follow_up)){
    # Filter out any censored
    event_st <- study_time[event_indicator]
    event_ot <- observed_time[event_indicator]
    subj_to_make_fu <- which(event_ot >= target_follow_up) #Get all individuals who have observed times at least to the minimum follow-up
    accrual_time = event_st-event_ot
    min_fu_time <- max(accrual_time[subj_to_make_fu])+target_follow_up
    analy_time <- if_else(min_fu_time < analy_time, min_fu_time, analy_time,
                          missing = analy_time)
  }

  analy_time

}
