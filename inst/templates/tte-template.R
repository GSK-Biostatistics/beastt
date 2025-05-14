# Simulation Overview ---------------------------------------------------------

# This script allows users to conduct a simulation study using the IPW BDB
# methodology for a time-to-event outcome as outlined in Psioda et al. (2025):
# https://doi.org/10.1080/10543406.2025.2489285
# Additional details relating to the methodology and the recommended
# simulation process are discussed in the paper.
# This script shows what each step would look like in the app while keeping
# the example runnable. We encourage users to read each commented step of the
# simulation study carefully and adjust the code as needed.

# Load packages
library(tidyverse)
library(beastt)
library(distributional)
library(furrr)
library(survival)

# Set up parallelization where "workers" indicates the number of cores
set.seed(123)
num_cores <- availableCores() - 1     # number of cores
plan(multisession, workers = num_cores)
para_opts <- furrr_options(seed = TRUE)

# Simulation Setup -------------------------------------------------------------

# Step 1: Read in external control data
external_dat <- beastt::ex_tte_df
# external_dat <- read_csv("Location of the external data")

# Model "true" regression coefficients corresponding to intercept and covariate
# effects using a Weibull proportional (PH) hazards regression model with the
# external data
weibull_ph_mod <- survreg(Surv(y, event) ~ cov1 + cov2 + cov3 + cov4, data = external_dat, dist = "weibull")
# weibull_ph_mod <- survreg(YOUR MODEL HERE, data = external_dat, dist = "weibull")

# Step 2: Define the underlying population characteristics to vary. Here, drift
# is the difference in the marginal control survival probabilities at a
# prespecified time t between the external and internal trials attributed to
# unmeasured confounding, and treatment effect is the difference in marginal
# survival probabilities at time t between the treated and control populations.
drift_surv_prob <- seq(-0.15, 0.15, by = 0.05)    # Internal vs External (positive drift = internal probability is higher)
trt_effect_surv_prob <- c(0, .15)                 # Treatment vs Control (positive TE = treatment probability is higher)

# Step 3: Convert the drift and treatment effects from marginal to conditional models

# Later, we will generate response data for the internal arms using a conditional
# outcome model (i.e., Weibull PH regression) that assumes a relationship between
# the covariates and the outcome. To account for the specified drift and treatment
# effect, we first need to convert these effects from the marginal scale to the
# conditional scale. For more information about this process, please refer to the
# function details for `calc_cond_weibull`.

# 3 a) Bootstrap populations corresponding to the internal trial under various
# scenarios (e.g., same covariate distributions as the external data, imbalanced
# covariate distributions compared to the external data). These scenario-specific
# populations will be used to identify the conditional drift and treatment effects
# and to later sample the covariate vectors for the internal trial arms.
pop_size <- 100000     # set "population" size to be very large (e.g., >=100,000)
ex_dat_cov <- external_dat |>
  select(cov1, cov2, cov3, cov4)     # keep only covariate columns

# Generate a population without covariate imbalance
no_imbal_pop <- bootstrap_cov(ex_dat_cov, n = pop_size)

# Generate populations that incorporate covariate imbalance with respect to a single
# binary covariate ("imbal_var"). Define the degree of imbalanced in the distribution
# by specifying the proportion of individuals with the reference level ("imbal_prop").
imbal_pop <- bootstrap_cov(ex_dat_cov, n = pop_size, imbal_var = cov2,
                           imbal_prop = 0.25, ref_val = 0)

# Combine all populations into a list
pop_ls <- list(no_imbal_pop, imbal_pop)
# Naming the populations will make them easier to identify later
names(pop_ls) <- c("no imbalance", "cov2: 0.25")

# 3 b) For each population, identify the conditional drift and treatment effect
# values that best correspond to the specified values of marginal drift and
# treatment effect
surv_prob_time <- 12      # time when marginal survival probabilities should be calculated (t)
pop_var <- pop_ls |>
  future_map(calc_cond_weibull, weibull_ph_mod, drift_surv_prob,
             trt_effect_surv_prob, surv_prob_time) |> # returns a list of data frames
  bind_rows(.id = "population")  # combines list into 1 data frame with a population column

# Step 4: Create a data frame of all possible simulation scenarios using your
# population variables (e.g., "drift_surv_prob", "trt_effect_surv_prob",
# "imbal_var", "imbal_prop") as a foundation. Add any additional design variables
# you would like to vary. If you would like to compare priors, you can make a
# vector of distributional objects. Simulation scenarios are defined by unique
# combinations of the population variables and design variables. For inputs that
# require a vector of information (e.g., accrual_periods), put the vector(s) in
# a list.
all_sims <- pop_var |>
  crossing(
    # Internal control sample size
    n_int_cont = 65,

    # Internal treatment sample size
    n_int_trt = 130,

    # Right time points defining accrual intervals
    # Include in a list if input is a vector (i.e., >1 accrual period defined)
    accrual_periods = list(c(0.5, 1.0, 2.0)),

    # Proportion of participants expected to be enrolled during each accrual period
    # Include in a list if input is a vector (i.e., >1 accrual period defined)
    accrual_props = list(c(0.1, 0.25, 0.65)),

    # Time-to-censorship hazard values assuming piecewise constant hazard model
    # Include in a list if input is a vector (i.e., >1 censorship periods)
    cns_hazard_values  = list(c(0.004, 0.005)),

    # Time-to-censorship hazard break points for piecewise constant hazard model
    # Include in a list if input is a vector (i.e., >1 break points)
    cns_hazard_periods = c(6),

    # Target number of events when analysis should be triggered. Leave commented
    # out if analysis time is determined only by target follow-up time.
    # target_events = NULL,

    # Target follow-up time when analysis should be triggered. Analysis is
    # conducted when all participants in the risk set reach the target follow-up
    # time OR when the target number of events is reached (if applicable),
    # whichever occurs first.
    target_follow_time = 12,

    # Initial prior for the intercept (log of the inverse scale hyperparameter)
    # in the Weibull PH regression model, incorporated into the power prior.
    # Must be a normal distributional object(s).
    initial_prior_intercept = dist_normal(0, 10),

    # Scale of the half-normal initial prior for the Weibull shape hyperparameter
    # incorporated into the power prior
    initial_shape_hyperpar = 50,

    # Prior mixture weight associated with the informative component (i.e.,
    # IPW power prior) in the robust mixture prior
    mix_weight = 0.5
  ) |>
  mutate(scenario = row_number()) |>    # add a scenario ID
  crossing(iter_id = c(1:1000))         # add an iteration ID (within scenario)

# Simulations ------------------------------------------------------------------

# We now iterate over all rows in the simulation data frame and calculate
# operating characteristics for each scenario. The pmap and list functions make
# it possible to refer to each column of the data frame by its name. To step
# through this code, add browser().
sim_output <- all_sims |>
  future_pmap(function(...){
    output <- with(list(...), {

      # Simulate data ----------------------------------------------------------
      # Sample covariates from the scenario-specific population for the internal arms
      int_cont_cov_df <- slice_sample(pop_ls[[population]], n = n_int_cont)  # control
      int_trt_cov_df <- slice_sample(pop_ls[[population]], n = n_int_trt)    # treatment

      # For each participant, simulate the following:
      #   (1) accrual time (uniform within each accrual period),
      #   (2) theoretical event time (Weibull PH regression, using the Weibull
      #       PH model fit previously to the external data while incorporating
      #       conditional drift and treatment effect),
      #   (3) theoretical censoring time (piecewise constant hazard model).
      # Calculate the observed time as the minimum of the theoretical event and
      # censoring times, and an event indicator (1 if obs. time is event time).
      int_cont_df <- int_cont_cov_df |>
        mutate(
          subjid = row_number(),
          accrual_time = sim_accrual(n = n_int_cont,
                                     accrual_periods = accrual_periods,
                                     accrual_props = accrual_props),
          sim_event_time = sim_weib_ph(weibull_ph_mod, samp_df = int_cont_cov_df,
                                       cond_drift = conditional_drift),
          sim_censor_time = sim_pw_const_haz(n = n_int_cont,
                                             hazard_periods = cns_hazard_periods,
                                             hazard_values = cns_hazard_values),
          obs_time = pmin(sim_event_time, sim_censor_time),
          event_ind = obs_time == sim_event_time,
          total_time = accrual_time + obs_time,
          trt = 0
        )

      # For the treatment arm, do the same as above, but add both conditional
      # drift and conditional treatment effects when simulating the theoretical
      # event times
      int_trt_df <- int_trt_cov_df |>
        mutate(
          subjid = row_number(),
          accrual_time = sim_accrual(n = n_int_trt,
                                     accrual_periods = accrual_periods,
                                     accrual_props = accrual_props),
          sim_event_time = sim_weib_ph(weibull_ph_mod, samp_df = int_trt_cov_df,
                                       cond_drift = conditional_drift,
                                       cond_trt_effect = conditional_trt_eff),
          sim_censor_time = sim_pw_const_haz(n = n_int_trt,
                                             hazard_periods = cns_hazard_periods,
                                             hazard_values = cns_hazard_values),
          obs_time = pmin(sim_event_time, sim_censor_time),
          event_ind = obs_time == sim_event_time,
          total_time = accrual_time + obs_time,
          trt = 1
        )

      # Combine data for internal control and treatment arms
      int_df <- bind_rows(int_cont_df, int_trt_df)

      # Determine when the analysis time will be based on the target number of
      # events and/or target follow-up time. Then administratively censor
      # participants in the risk set with observed event/censoring times after
      # the analysis time. Remove participants with accrual times after the
      # analysis time. Change the names of the response and event variables in
      # the internal data frame to match the names in the external dataset.
      int_df_admin_cen <- int_df |>
        mutate(
          analysis_time = calc_study_duration(
            study_time = total_time,
            observed_time = obs_time,
            event_indicator = event_ind,
            # target_events = target_events,   # uncomment if using target number of events
            target_follow_up = target_follow_time),
          event_ind = as.numeric(total_time <= analysis_time & event_ind == 1),
          total_time_admin_cen = if_else(total_time > analysis_time,
                                         analysis_time, total_time),
          obs_time = total_time_admin_cen - accrual_time
        ) |>
        filter(accrual_time < analysis_time) |>  # remove if enrolled after analysis time
        select(subjid,
               y = obs_time,       # renaming to match name in external data
               event = event_ind,
               trt, starts_with("cov"))

      # Analysis ---------------------------------------------------------------
      # Calculate the propensity scores and inverse probability weights for all control
      # participants (external and internal)
      ps_obj <- calc_prop_scr(internal_df = filter(int_df_admin_cen, trt == 0),
                              external_df = external_dat,
                              id_col = subjid,
                              # model = #YOUR MODEL HERE
                              model = ~ cov1 + cov2 + cov3 + cov4)

      # Approximate the inverse probability weighted (IPW) power prior for the log(shape)
      # and intercept of the Weibull PH regression model as a bivariate normal distribution
      pwr_prior <- calc_power_prior_weibull(ps_obj,
                                            response = y,
                                            event = event,
                                            intercept = initial_prior_intercept,
                                            shape = initial_shape_hyperpar,
                                            approximation = "Laplace")

      # "Robustify" the power prior by mixing it with the a vague prior to create an
      # inverse probability weighted robust mixture prior (IPW RMP). Weight the IPW
      # power prior (i.e., the informative component) using the specified mixture weight
      r_external <- sum(external_dat$event)   # number of observed events in external control data
      mix_prior <- robustify_mvnorm(pwr_prior, r_external,
                                    weights = c(mix_weight, 1 - mix_weight))   # IPW RMP

      # Using the RMP, sample from the posterior for the control survival probability at time t.
      # We first construct the posterior for the control log(shape) and intercept, and then
      # from this we derive the posterior distribution of the survival probability.
      post_control <- calc_post_weibull(filter(int_df_admin_cen, trt == 0),
                                        response = y,
                                        event = event,
                                        prior = mix_prior,
                                        analysis_time = surv_prob_time)
      samp_control <- unlist(rstan::extract(post_control, pars = c("survProb")))   # posterior sample
      mean_cont <- mean(samp_control)    # posterior mean of control survival prob at time t

      # Sample from the posterior for the treatment survival probability at time t. We extract
      # the vague portion of the RMP and use this as the prior for the log(shape) and intercept.
      vague_prior <- dist_multivariate_normal(mu = list(mix_means(mix_prior)[[2]]),
                                              sigma = list(mix_sigmas(mix_prior)[[2]]))
      post_treated <- calc_post_weibull(filter(int_df_admin_cen, trt == 1),
                                        response = y,
                                        event = event,
                                        prior = vague_prior,
                                        analysis_time = surv_prob_time)
      samp_trt <- unlist(rstan::extract(post_treated, pars = c("survProb")))   # posterior sample

      # Sample from the posterior distribution for the control survival probability at time t
      # without borrowing (needed for ESS calculation)
      post_control_no_borrow <- calc_post_weibull(filter(int_df_admin_cen, trt == 0),
                                                  response = y,
                                                  event = event,
                                                  prior = vague_prior,
                                                  analysis_time = surv_prob_time)
      samp_no_borrow <- unlist(rstan::extract(post_control_no_borrow, pars = c("survProb")))  # posterior sample

      # Obtain a posterior sample of the marginal treatment effect (difference in treatment
      # and control survival probabilities at time t)
      samp_trt_diff <- samp_trt - samp_control
      mean_trt_diff <- mean(samp_trt_diff)

      # Test H0: trt diff <= 0 vs. H1: trt diff > 0. Reject H0 if P(trt diff > 0|data) > 1 - alpha
      trt_diff_prob <- mean(samp_trt_diff > 0)   # posterior probability P(trt diff > 0|data)
      reject_H0_yes <- trt_diff_prob > .975      # H0 rejection indicator for alpha = 0.025

      # Hypothesis testing for no borrowing
      samp_trt_diff_no_borrow <- samp_trt - samp_no_borrow
      trt_diff_prob_no_borrow <- mean(samp_trt_diff_no_borrow > 0)   # posterior probability P(trt diff > 0|data)
      reject_H0_yes_no_borrow <- trt_diff_prob_no_borrow > .975      # H0 rejection indicator for alpha = 0.025

      # Calculate the effective sample size (ESS) of the posterior distribution of the control
      # survival probability at time t
      var_no_borrow <- variance(samp_no_borrow)         # post variance of control without borrowing
      var_borrow <- variance(samp_control)              # post variance of control with borrowing
      ess <- n_int_cont * var_no_borrow / var_borrow    # effective sample size

      # Add any iteration-specific summary statistic to this list of outputs
      list(
        "scenario" = scenario,                                     # scenario number
        "iter_id" = iter_id,                                       # iteration number
        "mean_post_cont" = mean_cont,                              # posterior mean of ctrl surv prob at t
        "median_post_cont" = median(samp_control),                 # posterior median of ctrl surv prob at t
        "variance_post_cont" = variance(samp_control),             # posterior variance of ctrl surv prob at t
        "q025_post_cont" = quantile(samp_control, .025),           # lower limit of 95% CrI of ctrl surv prob at t
        "p975_post_cont" = quantile(samp_control, .975),           # upper limit of 95% CrI of ctrl surv prob at t
        "mean_trt_diff" = mean_trt_diff,                           # posterior mean of the trt difference
        "trt_diff_prob" = trt_diff_prob,                           # post probability P(trt diff > 0|data)
        "reject_H0_yes" = reject_H0_yes,                           # H0 rejection indicator
        "no_borrowing_trt_diff_prob" = trt_diff_prob_no_borrow,    # post probability P(trt diff > 0|data) under bo borrowing
        "no_borrowing_reject_H0_yes" = reject_H0_yes_no_borrow,    # H0 rejection indicator under no borrowing
        "ess" = ess,                                               # posterior ESS of post dist for ctrl RR
        "irrt_bias_trteff" = mean_trt_diff - marg_trt_eff,         # contribution to bias of mean trt diff at t
        "irrt_mse_trteff" = (mean_trt_diff - marg_trt_eff)^2,      # contribution to MSE of mean trt diff at t
        "irrt_bias_cont" = mean_cont - true_control_surv_prob,     # contribution to bias of mean ctrl surv prob at t
        "irrt_mse_cont" = (mean_cont - true_control_surv_prob)^2,  # contribution to MSE of mean ctrl surv prob at t
        "pwr_prior" = pwr_prior,                                   # IPW power prior
        "mix_prior" = mix_prior                                    # mixture prior
      )

    }
    )
    output

  }, .options = para_opts, .progress = TRUE) |>
  bind_rows()

# Combine the output from the scenarios with the parameters of each simulation
combined_output <- all_sims |>
  left_join(sim_output, by = c("scenario", "iter_id"))

# Save results for all iterations
# save(combined_output, "Location where results should be saved.rda")

# Get the column names of everything that we want to summaries by (i.e.,
# everything but the iterations)
grouping_vars <- colnames(all_sims) |>
  discard(\(x) x == "iter_id")

# Final output
output <- combined_output |>
  summarise(across(mean_post_cont:irrt_mse_cont, mean),
            .by = all_of(grouping_vars))

# Save aggregate results for each scenario
# save(output, file = "Location where results should be saved.rda")


