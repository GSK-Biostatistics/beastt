# Simulation Overview ---------------------------------------------------------

# This script allows users to conduct a simulation study using the IPW BDB
# methodology for a binary outcome as outlined in Psioda et al. (2025):
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

# Set up parallelization where "workers" indicates the number of cores
plan(multisession, workers = availableCores()-1)
para_opts <- furrr_options(seed = TRUE)

# Simulation Setup -------------------------------------------------------------

# Step 1: Read in external control data
external_dat <- beastt::ex_binary_df
# external_dat <- read_csv("Location of the external data")

# Model "true" regression coefficients corresponding to intercept and covariates
# effects using a logistic regression model with the external data
logit_mod <- glm(y ~ cov1 + cov2 + cov3 + cov4, data = external_dat, family = binomial)
# logit_mod <- glm(#YOUR MODEL HERE,
#                  data = external_dat, family = binomial)

# Step 2: Define the underlying population characteristics to vary. Here, drift
# is the difference in the marginal control response rates (RR) between the
# external and internal trials attributed to unmeasured confounding, and
# treatment effect is the difference in marginal RRs between the treated and
# control populations.
drift_RR = seq(-0.16, 0.16, by = 0.02)   # Internal vs External (positive drift = internal RR is higher)
trt_effect_RR = c(0, .1, .15)            # Treatment vs Control (positive TE = treatment RR is higher)


# Step 3: Convert the drift and treatment effects from marginal to conditional models

# Later, we will generate response data for the internal arms using a conditional
# outcome model (i.e., logistic regression) that assumes a relationship between
# the covariates and the response. To account for the specified drift and treatment
# effect, we first need to convert these effects from the marginal scale to the
# conditional scale. For more information about this process, please refer to the
# function details for `calc_cond_binary`.

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
imbal_pop_1 <- bootstrap_cov(ex_dat_cov, pop_size, imbal_var = cov2, imbal_prop = c(0.25, 0.5, 0.75))
imbal_pop_2 <- bootstrap_cov(ex_dat_cov, pop_size, imbal_var = cov4, imbal_prop = 0.5)

# Combine all populations into a list
pop_ls <- c(list(no_imbal_pop, imbal_pop_2), imbal_pop_1)
# Naming the populations will make them easier to identify later
names(pop_ls) <- c("no imbalance", "cov4: 0.5", "cov2: 0.25", "cov2: 0.5", "cov2: 0.75")

# 3 b) For each population, identify the conditional drift and treatment effect values
# that best correspond to the specified values of marginal drift and treatment effect
pop_var <- pop_ls |>
  future_map(calc_cond_binary, logit_mod, drift_RR, trt_effect_RR) |> # returns a list of data frames
  bind_rows(.id = "population")  # combines list into 1 data frame with a population column


# Step 4: Create a data frame of all possible simulations using your population
# variables (e.g., "drift_RR", "trt_effect_RR", "imbal_var", "imbal_prop") as a
# foundation. Add any additional design variables you would like to vary. If you
# would like to compare priors, you can make a vector of distributional objects.
# Simulation scenarios are defined by unique combinations of the population
# variables and design variables.
all_sims <- pop_var |>
  crossing(
    # Internal control sample size
    n_int_cont = 65,

    # Internal treatment sample size
    n_int_trt = 130,

    # Initial beta prior incorporated into the power prior
    initial_prior = dist_beta(0.5, 0.5),

    # Vague prior incorporated into the RMP for the control RR and the posterior
    # distribution for the treatment RR
    vague_prior = dist_beta(0.5, 0.5),

    # Prior mixture weight associated with the informative component (i.e.,
    # IPW power prior) in the robust mixture prior
    mix_weight = 0.5
  ) |>
  mutate(scenario = row_number()) |>  # Add a scenario ID
  crossing(iter_id = c(1:1000))       # Add an iteration ID (within scenario)


# Simulations ------------------------------------------------------------------

# We now iterate over all rows in the simulation data frame and calculate
# operating characteristics for each scenario. The pmap and list functions make
# it possible to refer to each column of the data frame by its name. To step
# through this code, add browser()
sim_output <- all_sims |>
  future_pmap(function(...){
    output <- with(list(...), {

      # Simulate data ----------------------------------------------------------
      # Sample covariates from the scenario-specific population for the internal arms
      int_cont_df <- slice_sample(pop_ls[[population]], n = n_int_cont)  # control
      int_trt_df <- slice_sample(pop_ls[[population]], n = n_int_trt)    # treatment

      # Using the logistic model previously fit to the external control data, predict
      # the probability of response on the logit scale and adjust for the conditional
      # drift. Then, take the inverse logit to get the probability of response for
      # each simulated individual. Finally, sample a binary response.
      int_cont_df <- int_cont_df |>
        mutate(logit_pred = predict(logit_mod, int_cont_df) + conditional_drift,
               prob_response = inv_logit(logit_pred),
               y = rbinom(n_int_cont, 1, prob = prob_response),
               subjid = row_number())

      # For the treatment arm, do the same as above, but add both conditional
      # drift and conditional treatment effects
      int_trt_df <- int_trt_df |>
        mutate(logit_pred =
                 predict(logit_mod, int_trt_df) + conditional_drift + conditional_trt_eff,
               prob_response = inv_logit(logit_pred),
               y = rbinom(n_int_trt, 1, prob = prob_response),
               subjid = row_number())

      # Calculate the propensity scores and inverse probability weights for all control
      # participants (external and internal)
      ps_obj <- calc_prop_scr(internal_df = int_cont_df,
                              external_df = external_dat,
                              id_col = subjid,
                              # model = #YOUR MODEL HERE
                              model = ~ cov1 + cov2 + cov3 + cov4
      )

      # Calculate the inverse probability weighted power prior for the control RR using
      # the specified initial prior
      pwr_prior <- calc_power_prior_beta(ps_obj,
                                         response = y,
                                         prior = initial_prior)

      # "Robustify" the power prior by mixing it with the specified vague prior to create
      # an inverse probability weighted robust mixture prior (IPW RMP). Weight the power
      # prior (i.e., the informative component) using the specified mixture weight
      mix_prior <- dist_mixture(informative = pwr_prior,
                                vague = vague_prior,
                                weights = c(mix_weight, 1-mix_weight))

      # Calculate the posterior distribution for the control RR using the IPW RMP
      post_control <- calc_post_beta(int_cont_df,
                                     response = y,
                                     prior = mix_prior)
      mean_cont <- mean(post_control)     # posterior mean of the control RR

      # Calculate the posterior distribution for the control RR without borrowing
      post_control_no_borrow <- calc_post_beta(int_cont_df,
                                               response = y,
                                               prior = vague_prior)

      # Calculate the posterior distribution for the treatment RR
      post_trt <- calc_post_beta(int_trt_df,
                                 response = y,
                                 prior = vague_prior)

      # Obtain a posterior sample of the marginal treatment effect (risk difference)
      samp_control <- generate(x = post_control, times = 100000)[[1]]
      samp_trt <- generate(x = post_trt, times = 100000)[[1]]
      samp_trt_diff <- samp_trt - samp_control
      mean_trt_diff <- mean(samp_trt_diff)

      # Test H0: trt diff <= 0 vs. H1: trt diff > 0. Reject H0 if P(trt diff > 0|data) > 1 - alpha
      trt_diff_prob <- mean(samp_trt_diff > 0)   # posterior probability P(trt diff > 0|data)
      reject_H0_yes <- trt_diff_prob > .975      # H0 rejection indicator for alpha = 0.025

      # Calculate the effective sample size of the posterior distribution of the control RR
      var_no_borrow <- variance(post_control_no_borrow)   # post variance of control RR without borrowing
      var_borrow <- variance(post_control)                # post variance of control RR with borrowing
      ess <- n_int_cont * var_no_borrow / var_borrow      # effective sample size

      # Hypothesis testing for no borrowing
      samp_control_no_borrow <- generate(x = post_control_no_borrow, times = 100000)[[1]]
      samp_no_borrow_trt_diff <- samp_trt - samp_control_no_borrow
      no_borrowing_trt_diff_prob <- mean(samp_no_borrow_trt_diff > 0)   # posterior probability P(trt diff > 0|data)
      no_borrowing_reject_H0_yes <- no_borrowing_trt_diff_prob > .975   # H0 rejection indicator for alpha = 0.025


      # Add any iteration-specific summary statistic to this list of outputs
      list(
        "scenario" = scenario,                                     # scenario number
        "iter_id" = iter_id,                                       # iteration number
        "mean_post_cont" = mean_cont,                              # posterior mean of ctrl
        "median_post_cont" = median(samp_control),                 # posterior median of ctrl
        "variance_post_cont" = variance(samp_control),             # posterior variance of ctrl
        "q025_post_cont" = quantile(samp_control, .025),           # lower limit of 95% CrI of ctrl
        "p975_post_cont" = quantile(samp_control, .975),           # upper limit of 95% CrI of ctrl
        "mean_trt_diff" = mean_trt_diff,                           # posterior mean of the trt difference
        "post_mix_w" = parameters(post_control)$w[[1]]["informative"], #posterior control weight
        "trt_diff_prob" = trt_diff_prob,                           # post probability P(trt diff > 0|data)
        "reject_H0_yes" = reject_H0_yes,                           # H0 rejection indicator
        "no_borrowing_trt_diff_prob" =no_borrowing_trt_diff_prob,  # post probability P(trt diff > 0|data) under bo borrowing
        "no_borrowing_reject_H0_yes" = no_borrowing_reject_H0_yes, # H0 rejection indicator under bo borrowing
        "ess" = ess,                                               # posterior ESS of post dist for ctrl RR
        "irrt_bias_trteff" = mean_trt_diff - marg_trt_eff,         # contribution to bias of mean trt diff
        "irrt_mse_trteff" = (mean_trt_diff - marg_trt_eff)^2,      # contribution to MSE of mean trt diff
        "irrt_bias_cont" = mean_cont - true_control_RR,            # contribution to bias of mean ctrl
        "irrt_mse_cont" = (mean_cont - true_control_RR)^2,         # contribution to MSE of mean ctrl
        "pwr_prior" = pwr_prior                                    # IPW power prior
      )

    }
    )
    output

  }, .options = para_opts, .progress = TRUE) |>
  bind_rows()

# Combine the output from the scenarios with the parameters of each simulation
combined_output <- all_sims |>
  left_join(sim_output, by = c("scenario", "iter_id"))


# Get the column names of everything that we want to summaries by (i.e.

# everything but the iterations)
grouping_vars <- colnames(all_sims) |>
  discard(\(x) x == "iter_id")

# Final output
output <- combined_output |>
  summarise(across(mean_post_cont:irrt_mse_cont, mean),
            .by = all_of(grouping_vars))



