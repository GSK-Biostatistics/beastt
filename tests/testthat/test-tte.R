################################################################################
# Code to test various functions of beastt - TTE endpoint
################################################################################
skip_on_cran()
### Source R script with Stan code and compile Stan models
source("code_for_tte_tests_Stan.R")
# compile Stan model - power prior
stan_mod_pp <- rstan::stan_model(model_code = BDB_stan_pp)
# compile Stan model - posterior
stan_mod_post <- rstan::stan_model(model_code = BDB_stan_post)

################################################################################
# calc_power_prior_weibull
################################################################################

##### Test for calc_power_prior_weibull() with:
#####    (1) prop_scr_obj (not unweighted external data)
#####    (2) Laplace approximation
test_that("calc_power_prior_weibull returns correct values for first case", {
  ## Define values to be used for both beastt code and comparison code
  init_norm_prior <- dist_normal(0, 10)   # initial prior for intercept parameter
  init_hn_scale <- 50            # scale of the half-normal initial prior for the shape parameter
  init_prior <- dist_normal(mu = 0.5, sigma = 10)   # proper initial prior for power prior
  ps_obj <- calc_prop_scr(internal_df = filter(int_tte_df, trt == 0),
                          external_df = ex_tte_df,
                          id_col = subjid,
                          model = ~ cov1 + cov2 + cov3 + cov4)

  ## beastt code
  pwr_prior <- calc_power_prior_weibull(external_data = ps_obj,
                                        response = y,
                                        event = event,
                                        intercept = init_norm_prior,
                                        shape = init_hn_scale,
                                        approximation = "Laplace")
  pp_mean_beastt <- parameters(pwr_prior)$mu[[1]]
  pp_cov_beastt <- parameters(pwr_prior)$sigma[[1]]

  ## Comparison code
  ipws <- ps_obj$external_df$`___weight___`
  calc_neg_log_dens <- function(x, y_vec, event_vec, ipw_vec, beta0_mean, beta0_sd, alpha_scale){

    # Extract elements of x and save as parameters log(alpha) and beta0
    log_alpha <- x[1]
    beta0 <- x[2]

    # Log density
    log_dens <-
      sum( ipw_vec * event_vec * dweibull(x = y_vec, shape = exp(log_alpha), scale = exp(-beta0), log = TRUE) -
             ipw_vec * (1 - event_vec) * (y_vec / exp(-beta0))^exp(log_alpha) ) +   # log IWP likelihood of external data
      dnorm(x = beta0, mean = beta0_mean, sd = beta0_sd, log = TRUE) +       # log density of normal prior on beta0
      .5 * log(2/pi) - log(alpha_scale) - exp(2 * log_alpha) / (2 * alpha_scale^2)  # log dens of half norm prior on alpha
    -log_dens

  }
  inits <- c(0, 0)   # initial parameters for log(alpha) and beta0, respectively
  optim_pp_ctrl <- optim(par = inits,
                         fn = calc_neg_log_dens,
                         method = "Nelder-Mead",
                         y_vec = ex_tte_df$y,
                         event_vec = ex_tte_df$event,
                         ipw_vec = ipws,
                         beta0_mean = parameters(init_norm_prior)$mu,
                         beta0_sd = parameters(init_norm_prior)$sigma,
                         alpha_scale = init_hn_scale,
                         hessian = TRUE)
  pp_mean_comp <- optim_pp_ctrl$par              # mean vector
  pp_cov_comp <- solve(optim_pp_ctrl$hessian)    # covariance matrix

  ## Check that means and SDs of the normal power priors are equal using both methods
  expect_equal(pp_mean_beastt, pp_mean_comp)
  expect_equal(pp_cov_beastt, pp_cov_comp)

  ## Check that a distribution is returned
  expect_s3_class(pwr_prior, "distribution")
})

##### Test for calc_power_prior_weibull() with:
#####    (1) prop_scr_obj (not unweighted external data)
#####    (2) MCMC approximation
test_that("calc_power_prior_weibull returns correct values for second case", {
  ## Define values to be used for both beastt code and comparison code
  init_norm_prior <- dist_normal(0, 10)   # initial prior for intercept parameter
  init_hn_scale <- 50            # scale of the half-normal initial prior for the shape parameter
  init_prior <- dist_normal(mu = 0.5, sigma = 10)   # proper initial prior for power prior
  ps_obj <- calc_prop_scr(internal_df = filter(int_tte_df, trt == 0),
                          external_df = ex_tte_df,
                          id_col = subjid,
                          model = ~ cov1 + cov2 + cov3 + cov4)

  ## beastt code
  set.seed(123)
  pwr_prior <- calc_power_prior_weibull(external_data = ps_obj,
                                        response = y,
                                        event = event,
                                        intercept = init_norm_prior,
                                        shape = init_hn_scale,
                                        approximation = "MCMC")
  pp_mean_beastt <- parameters(pwr_prior)$mu[[1]]
  pp_cov_beastt <- parameters(pwr_prior)$sigma[[1]]

  ## Comparison code
  ipws <- ps_obj$external_df$`___weight___`
  data_input <- list(
    N = nrow(ex_tte_df),     # sample size of the external control arm
    y = ex_tte_df$y,         # N x 1 vector of observed times
    e = ex_tte_df$event,     # N x 1 vector of event indicators (1: event, 0: no event)
    wgt = ipws,              # N x 1 vector of IPWs
    beta_mean = parameters(init_norm_prior)$mu,    # mean for normal prior on the regression parameters
    beta_sd = parameters(init_norm_prior)$sigma,   # sd for normal prior on the regression parameters
    shape_mean = 0,                                # mean of half normal prior on weibull shape parameter
    shape_sd = init_hn_scale                       # sd of half normal prior on weibull shape parameter
  )
  set.seed(123)
  stan_fit <- sampling(stan_mod_pp, data = data_input, pars = c("beta", "shape"),
                       iter = 26000, warmup = 1000, chains = 4)
  pp_draws <- as.matrix(stan_fit)           # posterior samples of each parameter
  pp_mean_comp <- c(mean(log(pp_draws[,"shape"])), mean(pp_draws[,"beta"]))
  pp_cov_comp <- cov(cbind(log(pp_draws[,"shape"]), pp_draws[,"beta"]))
  names(pp_mean_comp) <- c("log_alpha", "beta0")
  dimnames(pp_cov_comp)[[1]] <- dimnames(pp_cov_comp)[[2]] <- c("log_alpha", "beta0")

  ## Check that means and SDs of the normal power priors are equal using both methods
  expect_equal(all(abs(pp_mean_beastt-pp_mean_comp) < 0.0025), TRUE)
  expect_equal(all(abs(pp_cov_beastt-pp_cov_comp) < 0.001), TRUE)

  ## Check that a distribution is returned
  expect_s3_class(pwr_prior, "distribution")
})

##### Test for calc_power_prior_weibull() with:
#####    (1) unweighted external data (read in external df, not prop_scr_obj)
#####    (2) Laplace approximation
test_that("calc_power_prior_weibull returns correct values for third case", {
  ## Define values to be used for both beastt code and comparison code
  init_norm_prior <- dist_normal(0, 10)   # initial prior for intercept parameter
  init_hn_scale <- 50            # scale of the half-normal initial prior for the shape parameter
  init_prior <- dist_normal(mu = 0.5, sigma = 10)   # proper initial prior for power prior

  ## beastt code
  pwr_prior <- calc_power_prior_weibull(external_data = ex_tte_df,
                                        response = y,
                                        event = event,
                                        intercept = init_norm_prior,
                                        shape = init_hn_scale,
                                        approximation = "Laplace")
  pp_mean_beastt <- parameters(pwr_prior)$mu[[1]]
  pp_cov_beastt <- parameters(pwr_prior)$sigma[[1]]

  ## Comparison code
  calc_neg_log_dens <- function(x, y_vec, event_vec, ipw_vec, beta0_mean, beta0_sd, alpha_scale){

    # Extract elements of x and save as parameters log(alpha) and beta0
    log_alpha <- x[1]
    beta0 <- x[2]

    # Log density
    log_dens <-
      sum( ipw_vec * event_vec * dweibull(x = y_vec, shape = exp(log_alpha), scale = exp(-beta0), log = TRUE) -
             ipw_vec * (1 - event_vec) * (y_vec / exp(-beta0))^exp(log_alpha) ) +   # log IWP likelihood of external data
      dnorm(x = beta0, mean = beta0_mean, sd = beta0_sd, log = TRUE) +       # log density of normal prior on beta0
      .5 * log(2/pi) - log(alpha_scale) - exp(2 * log_alpha) / (2 * alpha_scale^2)  # log dens of half norm prior on alpha
    -log_dens

  }
  inits <- c(0, 0)   # initial parameters for log(alpha) and beta0, respectively
  optim_pp_ctrl <- optim(par = inits,
                         fn = calc_neg_log_dens,
                         method = "Nelder-Mead",
                         y_vec = ex_tte_df$y,
                         event_vec = ex_tte_df$event,
                         ipw_vec = rep(1, nrow(ex_tte_df)),
                         beta0_mean = parameters(init_norm_prior)$mu,
                         beta0_sd = parameters(init_norm_prior)$sigma,
                         alpha_scale = init_hn_scale,
                         hessian = TRUE)
  pp_mean_comp <- optim_pp_ctrl$par              # mean vector
  pp_cov_comp <- solve(optim_pp_ctrl$hessian)    # covariance matrix

  ## Check that means and SDs of the normal power priors are equal using both methods
  expect_equal(pp_mean_beastt, pp_mean_comp)
  expect_equal(pp_cov_beastt, pp_cov_comp)

  ## Check that a distribution is returned
  expect_s3_class(pwr_prior, "distribution")
})

##### Test for calc_power_prior_weibull() with:
#####    (1) unweighted external data (read in external df, not prop_scr_obj)
#####    (2) MCMC approximation
test_that("calc_power_prior_weibull returns correct values for fourth case", {
  ## Define values to be used for both beastt code and comparison code
  init_norm_prior <- dist_normal(0, 10)   # initial prior for intercept parameter
  init_hn_scale <- 50            # scale of the half-normal initial prior for the shape parameter
  init_prior <- dist_normal(mu = 0.5, sigma = 10)   # proper initial prior for power prior

  ## beastt code
  set.seed(123)
  pwr_prior <- calc_power_prior_weibull(external_data = ex_tte_df,
                                        response = y,
                                        event = event,
                                        intercept = init_norm_prior,
                                        shape = init_hn_scale,
                                        approximation = "MCMC")
  pp_mean_beastt <- parameters(pwr_prior)$mu[[1]]
  pp_cov_beastt <- parameters(pwr_prior)$sigma[[1]]

  ## Comparison code
  data_input <- list(
    N = nrow(ex_tte_df),    # sample size of the external control arm
    y = ex_tte_df$y,         # N x 1 vector of observed times
    e = ex_tte_df$event,     # N x 1 vector of event indicators (1: event, 0: no event)
    wgt = rep(1, nrow(ex_tte_df)),                 # N x 1 vector of 1s
    beta_mean = parameters(init_norm_prior)$mu,    # mean for normal prior on the regression parameters
    beta_sd = parameters(init_norm_prior)$sigma,   # sd for normal prior on the regression parameters
    shape_mean = 0,                                # mean of half normal prior on weibull shape parameter
    shape_sd = init_hn_scale                       # sd of half normal prior on weibull shape parameter
  )
  set.seed(123)
  stan_fit <- sampling(stan_mod_pp, data = data_input, pars = c("beta", "shape"),
                       iter = 26000, warmup = 1000, chains = 4)
  pp_draws <- as.matrix(stan_fit)           # posterior samples of each parameter
  pp_mean_comp <- c(mean(log(pp_draws[,"shape"])), mean(pp_draws[,"beta"]))
  pp_cov_comp <- cov(cbind(log(pp_draws[,"shape"]), pp_draws[,"beta"]))
  names(pp_mean_comp) <- c("log_alpha", "beta0")
  dimnames(pp_cov_comp)[[1]] <- dimnames(pp_cov_comp)[[2]] <- c("log_alpha", "beta0")

  ## Check that means and SDs of the normal power priors are equal using both methods
  expect_equal(all(abs(pp_mean_beastt-pp_mean_comp) < 0.002), TRUE)
  expect_equal(all(abs(pp_cov_beastt-pp_cov_comp) < 0.001), TRUE)

  ## Check that a distribution is returned
  expect_s3_class(pwr_prior, "distribution")
})

# Test for invalid external data
test_that("calc_power_prior_weibull handles invalid external data", {
  expect_error(calc_power_prior_weibull(external_data = c(1, 2, 3),
                                        response = y,
                                        event = event,
                                        intercept = dist_normal(0, 10),
                                        shape = 50,
                                        approximation = "Laplace"))
  expect_error(calc_power_prior_weibull(external_data = "abc",
                                        response = y,
                                        event = event,
                                        intercept = dist_normal(0, 10),
                                        shape = 50,
                                        approximation = "Laplace"))
})

# Test for invalid intercept
test_that("calc_power_prior_weibull handles invalid intercept", {
  expect_error(calc_power_prior_weibull(external_data = ex_tte_df,
                                        response = y,
                                        event = event,
                                        intercept = "a",
                                        shape = 50,
                                        approximation = "Laplace"))
  expect_error(calc_power_prior_weibull(external_data = ex_tte_df,
                                        response = y,
                                        event = event,
                                        intercept = c(1, 2, 3),
                                        shape = 50,
                                        approximation = "Laplace"))
})

# Test for invalid shape
test_that("calc_power_prior_weibull handles invalid shape", {
  expect_error(calc_power_prior_weibull(external_data = ex_tte_df,
                                        response = y,
                                        event = event,
                                        intercept = dist_normal(0, 10),
                                        shape = "a",
                                        approximation = "Laplace"))
  expect_error(calc_power_prior_weibull(external_data = ex_tte_df,
                                        response = y,
                                        event = event,
                                        intercept = dist_normal(0, 10),
                                        shape = c(1, 2, 3),
                                        approximation = "Laplace"))
})

# Test for invalid approximation method
test_that("calc_power_prior_weibull handles invalid approximation method", {
  expect_error(calc_power_prior_weibull(external_data = ex_tte_df,
                                        response = y,
                                        event = event,
                                        intercept = dist_normal(0, 10),
                                        shape = 50,
                                        approximation = 12))
  expect_error(calc_power_prior_weibull(external_data = ex_tte_df,
                                        response = y,
                                        event = event,
                                        intercept = dist_normal(0, 10),
                                        shape = 50,
                                        approximation = c(5, 10)))
  expect_error(calc_power_prior_weibull(external_data = ex_tte_df,
                                        response = y,
                                        event = event,
                                        intercept = dist_normal(0, 10),
                                        shape = 50,
                                        approximation = "zero"))
})

################################################################################
# calc_post_weibull
################################################################################

##### Test for calc_post_weibull() with:
#####    (1) MVN prior (not a mixture of two MVN priors)
#####    (2) one analysis time
test_that("calc_post_weibull returns correct values for first case", {
  ## Define values to be used for both beastt code and comparison code
  MVN_prior <- dist_multivariate_normal(mu = list(c(0.5, 0.6)),
                                        sigma = list(matrix(c(10, -.5, -.5, 10), nrow = 2)))   # normal prior
  analysis_time <- 12

  ## beastt code
  set.seed(123)
  post_dist <- calc_post_weibull(internal_data = filter(int_tte_df, trt == 0),
                                 response = y,
                                 event = event,
                                 prior = MVN_prior,
                                 analysis_time = analysis_time)
  surv_prob_beastt <- as.data.frame(extract(post_dist, pars = c("survProb")))[,1]
  post_mean_beastt <- mean(surv_prob_beastt)
  post_var_beastt <- var(surv_prob_beastt)

  ## Comparison code
  analysis_time_comp <- c(analysis_time, analysis_time)  # making dimension >1 so code will run
  data_input <- list(
    N = nrow(filter(int_tte_df, trt == 0)),     # sample size of the internal control arm
    y = filter(int_tte_df, trt == 0)$y,         # N x 1 vector of observed times
    e = filter(int_tte_df, trt == 0)$event,     # N x 1 vector of event indicators (1: event, 0: no event)
    T = length(analysis_time_comp),             # length of analysis times
    times = c(analysis_time_comp),              # analysis time
    cov_inf = parameters(MVN_prior)$sigma[[1]],      # informative prior component covariance matrix
    cov_vague = parameters(MVN_prior)$sigma[[1]],    # vague prior component covariance matrix
    mu_inf = parameters(MVN_prior)$mu[[1]],          # informative prior component mean vector
    mu_vague = parameters(MVN_prior)$mu[[1]],        # vague prior component mean vector
    w0 = 1                                           # prior mixture weight
  )
  set.seed(123)
  stan_fit <- sampling(stan_mod_post, data = data_input, pars = "survProb",
                       iter = 26000, warmup = 1000, chains = 4)
  post_draws <- as.matrix(stan_fit)           # posterior samples of each parameter
  post_mean_comp <- mean(post_draws[,"survProb[1]"])
  post_var_comp <- var(post_draws[,"survProb[1]"])

  ## Check that means and SDs of the normal power priors are equal using both methods
  expect_equal(abs(post_mean_beastt-post_mean_comp) < 0.001, TRUE)
  expect_equal(abs(post_var_beastt-post_var_comp) < 0.0001, TRUE)

  ## Check that a stanfit is returned
  expect_s4_class(post_dist, "stanfit")
})

##### Test for calc_post_weibull() with:
#####    (1) MVN prior (not a mixture of two MVN priors)
#####    (2) two analysis times
test_that("calc_post_weibull returns correct values for second case", {
  ## Define values to be used for both beastt code and comparison code
  MVN_prior <- dist_multivariate_normal(mu = list(c(0.5, 0.6)),
                                        sigma = list(matrix(c(10, -.5, -.5, 10), nrow = 2)))   # normal prior
  analysis_times <- c(6, 8)

  ## beastt code
  set.seed(123)
  post_dist <- calc_post_weibull(internal_data = filter(int_tte_df, trt == 0),
                                 response = y,
                                 event = event,
                                 prior = MVN_prior,
                                 analysis_time = analysis_times)
  surv_prob_beastt <- as.data.frame(extract(post_dist, pars = c("survProb")))
  post_means_beastt <- colMeans(surv_prob_beastt)
  post_vars_beastt <- apply(surv_prob_beastt, 2, var)

  ## Comparison code
  data_input <- list(
    N = nrow(filter(int_tte_df, trt == 0)),     # sample size of the internal control arm
    y = filter(int_tte_df, trt == 0)$y,         # N x 1 vector of observed times
    e = filter(int_tte_df, trt == 0)$event,     # N x 1 vector of event indicators (1: event, 0: no event)
    T = length(analysis_times),                 # length of analysis times
    times = analysis_times,                     # analysis times
    cov_inf = parameters(MVN_prior)$sigma[[1]],      # informative prior component covariance matrix
    cov_vague = parameters(MVN_prior)$sigma[[1]],    # vague prior component covariance matrix
    mu_inf = parameters(MVN_prior)$mu[[1]],          # informative prior component mean vector
    mu_vague = parameters(MVN_prior)$mu[[1]],        # vague prior component mean vector
    w0 = 1                                           # prior mixture weight
  )
  set.seed(123)
  stan_fit <- sampling(stan_mod_post, data = data_input, pars = "survProb",
                       iter = 26000, warmup = 1000, chains = 4)
  post_draws <- as.matrix(stan_fit)           # posterior samples of each parameter
  post_means_comp <- colMeans(post_draws[,c("survProb[1]", "survProb[2]")])
  post_vars_comp <- apply(post_draws[,c("survProb[1]", "survProb[2]")], 2, var)

  ## Check that means and SDs of the normal power priors are equal using both methods
  expect_equal(all(abs(post_means_beastt-post_means_comp) < 0.001), TRUE)
  expect_equal(all(abs(post_vars_beastt-post_vars_comp) < 0.0001), TRUE)

  ## Check that a stanfit is returned
  expect_s4_class(post_dist, "stanfit")
})

##### Test for calc_post_weibull() with:
#####    (1) mixture prior (two MVN prior components)
#####    (2) one analysis time
test_that("calc_post_weibull returns correct values for third case", {
  ## Define values to be used for both beastt code and comparison code
  MVN_prior <- dist_multivariate_normal(mu = list(c(0.5, 0.6)),
                                        sigma = list(matrix(c(10, -.5, -.5, 10), nrow = 2)))   # normal prior
  mix_prior <- robustify_mvnorm(MVN_prior, n = sum(ex_tte_df$event), weights = c(.5, .5))
  analysis_time <- 12

  ## beastt code
  set.seed(123)
  post_dist <- calc_post_weibull(internal_data = filter(int_tte_df, trt == 0),
                                 response = y,
                                 event = event,
                                 prior = mix_prior,
                                 analysis_time = analysis_time)
  surv_prob_beastt <- as.data.frame(extract(post_dist, pars = c("survProb")))[,1]
  post_mean_beastt <- mean(surv_prob_beastt)
  post_var_beastt <- var(surv_prob_beastt)

  ## Comparison code
  analysis_time_comp <- c(analysis_time, analysis_time)  # making dimension >1 so code will run
  data_input <- list(
    N = nrow(filter(int_tte_df, trt == 0)),     # sample size of the internal control arm
    y = filter(int_tte_df, trt == 0)$y,         # N x 1 vector of observed times
    e = filter(int_tte_df, trt == 0)$event,     # N x 1 vector of event indicators (1: event, 0: no event)
    T = length(analysis_time_comp),             # length of analysis times
    times = c(analysis_time_comp),              # analysis time
    cov_inf = mix_sigmas(mix_prior)$informative,    # informative prior component covariance matrix
    cov_vague = mix_sigmas(mix_prior)$vague,        # vague prior component covariance matrix
    mu_inf = mix_means(mix_prior)$informative,      # informative prior component mean vector
    mu_vague = mix_means(mix_prior)$vague,          # vague prior component mean vector
    w0 = .5                                         # prior mixture weight
  )
  set.seed(123)
  stan_fit <- sampling(stan_mod_post, data = data_input, pars = "survProb",
                       iter = 26000, warmup = 1000, chains = 4)
  post_draws <- as.matrix(stan_fit)           # posterior samples of each parameter
  post_mean_comp <- mean(post_draws[,"survProb[1]"])
  post_var_comp <- var(post_draws[,"survProb[1]"])

  ## Check that means and SDs of the normal power priors are equal using both methods
  expect_equal(abs(post_mean_beastt-post_mean_comp) < 0.001, TRUE)
  expect_equal(abs(post_var_beastt-post_var_comp) < 0.0001, TRUE)

  ## Check that a stanfit is returned
  expect_s4_class(post_dist, "stanfit")
})

##### Test for calc_post_weibull() with:
#####    (1) mixture prior (two MVN prior components)
#####    (2) two analysis times
test_that("calc_post_weibull returns correct values for fourth case", {
  ## Define values to be used for both beastt code and comparison code
  MVN_prior <- dist_multivariate_normal(mu = list(c(0.5, 0.6)),
                                        sigma = list(matrix(c(10, -.5, -.5, 10), nrow = 2)))   # normal prior
  mix_prior <- robustify_mvnorm(MVN_prior, n = sum(ex_tte_df$event), weights = c(.5, .5))
  analysis_times <- c(6, 8)

  ## beastt code
  set.seed(123)
  post_dist <- calc_post_weibull(internal_data = filter(int_tte_df, trt == 0),
                                 response = y,
                                 event = event,
                                 prior = MVN_prior,
                                 analysis_time = analysis_times)
  surv_prob_beastt <- as.data.frame(extract(post_dist, pars = c("survProb")))
  post_means_beastt <- colMeans(surv_prob_beastt)
  post_vars_beastt <- apply(surv_prob_beastt, 2, var)

  ## Comparison code
  data_input <- list(
    N = nrow(filter(int_tte_df, trt == 0)),     # sample size of the internal control arm
    y = filter(int_tte_df, trt == 0)$y,         # N x 1 vector of observed times
    e = filter(int_tte_df, trt == 0)$event,     # N x 1 vector of event indicators (1: event, 0: no event)
    T = length(analysis_times),                 # length of analysis times
    times = analysis_times,                     # analysis times
    cov_inf = mix_sigmas(mix_prior)$informative,    # informative prior component covariance matrix
    cov_vague = mix_sigmas(mix_prior)$vague,        # vague prior component covariance matrix
    mu_inf = mix_means(mix_prior)$informative,      # informative prior component mean vector
    mu_vague = mix_means(mix_prior)$vague,          # vague prior component mean vector
    w0 = .5                                         # prior mixture weight
  )
  set.seed(123)
  stan_fit <- sampling(stan_mod_post, data = data_input, pars = "survProb",
                       iter = 26000, warmup = 1000, chains = 4)
  post_draws <- as.matrix(stan_fit)           # posterior samples of each parameter
  post_means_comp <- colMeans(post_draws[,c("survProb[1]", "survProb[2]")])
  post_vars_comp <- apply(post_draws[,c("survProb[1]", "survProb[2]")], 2, var)

  ## Check that means and SDs of the normal power priors are equal using both methods
  expect_equal(all(abs(post_means_beastt-post_means_comp) < 0.001), TRUE)
  expect_equal(all(abs(post_vars_beastt-post_vars_comp) < 0.0001), TRUE)

  ## Check that a stanfit is returned
  expect_s4_class(post_dist, "stanfit")
})

# Test for invalid internal data
test_that("calc_post_weibull handles invalid internal data", {
  expect_error(calc_post_weibull(internal_data = c(1, 2, 3),
                                 response = y,
                                 event = event,
                                 prior = dist_multivariate_normal(mu = list(c(0.5, 0.6)),
                                                                  sigma = list(matrix(c(10, -.5, -.5, 10), nrow = 2))),
                                 analysis_time = c(12, 24)))
  expect_error(calc_post_weibull(internal_data = "a",
                                 response = y,
                                 event = event,
                                 prior = dist_multivariate_normal(mu = list(c(0.5, 0.6)),
                                                                  sigma = list(matrix(c(10, -.5, -.5, 10), nrow = 2))),
                                 analysis_time = c(12, 24)))
})

# Test for invalid intercept
test_that("calc_post_weibull handles invalid prior", {
  expect_error(calc_post_weibull(internal_data = filter(int_tte_df, trt == 0),
                                 response = y,
                                 event = event,
                                 prior = c(1, 2, 3),
                                 analysis_time = c(12, 24)))
  expect_error(calc_post_weibull(internal_data = filter(int_tte_df, trt == 0),
                                 response = y,
                                 event = event,
                                 prior = "a",
                                 analysis_time = c(12, 24)))
})

# Test for invalid analysis time
test_that("calc_post_weibull handles invalid analysis_time", {
  expect_error(calc_post_weibull(internal_data = filter(int_tte_df, trt == 0),
                                 response = y,
                                 event = event,
                                 prior = dist_multivariate_normal(mu = list(c(0.5, 0.6)),
                                                                  sigma = list(matrix(c(10, -.5, -.5, 10), nrow = 2))),
                                 analysis_time = "a"))
  expect_error(calc_post_weibull(internal_data = filter(int_tte_df, trt == 0),
                                 response = y,
                                 event = event,
                                 prior = dist_multivariate_normal(mu = list(c(0.5, 0.6)),
                                                                  sigma = list(matrix(c(10, -.5, -.5, 10), nrow = 2))),
                                 analysis_time = dist_normal(0, 1)))
})
