################################################################################
# Code to test various functions of beastt - normal endpoint
################################################################################
skip_on_cran()
### Source R script with Stan code and compile Stan models
source("code_for_normal_tests_Stan.R")

# compile Stan model - sigma2 known
stan_mod_sigma2_known <- rstan::stan_model(model_code = BDB_stan_sigma2_known)

# compile Stan model - sigma2 unknown
stan_mod_sigma2_unknown <- rstan::stan_model(model_code = BDB_stan_sigma2_unknown)

################################################################################
# calc_power_prior_norm
################################################################################

##### Test for calc_power_prior_norm() with:
#####    (1) prop_scr_obj (not unweighted external data)
#####    (2) normal prior
#####    (3) external_sd known
test_that("calc_power_prior_norm returns correct values for first case", {
  ## Define values to be used for both beastt code and comparison code
  sd_external_control <- 0.15    # known SD for external control
  init_prior <- dist_normal(mu = 0.5, sigma = 10)   # proper initial prior for power prior
  ps_obj <- calc_prop_scr(internal_df = filter(int_norm_df, trt == 0),
                          external_df = ex_norm_df,
                          id_col = subjid,
                          model = ~ cov1 + cov2 + cov3 + cov4)

  ## beastt code
  pwr_prior <- calc_power_prior_norm(external_data = ps_obj,
                                     response = y,
                                     prior = init_prior,
                                     external_sd = sd_external_control)
  pp_mean_beastt <- parameters(pwr_prior)$mu
  pp_sd_beastt <- parameters(pwr_prior)$sigma

  ## Comparison code
  ipws <- ps_obj$external_df$`___weight___`
  pp_sd_comp <- sqrt( (1/parameters(init_prior)$sigma^2 + sum(ipws)/sd_external_control^2)^-1 )
  pp_mean_comp <- pp_sd_comp^2 * (parameters(init_prior)$mu/parameters(init_prior)$sigma^2 +
                                    sum(ipws * ex_norm_df$y)/sd_external_control^2)

  ## Check that means and SDs of the normal power priors are equal using both methods
  expect_equal(pp_mean_beastt, pp_mean_comp)
  expect_equal(pp_sd_beastt, pp_sd_comp)

  ## Check that a distribution is returned
  expect_s3_class(pwr_prior, "distribution")

})

##### Test for calc_power_prior_norm() with:
#####    (1) prop_scr_obj (not unweighted external data)
#####    (2) prior = NULL
#####    (3) external_sd known
test_that("calc_power_prior_norm returns correct values for second case", {
  ## Define values to be used for both beastt code and comparison code
  sd_external_control <- 0.15
  ps_obj <- calc_prop_scr(internal_df = filter(int_norm_df, trt == 0),
                          external_df = ex_norm_df,
                          id_col = subjid,
                          model = ~ cov1 + cov2 + cov3 + cov4)

  ## beastt code
  pwr_prior <- calc_power_prior_norm(external_data = ps_obj,
                                     response = y,
                                     prior = NULL,
                                     external_sd = sd_external_control)
  pp_mean_beastt <- parameters(pwr_prior)$mu
  pp_sd_beastt <- parameters(pwr_prior)$sigma

  ## Comparison code
  ipws <- ps_obj$external_df$`___weight___`
  pp_mean_comp <- sum(ipws * ex_norm_df$y) / sum(ipws)
  pp_sd_comp <- sqrt( sd_external_control^2 / sum(ipws) )

  ## Check that means and SDs of the normal power priors are equal using both methods
  expect_equal(pp_mean_beastt, pp_mean_comp)
  expect_equal(pp_sd_beastt, pp_sd_comp)

  ## Check that a distribution is returned
  expect_s3_class(pwr_prior, "distribution")

})

##### Test for calc_power_prior_norm() with:
#####    (1) prop_scr_obj (not unweighted external data)
#####    (2) prior = NULL
#####    (3) external_sd unknown
test_that("calc_power_prior_norm returns correct values for third case", {
  ## Define values to be used for both beastt code and comparison code
  ps_obj <- calc_prop_scr(internal_df = filter(int_norm_df, trt == 0),
                          external_df = ex_norm_df,
                          id_col = subjid,
                          model = ~ cov1 + cov2 + cov3 + cov4)

  ## beastt code
  pwr_prior <- calc_power_prior_norm(external_data = ps_obj,
                                     response = y,
                                     prior = NULL,
                                     external_sd = NULL)
  pp_mean_beastt <- parameters(pwr_prior)$mu
  pp_loc_beastt <- parameters(pwr_prior)$sigma
  pp_df_beastt <- parameters(pwr_prior)$df

  ## Comparison code
  ipws <- ps_obj$external_df$`___weight___`
  Z <- as.matrix(rep(1, nrow(ex_norm_df)))
  A <- diag(ipws)
  V <- Z %*% solve( t(Z) %*% A %*% Z) %*% t(Z) %*% A
  pp_mean_comp <- as.numeric(solve( t(Z) %*% A %*% Z) %*% t(Z) %*% A %*% as.matrix(ex_norm_df$y))
  pp_df_comp <- nrow(ex_norm_df) - 1
  pp_loc_comp <- as.numeric(sqrt( pp_df_comp^-1 * (t(as.matrix(ex_norm_df$y)) %*%
                                          (A - t(V) %*% A %*% V) %*% as.matrix(ex_norm_df$y)) %*%
                         solve( t(Z) %*% A %*% Z) ))

  ## Check that mean, location, and degrees of freedom of the t power priors are equal using both methods
  expect_equal(pp_mean_beastt, pp_mean_comp)
  expect_equal(pp_loc_beastt, pp_loc_comp)
  expect_equal(pp_df_beastt, pp_df_comp)

  ## Check that a distribution is returned
  expect_s3_class(pwr_prior, "distribution")
})

##### Test for calc_power_prior_norm() with:
#####    (1) unweighted external data (read in external df, not prop_scr_obj)
#####    (2) normal prior
#####    (3) external_sd known
test_that("calc_power_prior_norm returns correct values for fourth case", {
  ## Define values to be used for both beastt code and comparison code
  sd_external_control <- 0.15    # known SD for external control
  init_prior <- dist_normal(mu = 0.5, sigma = 10)   # proper initial prior for power prior

  ## beastt code
  pwr_prior <- calc_power_prior_norm(external_data = ex_norm_df,
                                     response = y,
                                     prior = init_prior,
                                     external_sd = sd_external_control)
  pp_mean_beastt <- parameters(pwr_prior)$mu
  pp_sd_beastt <- parameters(pwr_prior)$sigma

  ## Comparison code
  pp_sd_comp <- sqrt( (1/parameters(init_prior)$sigma^2 + nrow(ex_norm_df)/sd_external_control^2)^-1 )
  pp_mean_comp <- pp_sd_comp^2 * (parameters(init_prior)$mu/parameters(init_prior)$sigma^2 +
                                    sum(ex_norm_df$y)/sd_external_control^2)

  ## Check that means and SDs of the normal power priors are equal using both methods
  expect_equal(pp_mean_beastt, pp_mean_comp)
  expect_equal(pp_sd_beastt, pp_sd_comp)

  ## Check that a distribution is returned
  expect_s3_class(pwr_prior, "distribution")
})

##### Test for calc_power_prior_norm() with:
#####    (1) unweighted external data (read in external df, not prop_scr_obj)
#####    (2) prior = NULL
#####    (3) external_sd known
test_that("calc_power_prior_norm returns correct values for fifth case", {
  ## Define values to be used for both beastt code and comparison code
  sd_external_control <- 0.15

  ## beastt code
  pwr_prior <- calc_power_prior_norm(external_data = ex_norm_df,
                                     response = y,
                                     prior = NULL,
                                     external_sd = sd_external_control)
  pp_mean_beastt <- parameters(pwr_prior)$mu
  pp_sd_beastt <- parameters(pwr_prior)$sigma

  ## Comparison code
  pp_mean_comp <- mean(ex_norm_df$y)
  pp_sd_comp <- sqrt( sd_external_control^2 / nrow(ex_norm_df) )

  ## Check that means and SDs of the normal power priors are equal using both methods
  expect_equal(pp_mean_beastt, pp_mean_comp)
  expect_equal(pp_sd_beastt, pp_sd_comp)

  ## Check that a distribution is returned
  expect_s3_class(pwr_prior, "distribution")
})

##### Test for calc_power_prior_norm() with:
#####    (1) unweighted external data (read in external df, not prop_scr_obj)
#####    (2) prior = NULL
#####    (3) external_sd unknown
test_that("calc_power_prior_norm returns correct values for sixth case", {
  pwr_prior <- calc_power_prior_norm(external_data = ex_norm_df,
                                     response = y,
                                     prior = NULL,
                                     external_sd = NULL)
  pp_mean_beastt <- parameters(pwr_prior)$mu
  pp_loc_beastt <- parameters(pwr_prior)$sigma
  pp_df_beastt <- parameters(pwr_prior)$df

  ## Comparison code
  Z <- as.matrix(rep(1, nrow(ex_norm_df)))
  A <- diag(nrow(ex_norm_df))
  V <- Z %*% solve( t(Z) %*% A %*% Z) %*% t(Z) %*% A
  pp_mean_comp <- as.numeric(solve( t(Z) %*% A %*% Z) %*% t(Z) %*% A %*% as.matrix(ex_norm_df$y))
  pp_df_comp <- nrow(ex_norm_df) - 1
  pp_loc_comp <- as.numeric(pp_loc_comp <- sqrt( pp_df_comp^-1 * (t(as.matrix(ex_norm_df$y)) %*%
                                                         (A - t(V) %*% A %*% V) %*% as.matrix(ex_norm_df$y)) %*%
                                        solve( t(Z) %*% A %*% Z) ))

  ## Check that mean, location, and degrees of freedom of the t power priors are equal using both methods
  expect_equal(pp_mean_beastt, pp_mean_comp)
  expect_equal(pp_loc_beastt, pp_loc_comp)
  expect_equal(pp_df_beastt, pp_df_comp)

  ## Check that a distribution is returned
  expect_s3_class(pwr_prior, "distribution")
})

# Test for invalid external data
test_that("calc_power_prior_norm handles invalid external data", {
  expect_error(calc_power_prior_norm(c(1, 2, 3),
                                     response=y,
                                     prior=dist_normal(50, 10),
                                     external_sd = 0.15))
  expect_error(calc_power_prior_norm("abc",
                                     response=y,
                                     prior=dist_normal(50, 10),
                                     external_sd = 0.15))
})

# Test for invalid prior
test_that("calc_power_prior_norm handles invalid prior", {
  expect_error(calc_power_prior_norm(ex_norm_df,
                                     response=y,
                                     prior=5,
                                     external_sd = 0.15))
  expect_error(calc_power_prior_norm(ex_norm_df,
                                     response=y,
                                     prior="a",
                                     external_sd = 0.15))
  expect_error(calc_power_prior_norm(ex_norm_df,
                                     response=y,
                                     prior=dist_beta(5, 6),
                                     external_sd = 0.15))
})

# Test for invalid response variable
test_that("calc_power_prior_norm handles invalid response", {
  expect_error(calc_power_prior_norm(ex_norm_df,
                                     response=invalid,
                                     prior=dist_normal(50, 10),
                                     external_sd = 0.15))
})

# Test for invalid external sd
test_that("calc_power_prior_norm handles invalid external sd", {
  expect_error(calc_power_prior_norm(ex_norm_df,
                                     response=y,
                                     prior=dist_normal(50, 10),
                                     external_sd = -5))
  expect_error(calc_power_prior_norm(ex_norm_df,
                                     response=y,
                                     prior=dist_normal(50, 10),
                                     external_sd = "a"))
})

# Test for internal and external data with different response variable names
test_that("calc_power_prior_norm handles different response variable names", {
  int <- int_norm_df
  ex <- ex_norm_df |>
    dplyr::rename(y2 = y)
  init_prior <- dist_normal(mu = 0.5, sigma = 10)

  ps_obj <- calc_prop_scr(internal_df = filter(int, trt == 0),
                          external_df = ex,
                          id_col = subjid,
                          model = ~ cov1 + cov2 + cov3 + cov4)

  expect_error(calc_power_prior_norm(external_data = ps_obj,
                                     response = y,
                                     prior = init_prior))
})

################################################################################
# calc_post_norm
################################################################################

##### Test for calc_post_norm() with:
#####    (1) normal prior
#####    (2) internal_sd known
test_that("calc_post_norm returns the correct values for first case", {
  ## Define values to be used for both beastt code and comparison code
  sd_internal_control <- 0.15    # known SD for internal control
  norm_prior <- dist_normal(mu = 0.5, sigma = 10)   # normal prior

  ## beastt code
  post_dist <- calc_post_norm(filter(int_norm_df, trt == 0),
                              response = y,
                              prior = norm_prior,
                              internal_sd = sd_internal_control)
  post_mean_beastt <- parameters(post_dist)$mu
  post_sd_beastt <- parameters(post_dist)$sigma

  ## Comparison code
  post_sd_comp <- sqrt( (1/parameters(norm_prior)$sigma^2 +
                           nrow(filter(int_norm_df, trt == 0))/sd_internal_control^2)^-1 )
  post_mean_comp <- post_sd_comp^2 * (parameters(norm_prior)$mu/parameters(norm_prior)$sigma^2 +
                                        sum(filter(int_norm_df, trt == 0)$y)/sd_internal_control^2)

  ## Check that means and SDs of the normal posteriors are equal using both methods
  expect_equal(post_mean_beastt, post_mean_comp)
  expect_equal(post_sd_beastt, post_sd_comp)

  ## Check that a distribution is returned
  expect_s3_class(post_dist, "distribution")

})

##### Test for calc_post_norm() with:
#####    (1) t prior
#####    (2) internal_sd known
test_that("calc_post_norm returns the correct values for second case", {
  ## Define values to be used for both beastt code and comparison code
  sd_internal_control <- 0.15    # known SD for internal control
  t_prior <- dist_student_t(df = 10, mu = 0.5, sigma = 10)   # t prior

  ## beastt code
  post_dist <- calc_post_norm(filter(int_norm_df, trt == 0),
                              response = y,
                              prior = t_prior,
                              internal_sd = sd_internal_control)
  post_mean_beastt <- mean(post_dist)
  post_var_beastt <- variance(post_dist)

  ## Comparison code
  data_input <- list(
    N = nrow(filter(int_norm_df, trt == 0)),    # sample size of the internal control arm
    y = filter(int_norm_df, trt == 0)$y,        # N x 1 vector of binary responses (1: event, 0: no event)
    nu_pp = parameters(t_prior)$df,             # degrees of freedom for t prior component
    theta_pp = parameters(t_prior)$mu,          # mean hyperparameter for t prior component
    tau_pp = parameters(t_prior)$sigma,         # scale hyperparameter for t prior component
    theta_v = 0.5,      # mean hyperparameter for normal prior component
    tau_v = 10,         # standard deviation hyperparameter for normal prior component
    w = 1,              # prior weight associated with t prior component
    sigma_IC = sd_internal_control   # internal control SD (known)
  )
  set.seed(123)
  stan_fit <- sampling(stan_mod_sigma2_known, data = data_input, pars = "muC",
                       iter = 26000, warmup = 1000, chains = 4)
  post_draws <- as.matrix(stan_fit)           # posterior samples of each parameter
  post_mean_comp <- mean(post_draws[,"muC"])
  post_var_comp <- var(post_draws[,"muC"])

  ## Check that means and SDs of the normal posteriors are equal using both methods
  expect_equal(abs(post_mean_beastt-post_mean_comp) < 0.0001, TRUE)
  expect_equal(abs(post_var_beastt-post_var_comp) < 0.0001, TRUE)

  ## Check that a distribution is returned
  expect_s3_class(post_dist, "distribution")
})

##### Test for calc_post_norm() with:
#####    (1) mixture prior (t and normal)
#####    (2) internal_sd known
test_that("calc_post_norm returns the correct values for third case", {
  ## Define values to be used for both beastt code and comparison code
  sd_internal_control <- 0.15    # known SD for internal control
  norm_prior <- dist_normal(mu = 0.5, sigma = 10)   # normal prior
  t_prior <- dist_student_t(df = 10, mu = 0.5, sigma = 10)   # t prior
  mix_prior <- dist_mixture(t = t_prior, normal = norm_prior, weights = c(.5, .5))  # mixture prior

  ## beastt code
  post_dist <- calc_post_norm(filter(int_norm_df, trt == 0),
                              response = y,
                              prior = mix_prior,
                              internal_sd = sd_internal_control)
  post_mean_beastt <- mean(post_dist)
  post_var_beastt <- variance(post_dist)

  ## Comparison code
  data_input <- list(
    N = nrow(filter(int_norm_df, trt == 0)),    # sample size of the internal control arm
    y = filter(int_norm_df, trt == 0)$y,        # N x 1 vector of binary responses (1: event, 0: no event)
    nu_pp = parameters(t_prior)$df,             # degrees of freedom for t prior component
    theta_pp = parameters(t_prior)$mu,          # mean hyperparameter for t prior component
    tau_pp = parameters(t_prior)$sigma,         # scale hyperparameter for t prior component
    theta_v = parameters(t_prior)$mu,           # mean hyperparameter for normal prior component
    tau_v = parameters(t_prior)$sigma,          # standard deviation hyperparameter for normal prior component
    w = .5,                                     # prior weight associated with t prior component
    sigma_IC = sd_internal_control              # internal control SD (known)
  )
  set.seed(123)
  stan_fit <- sampling(stan_mod_sigma2_known, data = data_input, pars = "muC",
                       iter = 26000, warmup = 1000, chains = 4)
  post_draws <- as.matrix(stan_fit)           # posterior samples of each parameter
  post_mean_comp <- mean(post_draws[,"muC"])
  post_var_comp <- var(post_draws[,"muC"])

  ## Check that means and SDs of the normal posteriors are equal using both methods
  expect_equal(post_mean_beastt, post_mean_comp, tolerance=0.0002)
  # SDs very small so tolerance within expect_equal doesn't work
  expect_equal(abs(post_var_beastt-post_var_comp) < 0.0001, TRUE)

  ## Check that a distribution is returned
  expect_s3_class(post_dist, "distribution")
})

##### Test for calc_post_norm() with:
#####    (1) normal prior
#####    (2) internal_sd unknown
test_that("calc_post_norm returns the correct values for fourth case", {
  ## Define values to be used for both beastt code and comparison code
  norm_prior <- dist_normal(mu = 0.5, sigma = 10)   # normal prior

  ## beastt code
  post_dist <- calc_post_norm(filter(int_norm_df, trt == 0),
                              response = y,
                              prior = norm_prior,
                              internal_sd = NULL)
  post_mean_beastt <- mean(post_dist)
  post_var_beastt <- variance(post_dist)

  ## Comparison code
  data_input <- list(
    N = nrow(filter(int_norm_df, trt == 0)),    # sample size of the internal control arm
    y = filter(int_norm_df, trt == 0)$y,        # N x 1 vector of binary responses (1: event, 0: no event)
    nu_pp = 10,                                 # degrees of freedom for normal prior component
    theta_pp = .5,                              # mean hyperparameter for normal prior component
    tau_pp = 10,                                # scale hyperparameter for normal prior component
    theta_v = parameters(norm_prior)$mu,           # mean hyperparameter for normal prior component
    tau_v = parameters(norm_prior)$sigma,          # standard deviation hyperparameter for normal prior component
    w = .5                                      # prior weight associated with normal prior component
  )
  set.seed(123)
  stan_fit <- sampling(stan_mod_sigma2_unknown, data = data_input, pars = "muC",
                       iter = 26000, warmup = 1000, chains = 4)
  post_draws <- as.matrix(stan_fit)           # posterior samples of each parameter
  post_mean_comp <- mean(post_draws[,"muC"])
  post_var_comp <- var(post_draws[,"muC"])

  ## Check that means and SDs of the normal posteriors are equal using both methods
  expect_equal(abs(post_mean_beastt-post_mean_comp) < 0.0002, TRUE)
  expect_equal(abs(post_var_beastt-post_var_comp) < 0.0001, TRUE)

  ## Check that a distribution is returned
  expect_s3_class(post_dist, "distribution")
})

##### Test for calc_post_norm() with:
#####    (1) t prior
#####    (2) internal_sd unknown
test_that("calc_post_norm returns the correct values for fifth case", {
  ## Define values to be used for both beastt code and comparison code
  t_prior <- dist_student_t(df = 10, mu = 0.5, sigma = 10)   # t prior

  ## beastt code
  post_dist <- calc_post_norm(filter(int_norm_df, trt == 0),
                              response = y,
                              prior = t_prior,
                              internal_sd = NULL)
  post_mean_beastt <- mean(post_dist)
  post_var_beastt <- variance(post_dist)

  ## Comparison code
  data_input <- list(
    N = nrow(filter(int_norm_df, trt == 0)),    # sample size of the internal control arm
    y = filter(int_norm_df, trt == 0)$y,        # N x 1 vector of binary responses (1: event, 0: no event)
    nu_pp = parameters(t_prior)$df,             # degrees of freedom for t prior component
    theta_pp = parameters(t_prior)$mu,          # mean hyperparameter for t prior component
    tau_pp = parameters(t_prior)$sigma,         # scale hyperparameter for t prior component
    theta_v = 0.5,      # mean hyperparameter for normal prior component
    tau_v = 10,         # standard deviation hyperparameter for normal prior component
    w = 1               # prior weight associated with t prior component
  )
  set.seed(1234)
  stan_fit <- sampling(stan_mod_sigma2_unknown, data = data_input, pars = "muC",
                       iter = 26000, warmup = 1000, chains = 4)
  post_draws <- as.matrix(stan_fit)           # posterior samples of each parameter
  post_mean_comp <- mean(post_draws[,"muC"])
  post_var_comp <- var(post_draws[,"muC"])

  ## Check that means and SDs of the normal posteriors are equal using both methods
  expect_equal(post_mean_beastt, post_mean_comp, tolerance=0.0005)
  # SDs very small so tolerance within expect_equal doesn't work
  expect_equal(abs(post_var_beastt-post_var_comp) < 0.0001, TRUE)

  ## Check that a distribution is returned
  expect_s3_class(post_dist, "distribution")
})

##### Test for calc_post_norm() with:
#####    (1) mixture prior (t and normal)
#####    (2) internal_sd unknown
test_that("calc_post_norm returns the correct values for sixth case", {
  ## Define values to be used for both beastt code and comparison code
  norm_prior <- dist_normal(mu = 0.5, sigma = 10)   # normal prior
  t_prior <- dist_student_t(df = 10, mu = 0.5, sigma = 10)   # t prior
  mix_prior <- dist_mixture(t = t_prior, normal = norm_prior, weights = c(.5, .5))  # mixture prior

  ## beastt code
  post_dist <- calc_post_norm(filter(int_norm_df, trt == 0),
                              response = y,
                              prior = mix_prior,
                              internal_sd = NULL)
  post_mean_beastt <- mean(post_dist)
  post_var_beastt <- variance(post_dist)

  ## Comparison code
  data_input <- list(
    N = nrow(filter(int_norm_df, trt == 0)),    # sample size of the internal control arm
    y = filter(int_norm_df, trt == 0)$y,        # N x 1 vector of binary responses (1: event, 0: no event)
    nu_pp = parameters(t_prior)$df,             # degrees of freedom for t prior component
    theta_pp = parameters(t_prior)$mu,          # mean hyperparameter for t prior component
    tau_pp = parameters(t_prior)$sigma,         # scale hyperparameter for t prior component
    theta_v = parameters(t_prior)$mu,           # mean hyperparameter for normal prior component
    tau_v = parameters(t_prior)$sigma,          # standard deviation hyperparameter for normal prior component
    w = .5                                      # prior weight associated with t prior component
  )
  set.seed(123)
  stan_fit <- sampling(stan_mod_sigma2_unknown, data = data_input, pars = "muC",
                       iter = 26000, warmup = 1000, chains = 4)
  post_draws <- as.matrix(stan_fit)           # posterior samples of each parameter
  post_mean_comp <- mean(post_draws[,"muC"])
  post_var_comp <- var(post_draws[,"muC"])

  ## Check that means and SDs of the normal posteriors are equal using both methods
  expect_equal(abs(post_mean_beastt-post_mean_comp) < 0.0002, TRUE)
  expect_equal(abs(post_var_beastt-post_var_comp) < 0.0001, TRUE)

  ## Check that a distribution is returned
  expect_s3_class(post_dist, "distribution")
})

# Test for invalid internal data
test_that("calc_post_norm handles invalid internal data", {
  expect_error(calc_post_norm(c(5, 6),
                              response = y,
                              prior = pwr_prior <- dist_normal(0.71, 0.00035),
                              internal_sd = 0.15))
  expect_error(calc_post_norm("z",
                              response = y,
                              prior = pwr_prior <- dist_normal(0.71, 0.00035),
                              internal_sd = 0.15))
})

# Test for invalid prior
test_that("calc_post_norm handles invalid prior", {
  expect_error(calc_post_norm(filter(int_norm_df, trt == 0),
                              response = y,
                              prior = 5,
                              internal_sd = 0.15))
  expect_error(calc_post_norm(filter(int_norm_df, trt == 0),
                              response = y,
                              prior = "a",
                              internal_sd = 0.15))
  expect_error(calc_post_norm(filter(int_norm_df, trt == 0),
                              response = y,
                              prior = dist_beta(10, 15),
                              internal_sd = 0.15))
})

# Test for invalid response variable
test_that("calc_post_norm handles invalid response", {
  expect_error(calc_post_norm(filter(int_norm_df, trt == 0),
                              response = invalid,
                              prior = dist_normal(0.71, 0.00035),
                              internal_sd = 0.15))
})

# Test for invalid internal sd
test_that("calc_post_norm handles invalid internal sd", {
  expect_error(calc_post_norm(filter(int_norm_df, trt == 0),
                              response = y,
                              prior = dist_normal(0.71, 0.00035),
                              internal_sd = -2))
  expect_error(calc_post_norm(filter(int_norm_df, trt == 0),
                              response = y,
                              prior = dist_normal(0.71, 0.00035),
                              internal_sd = "c"))
})
