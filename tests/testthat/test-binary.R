###################################################################################################
# Code to test various functions of beastt - binary endpoint
###################################################################################################

###################################################################################################
# calc_power_prior_beta
###################################################################################################

##### Test for calc_power_prior_beta() with prop_scr_obj (not unweighted external data)
test_that("calc_power_prior_beta returns correct values when using prop_scr_obj", {
  ## Define values to be used for both beastt code and comparison code
  init_prior <- dist_beta(shape1 = 0.5, shape2 = 0.5)   # proper initial prior for power prior
  ps_obj <- calc_prop_scr(internal_df = filter(int_binary_df, trt == 0),
                          external_df = ex_binary_df,
                          id_col = subjid,
                          model = ~ cov1 + cov2 + cov3 + cov4)

  ## beastt code
  pwr_prior <- calc_power_prior_beta(external_data = ps_obj,
                                     response = y,
                                     prior = init_prior)
  pp_shape1_beastt <- parameters(pwr_prior)$shape1
  pp_shape2_beastt <- parameters(pwr_prior)$shape2

  ## Comparison code
  ipws <- ps_obj$external_df$`___weight___`
  pp_shape1_comp <- sum(ipws * ex_binary_df$y) + parameters(init_prior)$shape1
  pp_shape2_comp <- sum(ipws) - sum(ipws * ex_binary_df$y) + parameters(init_prior)$shape2

  ## Check that the shape 1 and shape 2 parameters of the beta power priors are equal using both methods
  expect_equal(pp_shape1_beastt, pp_shape1_comp)
  expect_equal(pp_shape2_beastt, pp_shape2_comp)

})

##### Test for calc_power_prior_beta() with unweighted external data (read in external df, not prop_scr_obj)
test_that("calc_power_prior_beta returns correct values when using external data", {
  ## Define values to be used for both beastt code and comparison code
  init_prior <- dist_beta(shape1 = 0.5, shape2 = 0.5)   # proper initial prior for power prior

  ## beastt code
  pwr_prior <- calc_power_prior_beta(external_data = ex_binary_df,
                                     response = y,
                                     prior = init_prior)
  pp_shape1_beastt <- parameters(pwr_prior)$shape1
  pp_shape2_beastt <- parameters(pwr_prior)$shape2

  ## Comparison code
  pp_shape1_comp <- sum(ex_binary_df$y) + parameters(init_prior)$shape1
  pp_shape2_comp <- nrow(ex_binary_df) - sum(ex_binary_df$y) + parameters(init_prior)$shape2

  ## Check that the shape 1 and shape 2 parameters of the beta power priors are equal using both methods
  expect_equal(pp_shape1_beastt, pp_shape1_comp)
  expect_equal(pp_shape2_beastt, pp_shape2_comp)

  # Check that a distribution is returned
  expect_s3_class(pwr_prior, "distribution")

})

# Test for invalid external data
test_that("calc_power_prior_beta handles invalid external data", {
  expect_error(calc_power_prior_beta(c(1, 2, 3), response=y, prior=dist_beta(0.5, 0.5)))
  expect_error(calc_power_prior_beta("abc", response=y, prior=dist_beta(0.5, 0.5)))
})

# Test for invalid prior
test_that("calc_power_prior_beta handles invalid prior", {
  expect_error(calc_power_prior_beta(ex_binary_df, response=y, prior=5))
  expect_error(calc_power_prior_beta(ex_binary_df, response=y, prior="a"))
  expect_error(calc_power_prior_beta(ex_binary_df, response=y, prior=dist_norm(5, 1)))
})

# Test for invalid response variable
test_that("calc_power_prior_beta handles invalid response", {
  expect_error(calc_power_prior_beta(ex_binary_df, response=invalid, prior=dist_beta(0.5, 0.5)))
})


###################################################################################################
# calc_posterior_beta
###################################################################################################

##### Test for calc_posterior_beta() with a single beta prior
test_that("calc_posterior_beta returns correct values with single beta prior", {
  ## Define values to be used for both beastt code and comparison code
  beta_prior <- dist_beta(shape1 = 0.5, shape2 = 0.5)   # beta prior

  ## beastt code
  post_dist <- calc_post_beta(filter(int_binary_df, trt == 0),
                              response = y,
                              prior = beta_prior)
  post_shape1_beastt <- parameters(post_dist)$shape1
  post_shape2_beastt <- parameters(post_dist)$shape2

  ## Comparison code
  post_shape1_comp <- sum(filter(int_binary_df, trt == 0)$y) + parameters(beta_prior)$shape1
  post_shape2_comp <- nrow(filter(int_binary_df, trt == 0)) - sum(filter(int_binary_df, trt == 0)$y) +
    parameters(beta_prior)$shape2

  ## Check that the shape 1 and shape 2 parameters of the beta posteriors are equal using both methods
  expect_equal(post_shape1_beastt, post_shape1_comp)
  expect_equal(post_shape2_beastt, post_shape2_comp)

  # Check that a distribution is returned
  expect_s3_class(post_dist, "distribution")
})

##### Test for calc_postior_beta() with a mixture prior with two beta components
test_that("calc_posterior_beta returns correct values with mixture prior with two beta components", {
  ## Define values to be used for both beastt code and comparison code
  beta1_shape1 <- .5
  beta1_shape2 <- .5
  beta2_shape1 <- 3
  beta2_shape2 <- 2
  w_beta1 <- .5
  beta_prior1 <- dist_beta(shape1 = beta1_shape1, shape2 = beta1_shape2)    # beta prior 1
  beta_prior2 <- dist_beta(shape1 = beta2_shape1, shape2 = beta2_shape2)    # beta prior 2
  mix_prior <- dist_mixture(beta1 = beta_prior1, beta2 = beta_prior2, weights = c(w_beta1, 1 - w_beta1))

  ## beastt code
  post_dist <- calc_post_beta(filter(int_binary_df, trt == 0),
                              response = y,
                              prior = mix_prior)
  post_mean_beastt <- mean(post_dist)
  post_var_beastt <- variance(post_dist)

  ## Comparison code
  mix_prior_rbest <- RBesT::mixbeta(beta1 = c(w_beta1, beta1_shape1, beta1_shape2),
                             beta2 = c(1 - w_beta1, beta2_shape1, beta2_shape2))
  post_beta_mix <- RBesT::postmix(mix_prior_rbest, n = nrow(filter(int_binary_df, trt == 0)),
                           r = sum(filter(int_binary_df, trt == 0)$y))
  post_beta1_w_comp <- post_beta_mix["w", "beta1"]
  post_beta1_shape1_comp <- post_beta_mix["a", "beta1"]
  post_beta1_shape2_comp <- post_beta_mix["b", "beta1"]
  post_beta2_w_comp <- post_beta_mix["w", "beta2"]
  post_beta2_shape1_comp <- post_beta_mix["a", "beta2"]
  post_beta2_shape2_comp <- post_beta_mix["b", "beta2"]
  post_mean1_comp <- post_beta1_shape1_comp / (post_beta1_shape1_comp + post_beta1_shape2_comp)
  post_var1_comp <- post_beta1_shape1_comp * post_beta1_shape2_comp /
    ((post_beta1_shape1_comp + post_beta1_shape2_comp)^2 * (post_beta1_shape1_comp + post_beta1_shape2_comp + 1))
  post_mean2_comp <- post_beta2_shape1_comp / (post_beta2_shape1_comp + post_beta2_shape2_comp)
  post_var2_comp <- post_beta2_shape1_comp * post_beta2_shape2_comp /
    ((post_beta2_shape1_comp + post_beta2_shape2_comp)^2 * (post_beta2_shape1_comp + post_beta2_shape2_comp + 1))
  post_mean_comp <- post_beta1_w_comp * post_mean1_comp + post_beta2_w_comp * post_mean2_comp
  post_var_comp <- post_beta1_w_comp * post_var1_comp + post_beta2_w_comp * post_var2_comp

  ## Check that the means and variances of the beta posteriors are equal using both methods
  expect_equal(abs(post_mean_beastt - post_mean_comp) < 0.000001, TRUE)
  expect_equal(abs(post_var_beastt - post_var_comp) < 0.00001, TRUE)

  # Check that a distribution is returned
  expect_s3_class(post_dist, "distribution")
})

# Test for invalid internal data
test_that("calc_post_beta handles invalid internal data", {
  expect_error(calc_post_beta(c(4, 8, 12), response = y, prior = dist_beta(38, 44)))
  expect_error(calc_post_beta("abc", response = y, prior = dist_beta(38, 44)))
})

# Test for invalid prior
test_that("calc_post_beta handles invalid prior", {
  expect_error(calc_post_beta(filter(int_binary_df, trt==0), response = y, prior = 5))
  expect_error(calc_post_beta(filter(int_binary_df, trt==0), response = y, prior = "a"))
  expect_error(calc_post_beta(filter(int_binary_df, trt==0), response = y, prior = dist_norm(5,1)))
})

# Test for invalid response variable
test_that("calc_post_beta handles invalid response", {
  expect_error(calc_post_beta(filter(int_binary_df, trt==0), response = invalid, prior = dist_beta(38, 44)))
})

