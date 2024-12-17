# Define test data
set.seed(1234)
internal_df <- data.frame(id_col = 1:20, cov1 = rnorm(10, 2), cov2 = rnorm(100, 20),
                          trt = rbinom(20, 1, 0.5), y = rbinom(20, 1, 0.3))
external_df <- data.frame(id_col = 21:40, cov1 = rnorm(10, 2), cov2 = rnorm(100, 18),
                          trt = rep(0, 20), y = rbinom(20, 1, 0.5))
model <- as.formula("~cov1 + cov2")
initial_prior <- dist_beta(0.5, 0.5)

# Test for valid inputs
test_that("calc_power_prior_beta handles valid inputs", {
  ps_obj <- calc_prop_scr(internal_df, external_df, id_col, model)
  pwr_prior_beta_prop_scr_result <- calc_power_prior_beta(ps_obj,
                                                 response = y,
                                                 prior = initial_prior)
  expect_s3_class(pwr_prior_beta_prop_scr_result, "distribution")
  pwr_prior_beta_df_result <- calc_power_prior_beta(external_df,
                                                    response = y,
                                                    prior = initial_prior)
  expect_s3_class(pwr_prior_beta_df_result, "distribution")
})

# Test for values
test_that("calc_power_prior_beta returns the correct value", {
  pwr_prior_beta_result <- calc_power_prior_beta(external_df,
                                                 response = y,
                                                 prior = initial_prior)
  expect_equal(parameters(pwr_prior_beta_result)$shape1, 65.5)
  expect_equal(parameters(pwr_prior_beta_result)$shape2, 35.5)
  expect_equal(family(pwr_prior_beta_result), "beta")
})

# Test for invalid external data
test_that("calc_power_prior_beta handles invalid external data", {
  expect_error(calc_power_prior_beta(c(1, 2, 3), response=y, prior=initial_prior))
  expect_error(calc_power_prior_beta("abc", response=y, prior=initial_prior))
})

# Test for invalid prior
test_that("calc_power_prior_beta handles invalid prior", {
  expect_error(calc_power_prior_beta(external_df, response=y, prior=5))
  expect_error(calc_power_prior_beta(external_df, response=y, prior="a"))
})

# Test for invalid response variable
test_that("calc_power_prior_beta handles invalid response", {
  expect_error(calc_power_prior_beta(external_df, response=invalid, prior=initial_prior))
})
