# Define test data
set.seed(1234)
internal_df <- data.frame(id_col = 1:20, cov1 = rnorm(10, 2), cov2 = rnorm(100, 20),
                          trt = rbinom(20, 1, 0.5), y = rnorm(20, 10, 0.5))
external_df <- data.frame(id_col = 21:40, cov1 = rnorm(10, 2), cov2 = rnorm(100, 18),
                          trt = rep(0, 20), y = rnorm(20, 8, 0.7))
model <- as.formula("~cov1 + cov2")
initial_prior <- dist_normal(50, 10)

# Test for valid inputs
test_that("calc_power_prior_norm handles valid inputs", {
  ps_obj <- calc_prop_scr(internal_df, external_df, id_col, model)
  pwr_prior_norm_prop_scr_result <- calc_power_prior_norm(ps_obj,
                                                          response = y,
                                                          prior = initial_prior,
                                                          external_sd = sd_external_control)
  expect_s3_class(pwr_prior_norm_prop_scr_result, "distribution")
  pwr_prior_norm_df_result <- calc_power_prior_norm(ps_obj,
                                                    response = y,
                                                    prior = initial_prior,
                                                    external_sd = sd_external_control)
  expect_s3_class(pwr_prior_norm_df_result, "distribution")
})

# Test for invalid external data
test_that("calc_power_prior_norm handles invalid external data", {
  expect_error(calc_power_prior_norm(c(1, 2, 3),
                                     response=y,
                                     prior=initial_prior,
                                     external_sd = sd_external_control))
  expect_error(calc_power_prior_norm("abc",
                                     response=y,
                                     prior=initial_prior,
                                     external_sd = sd_external_control))
})

# Test for invalid prior
test_that("calc_power_prior_norm handles invalid prior", {
  expect_error(calc_power_prior_norm(external_df,
                                     response=y,
                                     prior=5,
                                     external_sd = sd_external_control))
  expect_error(calc_power_prior_norm(external_df,
                                     response=y,
                                     prior="a",
                                     external_sd = sd_external_control))
})

# Test for invalid response variable
test_that("calc_power_prior_norm handles invalid response", {
  expect_error(calc_power_prior_norm(external_df,
                                     response=invalid,
                                     prior=initial_prior,
                                     external_sd = sd_external_control))
})

# Test for invalid external sd
test_that("calc_power_prior_norm handles invalid external sd", {
  expect_error(calc_power_prior_norm(external_df,
                                     response=y,
                                     prior=initial_prior,
                                     external_sd = -5))
  expect_error(calc_power_prior_norm(external_df,
                                     response=y,
                                     prior=initial_prior,
                                     external_sd = "a"))
})
