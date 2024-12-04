# Define test data
set.seed(1234)
internal_df <- data.frame(id_col = 1:20, cov1 = rnorm(10, 2), cov2 = rnorm(100, 20),
                          trt = rbinom(20, 1, 0.5), y = rnorm(20, 10, 0.5))
pwr_prior <- dist_normal(0.71, 0.00035)
sd_internal_control <- 0.15

# Test for valid inputs
test_that("calc_post_norm handles valid inputs", {
  post_norm_result <- calc_post_norm(filter(internal_df, trt == 0),
                                     response = y,
                                     prior = pwr_prior,
                                     internal_sd = sd_internal_control)
  expect_s3_class(post_norm_result, "distribution")
})

# Test for invalid internal data
test_that("calc_post_norm handles invalid internal data", {
  expect_error(calc_post_norm(c(5, 6),
                              response = y,
                              prior = pwr_prior,
                              internal_sd = sd_internal_control))
  expect_error(calc_post_norm("z",
                              response = y,
                              prior = pwr_prior,
                              internal_sd = sd_internal_control))
})

# Test for invalid prior
test_that("calc_post_norm handles invalid prior", {
  expect_error(calc_post_norm(filter(internal_df, trt == 0),
                              response = y,
                              prior = 5,
                              internal_sd = sd_internal_control))
  expect_error(calc_post_norm(filter(internal_df, trt == 0),
                              response = y,
                              prior = "a",
                              internal_sd = sd_internal_control))
})

# Test for invalid response variable
test_that("calc_post_norm handles invalid response", {
  expect_error(calc_post_norm(filter(internal_df, trt == 0),
                              response = invalid,
                              prior = pwr_prior,
                              internal_sd = sd_internal_control))
})

# Test for invalid internal sd
test_that("calc_post_norm handles invalid internal sd", {
  expect_error(calc_post_norm(filter(internal_df, trt == 0),
                              response = y,
                              prior = pwr_prior,
                              internal_sd = -2))
  expect_error(calc_post_norm(filter(internal_df, trt == 0),
                              response = y,
                              prior = pwr_prior,
                              internal_sd = "c"))
})
