# Define test data
set.seed(1234)
internal_df <- data.frame(id_col = 1:20, cov1 = rnorm(10, 2), cov2 = rnorm(100, 20),
                          trt = rbinom(20, 1, 0.5), y = rbinom(20, 1, 0.3))
pwr_prior <- dist_beta(38, 44)

# Test for valid inputs
test_that("calc_post_beta handles valid inputs", {
  post_beta_result <- calc_post_beta(filter(internal_df, trt==0),
                                     response = y,
                                     prior = pwr_prior)
  expect_s3_class(post_beta_result, "distribution")
})

# Test for values
test_that("calc_post_beta returns the correct value", {
  post_beta_result <- calc_post_beta(filter(internal_df, trt==0),
                                     response = y,
                                     prior = pwr_prior)
  expect_equal(parameters(post_beta_result)$shape1, 48)
  expect_equal(parameters(post_beta_result)$shape2, 94)
  expect_equal(family(post_beta_result), "beta")
})

# Test for invalid external data
test_that("calc_post_beta handles invalid internal data", {
  expect_error(calc_post_beta(c(4, 8, 12), response = y, prior = pwr_prior))
  expect_error(calc_post_beta("abc", response = y, prior = pwr_prior))
})

# Test for invalid prior
test_that("calc_post_beta handles invalid prior", {
  expect_error(calc_post_beta(filter(internal_df, trt==0), response = y, prior = 5))
  expect_error(calc_post_beta(filter(internal_df, trt==0), response = y, prior = "a"))
})

# Test for invalid response variable
test_that("calc_post_beta handles invalid response", {
  expect_error(calc_post_beta(filter(internal_df, trt==0), response = invalid, prior = pwr_prior))
})
