# Test for valid inputs
test_that("plot_dist_mvnorm handles valid inputs", {
  foo <- list(dist_multivariate_normal(mu = list(c(1, 3)),
                                       sigma = list(matrix(c(4, 2, 2, 3)))))
  plot_dist_result <- plot_dist_mvnorm(foo)
  expect_s3_class(plot_dist_result, "ggplot")
})

# Test for correct output
test_that("plot_dist_mvnorm produces correct output", {
  plot_dist_test <- plot_dist_mvnorm(dist_multivariate_normal(mu = list(c(1, 3)),
                                                              sigma = list(matrix(c(4, 2, 2, 3)))))
  vdiffr::expect_doppelganger("plot-dist-test", plot_dist_test)
})

# Test for invalid inputs
test_that("plot_dist_mvnorm handles invalid inputs", {
  expect_error(plot_dist_mvnorm("a"))
  expect_error(plot_dist_mvnorm(5))
  expect_error(plot_dist_mvnorm(dist_beta(5, 7)))
  expect_error(plot_dist_mvnorm(dist_multivariate_normal(mu = 0, sigma = diag(1))))
})
