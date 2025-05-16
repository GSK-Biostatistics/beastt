## Tests for plot_dist (and plot_dist_mvnorm) ##

# Test for valid inputs
test_that("plot_dist handles valid inputs", {
  plot_dist_result <- plot_dist(dist_normal(10, 2), dist_normal(5, 1))
  expect_s3_class(plot_dist_result, "ggplot")
})

test_that("plot_dist handles valid inputs for mvnorm", {
  plot_dist_mvn_result <- plot_dist(dist_multivariate_normal(mu = list(c(1, 2)),
                                                             sigma = list(matrix(c(4, 2, 2, 3), ncol=2))))
  expect_s3_class(plot_dist_mvn_result, "ggplot")
})

# Test for correct output
test_that("plot_dist produces correct output for one distribution", {
  plot_dist_one <- plot_dist(dist_normal(10, 2))
  vdiffr::expect_doppelganger("plot-dist-one", plot_dist_one)
})

test_that("plot_dist produces correct output for multiple distributions", {
  plot_dist_mult <- plot_dist(dist_normal(10, 2), dist_normal(5, 1))
  vdiffr::expect_doppelganger("plot-dist-mult", plot_dist_mult)
})

test_that("plot_dist produces correct output for duplicate distributions", {
  plot_dist_dup <- plot_dist(dist_normal(10, 2), dist_normal(10, 2))
  vdiffr::expect_doppelganger("plot-dist-dup", plot_dist_dup)
})

test_that("plot_dist produces correct output for one mvnorm", {
  plot_dist_mvn <- plot_dist(dist_multivariate_normal(mu = list(c(1, 2)),
                                                      sigma = list(matrix(c(4, 2, 2, 3), ncol=2))))
  vdiffr::expect_doppelganger("plot-dist-mvn", plot_dist_mvn)
})

test_that("plot_dist produces correct output for multiple mvnorm", {
  plot_dist_mvn2 <- plot_dist(dist_multivariate_normal(mu = list(c(1, 2)),
                                                       sigma = list(matrix(c(4, 2, 2, 3), ncol=2))),
                              dist_multivariate_normal(mu = list(c(3, 4)),
                                                       sigma = list(matrix(c(4, 2, 2, 3), ncol=2))))
  vdiffr::expect_doppelganger("plot-dist-mvn2", plot_dist_mvn2)
})

# Test for invalid inputs
test_that("plot_dist handles invalid inputs", {
  expect_error(plot_dist("a"))
  expect_error(plot_dist(5))
  expect_error(plot_dist(dist_multivariate_normal(mu = list(0), sigma = diag(1))))
})

## Tests for is_mvnorm ##

# Test for valid inputs
test_that("is_mvnorm handles inputs", {
  expect_equal(is_mvnorm(dist_beta(5, 3)), FALSE)
  expect_equal(is_mvnorm(dist_multivariate_normal(mu = list(c(1, 2)),
                                                  sigma = list(matrix(c(4, 2, 2, 3), ncol=2)))), TRUE)
})

## Tests for robustify_norm (and robustify_mvnorm) ##

# Test for valid inputs
test_that("robustify_norm handles valid inputs", {
  robustify_norm_test <- robustify_norm(dist_normal(0, 1), n = 15)
  expect_s3_class(robustify_norm_test, "distribution")

  robustify_mvnorm_test <- robustify_norm(dist_multivariate_normal(mu = list(c(1, 0)),
                                                                   sigma = list(c(10, 5))),
                                          n = 15)
  expect_s3_class(robustify_mvnorm_test, "distribution")
})

# Test for invalid inputs
test_that("robustify_norm handles invalid inputs", {
  expect_error(robustify_norm(dist_beta(5, 3), n = 15))
  expect_error(robustify_norm(dist_normal(0, 1), n = -5))
  expect_error(robustify_norm(dist_normal(0, 1), weights = c(0.4, 0.8)))
})
