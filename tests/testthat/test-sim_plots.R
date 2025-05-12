test_that("avg_dist correctly averages beta distributions", {
  # Create test beta distributions
  beta_dists <- c(
    distributional::dist_beta(shape1 = 2, shape2 = 4),
    distributional::dist_beta(shape1 = 4, shape2 = 2),
    distributional::dist_beta(shape1 = 3, shape2 = 3)
  )

  # Calculate average
  avg_beta <- avg_dist(beta_dists)

  # Extract parameters
  params <- distributional::parameters(avg_beta)

  # Check correct distribution family
  expect_equal(family(avg_beta), "beta")

  # Check parameter averaging is correct
  expect_equal(params$shape1, 3)
  expect_equal(params$shape2, 3)
})

test_that("avg_dist correctly averages normal distributions", {
  # Create test normal distributions
  norm_dists <- c(
    distributional::dist_normal(mu = 0, sigma = 1),
    distributional::dist_normal(mu = 2, sigma = 3),
    distributional::dist_normal(mu = -2, sigma = 2)
  )

  # Calculate average
  avg_norm <- avg_dist(norm_dists)

  # Extract parameters
  params <- distributional::parameters(avg_norm)

  # Check correct distribution family
  expect_equal(family(avg_norm), "normal")

  # Check parameter averaging is correct
  expect_equal(params$mu, 0)
  expect_equal(params$sigma, 2)
})

test_that("avg_dist correctly averages multivariate normal distributions", {
  # Create test multivariate normal distributions
  mvn_dists <- c(
    distributional::dist_multivariate_normal(
      mu = list(c(0, 0)),
      sigma = list(matrix(c(1, 0, 0, 1), nrow = 2))
    ),
    distributional::dist_multivariate_normal(
      mu = list(c(2, 2)),
      sigma = list(matrix(c(3, 1, 1, 3), nrow = 2))
    )
  )

  # Calculate average
  avg_mvn <- avg_dist(mvn_dists)

  # Extract parameters
  params <- parameters(avg_mvn)

  # Check correct distribution family
  expect_equal(family(avg_mvn), "mvnorm")

  # Check parameter averaging is correct
  expect_equal(params$mu, list(c(1, 1)))

  # Expected average of covariance matrices
  expected_sigma <- list(matrix(c(2, 0.5, 0.5, 2), nrow = 2))
  expect_equal(params$sigma, expected_sigma)
})

test_that("avg_dist handles single-element vectors", {
  # Single beta distribution
  single_beta <- distributional::dist_beta(shape1 = 2, shape2 = 5)
  avg_single_beta <- avg_dist(c(single_beta))

  # Should return identical distribution
  expect_equal(distributional::parameters(avg_single_beta)$shape1, 2)
  expect_equal(distributional::parameters(avg_single_beta)$shape2, 5)
})

test_that("avg_dist produces error for invalid inputs", {
  # Not a distribution object
  expect_error(
    avg_dist(c(1, 2, 3)),
    "`x` must be a distributional object"
  )

  # Mixed distribution types
  mixed_dists <- c(
    distributional::dist_beta(shape1 = 2, shape2 = 3),
    distributional::dist_normal(mu = 0, sigma = 1)
  )
  expect_error(
    avg_dist(mixed_dists),
    "All distributional objects in x must be of the same type"
  )

  # Unsupported distribution type
  poisson_dist <- c(distributional::dist_poisson(lambda = 3),
                    distributional::dist_poisson(lambda = 5))
  expect_error(
    avg_dist(poisson_dist))
})


test_that("approx_mvn_at_time correctly converts multivariate normal to beta", {
  # Create a simple multivariate normal distribution
  mvn_dist <- dist_multivariate_normal(
    mu = list(c(0, -1)),  # log(shape) = 0, log(scale) = -1
    sigma = list(matrix(c(0.1, 0, 0, 0.1), nrow = 2))
  )

  # Convert to beta at a specific time
  set.seed(123)
  beta_approx <- approx_mvn_at_time(mvn_dist, time = 12)

  # Check output type
  expect_s3_class(beta_approx, "distribution")
  expect_equal(family(beta_approx), "beta")

  # Check parameters are reasonable
  params <- parameters(beta_approx)
  expect_true(is.numeric(params$shape1))
  expect_true(is.numeric(params$shape2))
  expect_true(params$shape1 > 0)
  expect_true(params$shape2 > 0)

  # For this specific distribution at t=12, we expect survival probabilities
  # around the exp(-(12*exp(-1))^exp(0)) ≈ exp(-(12*0.368)^1) ≈ 0.012
  # So the beta distribution should be heavily concentrated near 0
  mean_survival <- params$shape1 / (params$shape1 + params$shape2)
  expect_lt(mean_survival, 0.1)  # Mean should be quite small
})

test_that("approx_mvn_at_time correctly handles vectors of distributions", {
  # Create two multivariate normal distributions
  mvn_dist1 <- dist_multivariate_normal(
    mu = list(c(0, -1)),  # log(shape) = 0, log(scale) = -1
    sigma = list(matrix(c(0.1, 0, 0, 0.1), nrow = 2))
  )

  mvn_dist2 <- dist_multivariate_normal(
    mu = list(c(0.5, -0.5)),  # log(shape) = 0.5, log(scale) = -0.5
    sigma = list(matrix(c(0.1, 0, 0, 0.1), nrow = 2))
  )

  # Combine into a vector
  mvn_vector <- c(mvn_dist1, mvn_dist2)

  # Convert to beta at a specific time
  set.seed(123)
  beta_vector <- approx_mvn_at_time(mvn_vector, time = 12)

  # Check output type
  expect_s3_class(beta_vector, "distribution")
  expect_equal(length(beta_vector), 2)
  expect_equal(unique(family(beta_vector)), "beta")

})

test_that("approx_mvn_at_time correctly handles mixture distributions", {
  # Create two multivariate normal distributions
  mvn_dist1 <- dist_multivariate_normal(
    mu = list(c(0, -1)),
    sigma = list(matrix(c(0.1, 0, 0, 0.1), nrow = 2))
  )

  mvn_dist2 <- dist_multivariate_normal(
    mu = list(c(0.5, -0.5)),
    sigma = list(matrix(c(0.1, 0, 0, 0.1), nrow = 2))
  )

  # Create a mixture distribution
  mix_dist <- dist_mixture(mvn_dist1, mvn_dist2, weights = c(0.7, 0.3))

  # Convert to beta at a specific time
  set.seed(123)
  beta_mix <- approx_mvn_at_time(mix_dist, time = 12)

  # Check output type
  expect_s3_class(beta_mix, "distribution")
  expect_equal(family(beta_mix), "mixture")

  # Check the mixture has the expected components
  mix_params <- parameters(beta_mix)
  expect_length(mix_params$dist[[1]], 2)  # Should have 2 component distributions
  expect_equal(mix_params$w[[1]], c(0.7, 0.3))  # Weights should be preserved


})

