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


test_that("Test highlight warning", {

  # Subset population
  binary_sim_rev <- binary_sim_df |>
    dplyr::filter(population == "no imbalance" & marg_trt_eff %in% c(0, .15)) |>
    dplyr::group_by(marg_trt_eff) |>
    dplyr::mutate(
      reject_H0_yes = dplyr::case_when(
        marg_trt_eff == 0 ~ reject_H0_yes,
        TRUE ~ rev(reject_H0_yes)
      ),
      no_borrowing_reject_H0_yes =
        dplyr::case_when(
          marg_trt_eff == 0 ~ no_borrowing_reject_H0_yes,
          TRUE ~ rev(no_borrowing_reject_H0_yes))

    )

  expect_warning(
    sweet_spot_plot(
      .data = binary_sim_rev,
      scenario_vars = c("marg_trt_eff"),
      trt_diff = marg_trt_eff,
      control_marg_param = true_control_RR,
      design_prior = pwr_prior,
      h0_prob = reject_H0_yes,
      h0_prob_no_borrowing = no_borrowing_reject_H0_yes
    )
  )

})

test_that("sweet_spot_plot errors with non-distributional prior column", {
  # Create test data with non-distributional prior column
  test_data <- data.frame(
    drift = rep(c(-0.1, 0, 0.1), each = 100),
    treatment_effect = rep(c(0, 0.1), each = 50, times = 3),
    control_param = runif(300, 0.1, 0.5),
    prob_accept_h0 = runif(300),
    prob_accept_h0_no_borrow = runif(300),
    prior_not_dist = runif(300)  # Non-distributional prior
  )

  # Should error because prior is not a distributional object
  expect_error(
    sweet_spot_plot(
      .data = test_data,
      scenario_vars = c(drift),
      trt_diff = treatment_effect,
      control_marg_param = control_param,
      design_prior = prior_not_dist,
      h0_prob = prob_accept_h0,
      h0_prob_no_borrowing = prob_accept_h0_no_borrow
    ),
    "`design_prior` must be a column of distributional objects"
  )
})

test_that("sweet_spot_plot can handel a null prior",{
  plots <- sweet_spot_plot(.data = binary_sim_df,
                           scenario_vars = c("population", "marg_trt_eff"),
                           trt_diff = marg_trt_eff,
                           control_marg_param = true_control_RR,
                           design_prior = NULL,
                           h0_prob = reject_H0_yes,
                           h0_prob_no_borrowing = no_borrowing_reject_H0_yes
  )

})

test_that("sweet_spot_plot errors with multivariate normal prior without approximation", {
  # Create test data with multivariate normal prior
  set.seed(123)
  mvn_prior <- replicate(300,
                         dist_multivariate_normal(
                           mu = list(c(0, -1)),
                           sigma = list(matrix(c(0.1, 0, 0, 0.1), nrow = 2))
                         ))
  class(mvn_prior) <- c("distribution", "vctrs_vctr", "list")

  test_data <- dplyr::tibble(
    drift = rep(c(-0.1, 0, 0.1), each = 100),
    treatment_effect = rep(c(0, 0.1), each = 50, times = 3),
    control_param = runif(300, 0.1, 0.5),
    prob_accept_h0 = runif(300),
    prob_accept_h0_no_borrow = runif(300),
    mvn_prior = mvn_prior
  )

  # Should error because multivariate normal priors need to be approximated
  expect_error(
    sweet_spot_plot(
      .data = test_data,
      scenario_vars = c(drift),
      trt_diff = treatment_effect,
      control_marg_param = control_param,
      design_prior = mvn_prior,
      h0_prob = prob_accept_h0,
      h0_prob_no_borrowing = prob_accept_h0_no_borrow
    )
  )
})

test_that("sweet_spot_plot handles case where no scenarios have trt_diff = 0", {
  # Create test data without any null scenarios
  binary_sim_no_0 <- binary_sim_df |>
    dplyr::filter(marg_trt_eff == .15)

  expect_error(
    sweet_spot_plot(
      .data = binary_sim_no_0,
      scenario_vars = c("marg_trt_eff"),
      trt_diff = marg_trt_eff,
      control_marg_param = true_control_RR,
      design_prior = pwr_prior,
      h0_prob = reject_H0_yes,
      h0_prob_no_borrowing = no_borrowing_reject_H0_yes
    ),
    "Unable to calculate Type 1 Error without a scenario where `trt_diff` equals 0"
  )
})

test_that("snapshot plot",{
  # Assuming binary_sim_df is a data frame with simulation results in the shape of binary template code  plots <- sweet_spot_plot(
  plots <- sweet_spot_plot(.data = binary_sim_df,
    scenario_vars = c("population", "marg_trt_eff"),
    trt_diff = marg_trt_eff,
    control_marg_param = true_control_RR,
    design_prior = pwr_prior,
    h0_prob = reject_H0_yes,
    h0_prob_no_borrowing = no_borrowing_reject_H0_yes
  )

  # test the binary plot
  vdiffr::expect_doppelganger("plot-sweet_spot_bin", plots[[1]])

  tte_plots <- tte_sim_df |>
   mutate(beta_appox = approx_mvn_at_time(mix_prior, time = 12)) |>
   sweet_spot_plot(
     scenario_vars = c("population", "marg_trt_eff"),
     trt_diff = marg_trt_eff,
     control_marg_param = true_control_surv_prob,
     design_prior = beta_appox,
     h0_prob = reject_H0_yes,
     h0_prob_no_borrowing = no_borrowing_reject_H0_yes
   )

  vdiffr::expect_doppelganger("plot-sweet_spot_tte", tte_plots[[1]])

  # Test without highlight
  plots_no_highlight <- sweet_spot_plot(.data = binary_sim_df,
                           scenario_vars = c("population", "marg_trt_eff"),
                           trt_diff = marg_trt_eff,
                           control_marg_param = true_control_RR,
                           design_prior = pwr_prior,
                           h0_prob = reject_H0_yes,
                           h0_prob_no_borrowing = no_borrowing_reject_H0_yes,
                           highlight = FALSE
  )

  # test the binary plot
  vdiffr::expect_doppelganger("plot-sweet_spot_bin_no_high", plots_no_highlight[[1]])

  # Test negative plot
  # reverse everything
  # Subset population
  binary_sim_all_rev <- binary_sim_df |>
    dplyr::filter(population == "no imbalance" & marg_trt_eff %in% c(0, .15)) |>
    dplyr::group_by(marg_trt_eff) |>
    dplyr::mutate(
      reject_H0_yes = rev(reject_H0_yes),
      no_borrowing_reject_H0_yes = rev(no_borrowing_reject_H0_yes))

  neg_plot <- sweet_spot_plot(
    .data = binary_sim_all_rev,
    scenario_vars = c("marg_trt_eff"),
    trt_diff = marg_trt_eff,
    control_marg_param = true_control_RR,
    design_prior = pwr_prior,
    h0_prob = reject_H0_yes,
    h0_prob_no_borrowing = no_borrowing_reject_H0_yes
  )

  vdiffr::expect_doppelganger("plot-sweet_spot_bin_neg_plot", neg_plot)
})

