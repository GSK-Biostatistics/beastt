
# simulate_accrual --------------------------------------------------------

test_that("simulate_accrual returns correct number of values", {
  n <- 100
  result <- simulate_accrual(n = n, accrual_periods = c(6, 8), accrual_props = c(0.5, 0.5))

  expect_equal(length(result), n)
  expect_true(is.numeric(result))
})

test_that("simulate_accrual returns values within expected range", {
  # Test with multiple periods
  result <- simulate_accrual(n = 1000,
                             accrual_periods = c(6, 8),
                             accrual_props = c(0.5, 0.5))

  expect_true(all(result >= 0))
  expect_true(all(result <= 8))

  # Values should be distributed across both periods
  period1_count <- sum(result < 6)
  period2_count <- sum(result >= 6)

  # Allow some flexibility due to randomness, but should be roughly 50/50
  expect_true(period1_count > 475)
  expect_true(period2_count > 475)
})

test_that("simulate_accrual works with single accrual period", {
  result <- simulate_accrual(n = 100,
                             accrual_periods = c(10),
                             accrual_props = c(1))

  expect_equal(length(result), 100)
  expect_true(all(result >= 0))
  expect_true(all(result <= 10))
})

test_that("simulate_accrual works with unequal accrual proportions", {
  set.seed(123) # For reproducibility
  result <- simulate_accrual(n = 1000,
                             accrual_periods = c(5, 10, 15),
                             accrual_props = c(0.2, 0.3, 0.5))

  # Roughly 20% should be in first period
  period1_count <- sum(result < 5)
  # Roughly 30% should be in second period
  period2_count <- sum(result >= 5 & result < 10)
  # Roughly 50% should be in third period
  period3_count <- sum(result >= 10)

  # Allow flexibility due to randomness
  expect_true(abs(period1_count/1000 - 0.2) < 0.05)
  expect_true(abs(period2_count/1000 - 0.3) < 0.05)
  expect_true(abs(period3_count/1000 - 0.5) < 0.05)
})


test_that("simulate_accrual handles extreme probability distributions", {
  # Test with very skewed probability
  result <- simulate_accrual(n = 1000,
                             accrual_periods = c(5, 10, 15),
                             accrual_props = c(0.01, 0.01, 0.98))

  # Most values should be in the third period (if perfect would be 980, there is wiggle room)
  expect_true(sum(result >= 10) > 950)
})


# simulate_tte_pwch -------------------------------------------------------


test_that("simulate_tte_pwch returns correct number and type of values", {
  n <- 100
  hazard_periods <- c(5, 10)
  hazard_values <- c(0.1, 0.2, 0.3)

  # Run the function
  set.seed(123)
  result <- simulate_tte_pwch(n = n, hazard_periods = hazard_periods, hazard_values = hazard_values)

  # Check results
  expect_equal(length(result), n)
  expect_true(is.numeric(result))
  expect_true(all(result >= 0))
})

test_that("simulate_tte_pwch handles constant hazard case", {
  n <- 1000
  hazard_value <- 0.1

  # Run with single hazard value (constant hazard)
  set.seed(123)
  result <- simulate_tte_pwch(n = n, hazard_values = hazard_value)

  # For exponential distribution with rate 0.1, mean should be approximately 10
  expect_gt(mean(result), 9)
  expect_lt(mean(result), 11)
})

test_that("simulate_tte_pwch handles piecewise constant hazard", {
  n <- 1000
  hazard_periods <- c(5)
  hazard_values <- c(0.05, 0.8)  # Low hazard then high hazard

  # Run the function
  set.seed(123)
  result <- simulate_tte_pwch(n = n, hazard_periods = hazard_periods, hazard_values = hazard_values)

  # Count events in different periods
  early_events <- sum(result <= 5)
  late_events <- sum(result > 5)

  # With higher hazard in second period, expect more events per unit time
  # Calculate event density (events per unit time)
  max_time <- max(result)
  early_density <- early_events / 5
  late_density <- late_events / (max_time - 5)

  # Higher hazard should lead to higher event density
  expect_gt(late_density, early_density)
})

test_that("simulate_tte_pwch handles edge cases", {
  # Very high hazard (should result in early events)
  high_hazard_result <- simulate_tte_pwch(n = 100, hazard_values = 10)
  expect_true(mean(high_hazard_result) < 1)  # Mean should be less than 1 with hazard=10

  # Very low hazard (should result in later events)
  low_hazard_result <- simulate_tte_pwch(n = 100, hazard_values = 0.01)
  expect_true(mean(low_hazard_result) > 50)  # Mean should be greater than 50 with hazard=0.01
})

test_that("simulate_tte_pwch validates inputs", {
  # Number of hazard values should match number of periods + 1
  expect_error(
    simulate_tte_pwch(n = 100, hazard_periods = c(5, 10), hazard_values = c(0.1, 0.2)),
    "`hazard_values` should have length equal to one more than the length of `hazard_periods`"
  )

  # Negative hazard values should error
  expect_error(
    simulate_tte_pwch(n = 100, hazard_values = -0.1),
    "`hazard_values` and `hazard_periods` can only have positive values"
  )
})


# simulate_tte_weib_ph ----------------------------------------------------


test_that("simulate_tte_weib_ph returns correct format", {
  # Create a simple survreg Weibull model
  set.seed(123)
  n <- 50
  df <- data.frame(
    time = rweibull(n, shape = 1.5, scale = 5),
    event = sample(0:1, n, replace = TRUE, prob = c(0.2, 0.8)),
    cov1 = rnorm(n),
    cov2 = rbinom(n, 1, 0.5)
  )
  weib_model <- survival::survreg(survival::Surv(time, event) ~ cov1 + cov2, data = df, dist = "weibull")

  # Create sample data frame
  samp_df <- data.frame(
    cov1 = rnorm(10),
    cov2 = rbinom(10, 1, 0.5)
  )

  # Run the function
  result <- simulate_tte_weib_ph(weibull_ph_mod = weib_model, samp_df = samp_df)

  # Check results
  expect_length(result, nrow(samp_df))
  expect_true(all(result > 0))
  expect_type(result, "double")
})

test_that("simulate_tte_weib_ph handles conditional effects correctly", {
  # Create a simple survreg Weibull model
  set.seed(123)
  n <- 50
  df <- data.frame(
    time = rweibull(n, shape = 1.5, scale = 5),
    event = sample(0:1, n, replace = TRUE, prob = c(0.2, 0.8)),
    cov1 = rnorm(n),
    cov2 = rbinom(n, 1, 0.5)
  )
  weib_model <- survival::survreg(survival::Surv(time, event) ~ cov1 + cov2, data = df, dist = "weibull")

  # Create sample data frame
  set.seed(456)
  samp_df <- data.frame(
    cov1 = rnorm(100),
    cov2 = rbinom(100, 1, 0.5)
  )

  # Generate times with no drift/treatment effect
  base_times <- simulate_tte_weib_ph(weibull_ph_mod = weib_model, samp_df = samp_df)

  # Generate times with positive drift (should increase survival time)
  pos_times <- simulate_tte_weib_ph(
    weibull_ph_mod = weib_model,
    samp_df = samp_df,
    cond_drift = 0.5
  )

  # Generate times with negative treatment effect (should decrease survival time)
  neg_times <- simulate_tte_weib_ph(
    weibull_ph_mod = weib_model,
    samp_df = samp_df,
    cond_trt_effect = -0.5
  )

  # Positive drift should increase survival time (higher median)
  expect_gt(median(neg_times), median(base_times))

  # Negative treatment effect should decrease survival time
  expect_lt(median(pos_times), median(base_times))
})



test_that("simulate_tte_weib_ph handles empty data correctly", {
  # Create a simple survreg Weibull model
  set.seed(123)
  n <- 50
  df <- data.frame(
    time = rweibull(n, shape = 1.5, scale = 5),
    event = sample(0:1, n, replace = TRUE, prob = c(0.2, 0.8)),
    cov1 = rnorm(n),
    cov2 = rbinom(n, 1, 0.5)
  )
  weib_model <- survival::survreg(survival::Surv(time, event) ~ cov1 + cov2, data = df, dist = "weibull")

  # Create empty sample data frame with correct columns
  samp_df <- data.frame(
    cov1 = numeric(0),
    cov2 = numeric(0)
  )

  # Function should return empty vector
  result <- simulate_tte_weib_ph(weibull_ph_mod = weib_model, samp_df = samp_df)
  expect_length(result, 0)
})


test_that("simulate_tte_weib_ph validates cond_drift input", {
  # Create a simple survreg Weibull model
  set.seed(123)
  n <- 50
  df <- data.frame(
    time = rweibull(n, shape = 1.5, scale = 5),
    event = sample(0:1, n, replace = TRUE, prob = c(0.2, 0.8)),
    cov1 = rnorm(n),
    cov2 = rbinom(n, 1, 0.5)
  )
  weib_model <- survival::survreg(survival::Surv(time, event) ~ cov1 + cov2, data = df, dist = "weibull")

  # Sample data
  samp_df <- data.frame(
    cov1 = rnorm(10),
    cov2 = rbinom(10, 1, 0.5)
  )

  # Test error when cond_drift is not numeric
  expect_error(
    simulate_tte_weib_ph(weib_model, samp_df, cond_drift = "not a number"),
    "`cond_drift` must be a single number"
  )

  # Test error when cond_drift is a vector
  expect_error(
    simulate_tte_weib_ph(weib_model, samp_df, cond_drift = c(0.1, 0.2)),
    "`cond_drift` must be a single number"
  )
})

test_that("simulate_tte_weib_ph validates cond_trt_effect input", {
  # Create a simple survreg Weibull model
  set.seed(123)
  n <- 50
  df <- data.frame(
    time = rweibull(n, shape = 1.5, scale = 5),
    event = sample(0:1, n, replace = TRUE, prob = c(0.2, 0.8)),
    cov1 = rnorm(n),
    cov2 = rbinom(n, 1, 0.5)
  )
  weib_model <- survival::survreg(survival::Surv(time, event) ~ cov1 + cov2, data = df, dist = "weibull")

  # Sample data
  samp_df <- data.frame(
    cov1 = rnorm(10),
    cov2 = rbinom(10, 1, 0.5)
  )

  # Test error when cond_trt_effect is not numeric
  expect_error(
    simulate_tte_weib_ph(weib_model, samp_df, cond_trt_effect = "not a number"),
    "`cond_trt_effect` must be a single number"
  )

  # Test error when cond_trt_effect is a vector
  expect_error(
    simulate_tte_weib_ph(weib_model, samp_df, cond_trt_effect = c(0.1, 0.2)),
    "`cond_trt_effect` must be a single number"
  )
})

test_that("simulate_tte_weib_ph validates weibull_ph_mod input", {
  # Create sample data
  samp_df <- data.frame(
    cov1 = rnorm(10),
    cov2 = rbinom(10, 1, 0.5)
  )

  # Test error with non-survreg object
  expect_error(
    simulate_tte_weib_ph(weibull_ph_mod = "not a survreg object", samp_df = samp_df),
    "`weibull_ph_mod` must be a `survreg` object with a Wiebull distribution"
  )

  # Test error with non-Weibull survreg object
  set.seed(123)
  n <- 50
  df <- data.frame(
    time = rexp(n, 0.1),
    event = sample(0:1, n, replace = TRUE, prob = c(0.2, 0.8)),
    cov1 = rnorm(n),
    cov2 = rbinom(n, 1, 0.5)
  )

  # Create an exponential model (not Weibull)
  exp_model <- survival::survreg(survival::Surv(time, event) ~ cov1 + cov2, data = df, dist = "exponential")

  # This should error because it's not a Weibull model
  expect_error(
    simulate_tte_weib_ph(weibull_ph_mod = exp_model, samp_df = samp_df),
    "`weibull_ph_mod` must be a `survreg` object with a Wiebull distribution"
  )
})

test_that("simulate_tte_weib_ph validates cond_drift input", {
  # Create a simple survreg Weibull model
  set.seed(123)
  n <- 50
  df <- data.frame(
    time = rweibull(n, shape = 1.5, scale = 5),
    event = sample(0:1, n, replace = TRUE, prob = c(0.2, 0.8)),
    cov1 = rnorm(n),
    cov2 = rbinom(n, 1, 0.5)
  )
  weib_model <- survival::survreg(survival::Surv(time, event) ~ cov1 + cov2, data = df, dist = "weibull")

  # Sample data
  samp_df <- data.frame(
    cov1 = rnorm(10),
    cov2 = rbinom(10, 1, 0.5)
  )

  # Test error when cond_drift is not numeric
  expect_error(
    simulate_tte_weib_ph(weib_model, samp_df, cond_drift = "not a number"),
    "`cond_drift` must be a single number"
  )

  # Test error when cond_drift is a vector
  expect_error(
    simulate_tte_weib_ph(weib_model, samp_df, cond_drift = c(0.1, 0.2)),
    "`cond_drift` must be a single number"
  )
})

test_that("simulate_tte_weib_ph validates cond_trt_effect input", {
  # Create a simple survreg Weibull model
  set.seed(123)
  n <- 50
  df <- data.frame(
    time = rweibull(n, shape = 1.5, scale = 5),
    event = sample(0:1, n, replace = TRUE, prob = c(0.2, 0.8)),
    cov1 = rnorm(n),
    cov2 = rbinom(n, 1, 0.5)
  )
  weib_model <- survival::survreg(survival::Surv(time, event) ~ cov1 + cov2, data = df, dist = "weibull")

  # Sample data
  samp_df <- data.frame(
    cov1 = rnorm(10),
    cov2 = rbinom(10, 1, 0.5)
  )

  # Test error when cond_trt_effect is not numeric
  expect_error(
    simulate_tte_weib_ph(weib_model, samp_df, cond_trt_effect = "not a number"),
    "`cond_trt_effect` must be a single number"
  )

  # Test error when cond_trt_effect is a vector
  expect_error(
    simulate_tte_weib_ph(weib_model, samp_df, cond_trt_effect = c(0.1, 0.2)),
    "`cond_trt_effect` must be a single number"
  )
})

test_that("simulate_tte_weib_ph validates weibull_ph_mod input", {
  # Create sample data
  samp_df <- data.frame(
    cov1 = rnorm(10),
    cov2 = rbinom(10, 1, 0.5)
  )

  # Test error with non-survreg object
  expect_error(
    simulate_tte_weib_ph(weibull_ph_mod = "not a survreg object", samp_df = samp_df),
    "`weibull_ph_mod` must be a `survreg` object with a Wiebull distribution"
  )

  # Test error with non-Weibull survreg object
  set.seed(123)
  n <- 50
  df <- data.frame(
    time = rexp(n, 0.1),
    event = sample(0:1, n, replace = TRUE, prob = c(0.2, 0.8)),
    cov1 = rnorm(n),
    cov2 = rbinom(n, 1, 0.5)
  )

  # Create an exponential model (not Weibull)
  exp_model <- survival::survreg(survival::Surv(time, event) ~ cov1 + cov2, data = df, dist = "exponential")

  # This should error because it's not a Weibull model
  expect_error(
    simulate_tte_weib_ph(weibull_ph_mod = exp_model, samp_df = samp_df),
    "`weibull_ph_mod` must be a `survreg` object with a Wiebull distribution"
  )
})


