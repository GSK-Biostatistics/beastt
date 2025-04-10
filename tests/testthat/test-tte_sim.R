
# sim_accrual --------------------------------------------------------

test_that("sim_accrual returns correct number of values", {
  n <- 100
  result <- sim_accrual(n = n, accrual_periods = c(6, 8), accrual_props = c(0.5, 0.5))

  expect_equal(length(result), n)
  expect_true(is.numeric(result))
})

test_that("sim_accrual returns values within expected range", {
  # Test with multiple periods
  result <- sim_accrual(n = 1000,
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

test_that("sim_accrual works with single accrual period", {
  result <- sim_accrual(n = 100,
                             accrual_periods = c(10),
                             accrual_props = c(1))

  expect_equal(length(result), 100)
  expect_true(all(result >= 0))
  expect_true(all(result <= 10))
})

test_that("sim_accrual works with unequal accrual proportions", {
  set.seed(123) # For reproducibility
  result <- sim_accrual(n = 1000,
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


test_that("sim_accrual handles extreme probability distributions", {
  # Test with very skewed probability
  result <- sim_accrual(n = 1000,
                             accrual_periods = c(5, 10, 15),
                             accrual_props = c(0.01, 0.01, 0.98))

  # Most values should be in the third period (if perfect would be 980, there is wiggle room)
  expect_true(sum(result >= 10) > 950)
})

test_that("sim_accrual handles uneven inputs", {
  expect_error(sim_accrual(n = 1000,
                             accrual_periods = c(5, 10, 15),
                             accrual_props = c(0.01, 0.01)),
               "`accrual_periods` and `accrual_props` should have equal lengths")
})


# sim_pw_const_haz -------------------------------------------------------


test_that("sim_pw_const_haz returns correct number and type of values", {
  n <- 100
  hazard_periods <- c(5, 10)
  hazard_values <- c(0.1, 0.2, 0.3)

  # Run the function
  set.seed(123)
  result <- sim_pw_const_haz(n = n, hazard_periods = hazard_periods, hazard_values = hazard_values)

  # Check results
  expect_equal(length(result), n)
  expect_true(is.numeric(result))
  expect_true(all(result >= 0))
})

test_that("sim_pw_const_haz handles constant hazard case", {
  n <- 1000
  hazard_value <- 0.1

  # Run with single hazard value (constant hazard)
  set.seed(123)
  result <- sim_pw_const_haz(n = n, hazard_values = hazard_value)

  # For exponential distribution with rate 0.1, mean should be approximately 10
  expect_gt(mean(result), 9)
  expect_lt(mean(result), 11)
})

test_that("sim_pw_const_haz handles piecewise constant hazard", {
  n <- 1000
  hazard_periods <- c(5)
  hazard_values <- c(0.05, 0.8)  # Low hazard then high hazard

  # Run the function
  set.seed(123)
  result <- sim_pw_const_haz(n = n, hazard_periods = hazard_periods, hazard_values = hazard_values)

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

test_that("sim_pw_const_haz handles edge cases", {
  # Very high hazard (should result in early events)
  high_hazard_result <- sim_pw_const_haz(n = 100, hazard_values = 10)
  expect_true(mean(high_hazard_result) < 1)  # Mean should be less than 1 with hazard=10

  # Very low hazard (should result in later events)
  low_hazard_result <- sim_pw_const_haz(n = 100, hazard_values = 0.01)
  expect_true(mean(low_hazard_result) > 50)  # Mean should be greater than 50 with hazard=0.01
})

test_that("sim_pw_const_haz validates inputs", {
  # Number of hazard values should match number of periods + 1
  expect_error(
    sim_pw_const_haz(n = 100, hazard_periods = c(5, 10), hazard_values = c(0.1, 0.2)),
    "`hazard_values` should have length equal to one more than the length of `hazard_periods`"
  )

  # Negative hazard values should error
  expect_error(
    sim_pw_const_haz(n = 100, hazard_values = -0.1),
    "`hazard_values` and `hazard_periods` can only have positive values"
  )
})


# sim_weib_ph ----------------------------------------------------


test_that("sim_weib_ph returns correct format", {
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
  result <- sim_weib_ph(weibull_ph_mod = weib_model, samp_df = samp_df)

  # Check results
  expect_length(result, nrow(samp_df))
  expect_true(all(result > 0))
  expect_type(result, "double")
})

test_that("sim_weib_ph handles conditional effects correctly", {
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
  base_times <- sim_weib_ph(weibull_ph_mod = weib_model, samp_df = samp_df)

  # Generate times with positive drift (should increase survival time)
  pos_times <- sim_weib_ph(
    weibull_ph_mod = weib_model,
    samp_df = samp_df,
    cond_drift = 0.5
  )

  # Generate times with negative treatment effect (should decrease survival time)
  neg_times <- sim_weib_ph(
    weibull_ph_mod = weib_model,
    samp_df = samp_df,
    cond_trt_effect = -0.5
  )

  # Positive drift should increase survival time (higher median)
  expect_gt(median(neg_times), median(base_times))

  # Negative treatment effect should decrease survival time
  expect_lt(median(pos_times), median(base_times))
})



test_that("sim_weib_ph handles empty data correctly", {
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
  result <- sim_weib_ph(weibull_ph_mod = weib_model, samp_df = samp_df)
  expect_length(result, 0)
})


test_that("sim_weib_ph validates cond_drift input", {
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
    sim_weib_ph(weib_model, samp_df, cond_drift = "not a number"),
    "`cond_drift` must be a single number"
  )

  # Test error when cond_drift is a vector
  expect_error(
    sim_weib_ph(weib_model, samp_df, cond_drift = c(0.1, 0.2)),
    "`cond_drift` must be a single number"
  )
})

test_that("sim_weib_ph validates cond_trt_effect input", {
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
    sim_weib_ph(weib_model, samp_df, cond_trt_effect = "not a number"),
    "`cond_trt_effect` must be a single number"
  )

  # Test error when cond_trt_effect is a vector
  expect_error(
    sim_weib_ph(weib_model, samp_df, cond_trt_effect = c(0.1, 0.2)),
    "`cond_trt_effect` must be a single number"
  )
})

test_that("sim_weib_ph validates weibull_ph_mod input", {
  # Create sample data
  samp_df <- data.frame(
    cov1 = rnorm(10),
    cov2 = rbinom(10, 1, 0.5)
  )

  # Test error with non-survreg object
  expect_error(
    sim_weib_ph(weibull_ph_mod = "not a survreg object", samp_df = samp_df),
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
    sim_weib_ph(weibull_ph_mod = exp_model, samp_df = samp_df),
    "`weibull_ph_mod` must be a `survreg` object with a Wiebull distribution"
  )
})


# calc_cond_weibull -------------------------------------------------------

test_that("calc_cond_weibull produces expected results format and validates inputs", {
  # Setup test data
  set.seed(123)
  n <- 50
  df <- data.frame(
    y = rweibull(n, shape = 1.5, scale = 5),
    event = sample(0:1, n, replace = TRUE, prob = c(0.2, 0.8)),
    cov1 = rnorm(n),
    cov2 = rbinom(n, 1, 0.5),
    cov3 = rnorm(n, mean = 0.5),
    cov4 = rbinom(n, 1, 0.7)
  )

  # Create model for testing
  weib_model <- survival::survreg(survival::Surv(y, event) ~ cov1 + cov2 + cov3 + cov4, data = df, dist = "weibull")

  # Create small test population - small to make test run faster
  pop <- data.frame(
    cov1 = rnorm(1000),
    cov2 = rbinom(1000, 1, 0.5),
    cov3 = rnorm(1000, mean = 0.5),
    cov4 = rbinom(1000, 1, 0.7)
  )

  # Test with wrong model class
  expect_error(
    calc_cond_weibull(
      population = pop,
      weibull_ph_mod = "not a survreg object",
      marg_drift = c(0),
      marg_trt_eff = c(0.1),
      analysis_time = 12
    ),
    "`weibull_ph_mod` must be a survreg object"
  )

  # Test with wrong model distribution
  exp_model <- survival::survreg(survival::Surv(y, event) ~ cov1 + cov2 + cov3 + cov4, data = df, dist = "exponential")
  expect_error(
    calc_cond_weibull(
      population = pop,
      weibull_ph_mod = exp_model,
      marg_drift = c(0),
      marg_trt_eff = c(0.1),
      analysis_time = 12
    ),
    "`weibull_ph_mod` must use a weibull distribution"
  )

  # Test with wrong population format
  expect_error(
    calc_cond_weibull(
      population = as.list(pop),
      weibull_ph_mod = weib_model,
      marg_drift = c(0),
      marg_trt_eff = c(0.1),
      analysis_time = 12
    ),
    "`population` must be a tibble or dataframe"
  )

  # Test with missing covariates in population
  pop_missing <- pop |> select(-cov4)
  expect_error(
    calc_cond_weibull(
      population = pop_missing,
      weibull_ph_mod = weib_model,
      marg_drift = c(0),
      marg_trt_eff = c(0.1),
      analysis_time = 12
    ),
    "Not all covariates in `weibull_ph_mod` are in the population"
  )

  expect_error(
    calc_cond_weibull(
      population = pop,
      weibull_ph_mod = weib_model,
      marg_drift = c(0),
      marg_trt_eff = c(0.1),
      analysis_time = c(6, 12)
    ),
    "`analysis_time` must be a single number"
  )

})

test_that("calc_cond_weibull correctly calculates conditional effects", {
  # Setup test data - to avoid recomputing
  skip_on_cran() # Skip on CRAN since this is a slower test

  set.seed(123)
  n <- 100
  df <- data.frame(
    y = rweibull(n, shape = 1.5, scale = 5),
    event = sample(0:1, n, replace = TRUE, prob = c(0.2, 0.8)),
    cov1 = rnorm(n),
    cov2 = rbinom(n, 1, 0.5)
  )

  weib_model <- survival::survreg(survival::Surv(y, event) ~ cov1 + cov2, data = df, dist = "weibull")

  # Create a larger population for more accurate calculation
  set.seed(456)
  pop <- data.frame(
    cov1 = rnorm(10000),
    cov2 = rbinom(10000, 1, 0.5)
  )

  # Run function once and save results for multiple assertions
  result <- calc_cond_weibull(
    population = pop,
    weibull_ph_mod = weib_model,
    marg_drift = c(-0.1, 0, 0.1),
    marg_trt_eff = c(0, 0.1),
    analysis_time = 12
  )

  # Test structure
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 5) # 3 drift values × 2 treatment effects (-the case with 0 trt effect and neg drift)
  expect_equal(ncol(result), 6) # 6 columns in output

  # Check column names
  expected_cols <- c("marg_drift", "marg_trt_eff", "conditional_drift",
                     "true_control_surv_prob", "conditional_trt_eff",
                     "true_trt_surv_prob")
  expect_true(all(expected_cols %in% colnames(result)))

})

test_that("calc_cond_weibull survival probabilities match marginal effects", {
  # Use the same results object from previous test to avoid recomputation
  skip_on_cran() # Skip on CRAN since this is a slower test

  set.seed(123)
  n <- 100
  df <- data.frame(
    y = rweibull(n, shape = 1.5, scale = 5),
    event = sample(0:1, n, replace = TRUE, prob = c(0.2, 0.8)),
    cov1 = rnorm(n),
    cov2 = rbinom(n, 1, 0.5)
  )

  weib_model <- survival::survreg(survival::Surv(y, event) ~ cov1 + cov2, data = df, dist = "weibull")

  # Create a larger population for more accurate calculation
  set.seed(456)
  pop <- data.frame(
    cov1 = rnorm(5000), # Smaller to run faster in tests
    cov2 = rbinom(5000, 1, 0.5)
  )

  # Run function and save results
  result <- calc_cond_weibull(
    population = pop,
    weibull_ph_mod = weib_model,
    marg_drift = c(0),    # Testing just one value to speed up test
    marg_trt_eff = c(0.1),
    analysis_time = 12
  )

  # For drift = 0, treatment effect = 0.1
  trt_effect_row <- result[1, ]

  # The difference between treatment and control survival probabilities should be close to the marginal effect
  observed_effect <- trt_effect_row$true_trt_surv_prob - trt_effect_row$true_control_surv_prob
  expect_equal(observed_effect, trt_effect_row$marg_trt_eff, tolerance = 0.01)
})

