###################################################################################################
# Code to test various functions of beastt - simulation functions for binary endpoint
###################################################################################################

### Load libraries
library(beastt)
library(distributional)
library(dplyr)


###################################################################################################
# bootstrap_cov
###################################################################################################

##### Test for bootstrap_cov() without stratification (i.e., balanced covariate distributions)
test_that("bootstrap_cov returns dataset balanced with external data", {
  ## Remove all columns that aren't covariates
  ex_dat_cov <- ex_binary_df |>
    select(-subjid, -y)

  ## Obtain a very large sample
  set.seed(123)
  samp <- bootstrap_cov(external_dat = ex_dat_cov, n = 100000)

  ## Check that the means of all covariates are similar
  expect_equal(all(abs(colMeans(samp)[-1]-colMeans(ex_dat_cov)[-1]) < 0.01), TRUE)  # binary covs
  expect_equal(abs(mean(samp$cov1)-mean(ex_dat_cov$cov1)) < 0.1, TRUE)    # continuous covs
})

##### Test for bootstrap_cov() with stratification (i.e., covariate imbalance)
test_that("bootstrap_cov returns dataset with imbalanced covariate distributions", {
  ## Remove all columns that aren't covariates
  ex_dat_cov <- ex_binary_df |>
    select(-subjid, -y)

  ## Obtain a very large sample
  set.seed(123)
  samp <- bootstrap_cov(external_dat = ex_dat_cov, n = 100,
                        imbal_var = cov2, imbal_prop = .25, ref_val = 0)

  ## Check that the proportion participants with baseline cov2 is 0.25
  expect_equal(mean(samp$cov2 == 0), .25)
})

##### Test for bootstrap_cov() with multiple imbalance proportions
test_that("bootstrap_cov returns dataset with imbalanced covariate distributions", {
  ## Remove all columns that aren't covariates
  ex_dat_cov <- ex_binary_df |>
    select(-subjid, -y)

  ## Obtain a very large sample
  set.seed(123)
  samp <- bootstrap_cov(external_dat = ex_dat_cov, n = 100,
                        imbal_var = cov2, imbal_prop = c(.2, .5, .7), ref_val = 0)

  ## Check that the proportion participants with baseline cov2 is .2, .5, .7 for three samples
  expect_equal(mean(samp[[1]]$cov2 == 0), .2)
  expect_equal(mean(samp[[2]]$cov2 == 0), .5)
  expect_equal(mean(samp[[3]]$cov2 == 0), .7)
})

### Test for invalid external data
test_that("bootstrap_cov handles invalid external data", {
  expect_error(bootstrap_cov(c(1,2,3), n = 100))
  expect_error(bootstrap_cov("abc", n = 100))
})

### Test for invalid imbalance covariates
test_that("bootstrap_cov handles invalid imbalance covariates", {
  expect_error(bootstrap_cov(ex_binary_df, n = 100, imbal_var = "cov2", imbal_prop = .2))
  expect_error(bootstrap_cov(ex_binary_df, n = 100, imbal_var = cov8, imbal_prop = .2))
})

### Test for invalid imbalance proportions
test_that("bootstrap_cov handles invalid imbalance proportions", {
  expect_error(bootstrap_cov(ex_binary_df, n = 100, imbal_var = cov2, imbal_prop = -.2))
  expect_error(bootstrap_cov(ex_binary_df, n = 100, imbal_var = cov2, imbal_prop = 1.2))
})

### Test for invalid reference level
test_that("bootstrap_cov handles invalid reference level", {
  expect_error(bootstrap_cov(ex_binary_df, n = 100, imbal_var = cov2, imbal_prop = .2, ref_val = 2))
})


###################################################################################################
# calc_cond_binary
###################################################################################################

##### Test for calc_cond_binary() with a single population
test_that("calc_cond_binary returns data frame for a single population", {
  ## Fit logistic regression model using external data
  logit_mod <- glm(y ~ cov1 + cov2 + cov3 + cov4, data = ex_binary_df, family = binomial)

  ## Create a population with imbalance
  ex_dat_cov <- ex_binary_df |>
    select(-subjid, -y)
  set.seed(123)
  popn <- bootstrap_cov(external_dat = ex_dat_cov, n = 100000, imbal_var = cov2, imbal_prop = .25)

  ## Calculate conditional drift and treatment effect values
  pop_var <- calc_cond_binary(population = popn, glm = logit_mod,
                              marg_drift = c(-.1, 0, .1), marg_trt_eff = c(0, .1))

  ## Response rate of external control after standardizing covariate distributions to match the internal trial
  RR_EC <- mean( inv_logit(as.matrix(cbind(1, popn)) %*% logit_mod$coefficients) )

  ## Check that the internal and external control marginal RRs differ by the marginal drift (with tolerance)
  expect_equal(all(abs((pop_var$true_control_RR - RR_EC) - pop_var$marg_drift) < 0.00001), TRUE)
  ## Check that the internal control and trt marginal RRs differ by the marginal trt effect (with tolerance)
  expect_equal(all(abs((pop_var$true_trt_RR - pop_var$true_control_RR) - pop_var$marg_trt_eff) < 0.00001), TRUE)
})

##### Test for calc_cond_binary() with a list of populations
test_that("calc_cond_binary returns data frame for a list of populations", {
  ## Fit logistic regression model using external data
  logit_mod <- glm(y ~ cov1 + cov2 + cov3 + cov4, data = ex_binary_df, family = binomial)

  ## Create a list of populations
  ex_dat_cov <- ex_binary_df |>
    select(-subjid, -y)
  set.seed(123)
  pop1 <- bootstrap_cov(ex_dat_cov, n = 100000)
  pop23 <- bootstrap_cov(external_dat = ex_dat_cov, n = 100000, imbal_var = cov2, imbal_prop = c(.2, .5))
  pop_ls <- c(list(pop1), pop23)
  names(pop_ls) <- c("pop1", "pop2", "pop3")

  ## Calculate conditional drift and treatment effect values
  pop_var <- pop_ls |>
    map(calc_cond_binary, logit_mod, marg_drift = c(-.1, 0, .1), marg_trt_eff = c(0, .1)) |>
    bind_rows(.id = "population")

  ## Response rate of external control after standardizing covariate distributions to match the internal trial
  RR_EC <- pop_ls |>
    map(function(popn, glm){mean( inv_logit(as.matrix(cbind(1, popn)) %*% glm$coefficients) )},
        glm = logit_mod) |>
    unlist()

  ## Check that the internal and external control marginal RRs differ by the marginal drift (with tolerance)
  expect_equal(all(abs((pop_var$true_control_RR - RR_EC[pop_var$population]) - pop_var$marg_drift) < 0.00001), TRUE)
  ## Check that the internal control and trt marginal RRs differ by the marginal trt effect (with tolerance)
  expect_equal(all(abs((pop_var$true_trt_RR - pop_var$true_control_RR) - pop_var$marg_trt_eff) < 0.00001), TRUE)
})

##### Test for calc_cond_binary() with a single marginal drift and marginal trt effect value
test_that("calc_cond_binary returns data frame for single marg_drift and marg_trt_eff values", {
  ## Fit logistic regression model using external data
  logit_mod <- glm(y ~ cov1 + cov2 + cov3 + cov4, data = ex_binary_df, family = binomial)

  ## Create a population with imbalance
  ex_dat_cov <- ex_binary_df |>
    select(-subjid, -y)
  set.seed(123)
  popn <- bootstrap_cov(external_dat = ex_dat_cov, n = 100000, imbal_var = cov2, imbal_prop = .25)

  ## Calculate conditional drift and treatment effect values
  pop_var <- calc_cond_binary(population = popn, glm = logit_mod,
                              marg_drift = 0, marg_trt_eff = 0)

  ## Response rate of external control after standardizing covariate distributions to match the internal trial
  RR_EC <- mean( inv_logit(as.matrix(cbind(1, popn)) %*% logit_mod$coefficients) )

  ## Check that the internal and external control marginal RRs differ by the marginal drift (with tolerance)
  expect_equal(all(abs((pop_var$true_control_RR - RR_EC) - pop_var$marg_drift) < 0.00001), TRUE)
  ## Check that the internal control and trt marginal RRs differ by the marginal trt effect (with tolerance)
  expect_equal(all(abs((pop_var$true_trt_RR - pop_var$true_control_RR) - pop_var$marg_trt_eff) < 0.00001), TRUE)
})

##### Test for calc_cond_binary() with some scenarios filtered out if RR_EC + marg_drift + marg_trt_eff > 1
test_that("calc_cond_binary returns data frame for a single population", {
  ## Fit logistic regression model using external data
  logit_mod <- glm(y ~ cov1 + cov2 + cov3 + cov4, data = ex_binary_df, family = binomial)

  ## Create a population with imbalance
  ex_dat_cov <- ex_binary_df |>
    select(-subjid, -y)
  set.seed(123)
  popn <- bootstrap_cov(external_dat = ex_dat_cov, n = 100000, imbal_var = cov2, imbal_prop = .25)

  ## Calculate conditional drift and treatment effect values
  pop_var <- calc_cond_binary(population = popn, glm = logit_mod,
                              marg_drift = c(0, .5), marg_trt_eff = c(0, .5))

  ## Response rate of external control after standardizing covariate distributions to match the internal trial
  RR_EC <- mean( inv_logit(as.matrix(cbind(1, popn)) %*% logit_mod$coefficients) )

  ## Check that the internal and external control marginal RRs differ by the marginal drift (with tolerance)
  expect_equal(all(abs((pop_var$true_control_RR - RR_EC) - pop_var$marg_drift) < 0.00001), TRUE)
  ## Check that the internal control and trt marginal RRs differ by the marginal trt effect (with tolerance)
  expect_equal(all(abs((pop_var$true_trt_RR - pop_var$true_control_RR) - pop_var$marg_trt_eff) < 0.00001), TRUE)
})

### Test for invalid population
test_that("calc_cond_binary handles invalid population", {
  logit_mod <- glm(y ~ cov1 + cov2 + cov3 + cov4, data = ex_binary_df, family = binomial)
  expect_error(calc_cond_binary(population = "abc", glm = logit_mod, marg_drift = 0, marg_trt_eff = 0))
})

### Test for mismatched covariate names
test_that("calc_cond_binary handles mismatched covariate names", {
  logit_mod <- glm(y ~ cov1 + cov2 + cov3 + cov4, data = ex_binary_df, family = binomial)

  ## Create a population with different covariate names
  ex_dat_cov <- ex_binary_df |>
    select(-subjid, -y)
  set.seed(123)
  popn <- bootstrap_cov(external_dat = ex_dat_cov, n = 100000)
  colnames(popn) <- c("cov1", "cov2", "cov3", "abc")

  expect_error(calc_cond_binary(population = popn, glm = logit_mod, marg_drift = 0, marg_trt_eff = 0))
})



###################################################################################################
# inv_logit
###################################################################################################

##### Test for inv_logit() with a single number
test_that("inv_logit returns a single probability if given a single number", {
  val <- .3
  expect_equal(inv_logit(val), exp(val)/(1 + exp(val)))
})

##### Test for inv_logit() with a vector of numbers
test_that("inv_logit returns a single probability if given a vector of numbers", {
  val <- c(-1.5, 0, 1.5)
  expect_equal(inv_logit(val), exp(val)/(1 + exp(val)))
})

### Test for invalid input
test_that("inv_logit handles invalid input", {
  val <- "abc"
  expect_error(inv_logit(val))
})
