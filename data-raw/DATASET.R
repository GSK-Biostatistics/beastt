
# Normal Data -------------------------------------------------------------

set.seed(123)

## Make the external dataset
n_external <- 150
sd_external_control = 5

ex_norm_df <- tibble(subjid = 1:n_external,   # participant ID
                      cov1 = round(rnorm(n_external, 50, 10)),   # covariate 1
                      cov2 = rbinom(n_external, 1, .2),          # covariate 2  (binary): ~20% level 1
                      cov3 = rbinom(n_external, 1, .6),          # covariate 3  (binary): ~60% level 1
                      cov4 = rbinom(n_external, 1, .3),          # covariate 4  (binary): ~30% level 1
                      y = rnorm(n_external, mean = 45, sd = sd_external_control))

usethis::use_data(ex_norm_df, overwrite = TRUE)
## Make the internal dataset

# Sample sizes
n_inter_control <- 60
n_inter_treat <- 60
n_internal <- n_inter_control + n_inter_treat     # total internal sample size
sd_internal_control = 3

int_norm_df <- tibble(subjid = 1:n_internal,                  # participant ID
                      cov1 = round(rnorm(n_internal, 55, 8)),          # covariate 1: mean=55, SD=8 (not mean=50, SD=10)
                      cov2 = rbinom(n_internal, 1, .3),                # covariate 2: ~30% level 1 (not 20%)
                      cov3 = rbinom(n_internal, 1, .5),                # covariate 3: ~50% level 1 (not 60%)
                      cov4 = rbinom(n_internal, 1, .3),                # covariate 4: ~30% level 1 (same)
                      trt = rep(0:1, c(n_inter_control, n_inter_treat)), # treatment (0: control, 1: active treatment)
                      y = rnorm(n_internal, mean = 50, sd = sd_internal_control))    # response

usethis::use_data(int_norm_df, overwrite = TRUE)


# Binary Data -------------------------------------------------------------

### Simulate external control (EC) data

# Sample size
nEC <- 150

# Sample baseline covariates
set.seed(123)
cov1_EC <- round(rnorm(nEC, 65, 15))  # covariate 1 (continuous): mean=65, SD=15
cov2_EC <- rbinom(nEC, 1, .3)         # covariate 2 (binary): ~30% level 1
cov3_EC <- rbinom(nEC, 1, .4)         # covariate 3 (binary): ~40% level 1
cov4_EC <- rbinom(nEC, 1, .5)         # covariate 4 (binary): ~50% level 1

# Combine baseline covariates into data frame
EC_dat <- tibble(subjid = 1:nEC,         # participant ID
                 cov1 = cov1_EC,     # covariate 1
                 cov2 = cov2_EC,     # covariate 2
                 cov3 = cov3_EC,     # covariate 3
                 cov4 = cov4_EC)     # covariate 4
X_EC <- as.matrix(cbind(int = 1, EC_dat[,2:5]))    # design matrix (with intercept)

# Underlying regression effects for simulating responses - coefficients for intercept and covariates 1-4
beta_coefs <- matrix(c(.1, .005, .3, -.2, -.4), ncol = 1)

# Simulate response probabilities and observed responses (y, 1: response, 0: no response)
prob_rspns_EC <- exp(X_EC %*% beta_coefs) / (1 + exp(X_EC %*% beta_coefs))    # response probabilities
EC_dat$y <- rbinom(nEC, 1, prob = prob_rspns_EC)                              # observed responses
ex_binary_df <- EC_dat
usethis::use_data(ex_binary_df, overwrite = TRUE)



### Simulate internal data for internal control (IC) and active treated group (TG)

# Sample sizes
nIC <- 80              # internal control
nTG <- 80              # treated group
n_int <- nIC + nTG     # total internal sample size

# Underlying treatment effect on the linear predictor scale (corresponds to approximately a 20%)
# difference in response rates between TG and IC)
TE_lp <- .75

# Sample baseline covariates (some variations from EC data in underlying covariate distributions)
set.seed(456)
cov1_int <- round(rnorm(n_int, 62, 5))   # covariate 1: mean=65, SD=2
cov2_int <- rbinom(n_int, 1, .6)         # covariate 2: ~60% level 1
cov3_int <- rbinom(n_int, 1, .4)         # covariate 3: ~40% level 1
cov4_int <- rbinom(n_int, 1, .3)         # covariate 4: ~30% level 1

# Combine baseline covariates into data frame -
int_dat <- tibble(ID = 1:n_int,                   # participant ID
                  cov1 = cov1_int,                # covariate 1
                  cov2 = cov2_int,                # covariate 2
                  cov3 = cov3_int,                # covariate 3
                  cov4 = cov4_int,                # covariate 4
                  trt = rep(0:1, c(nIC, nTG)))    # treatment (0: control, 1: active treatment)
X_int <- as.matrix(cbind(int = 1, int_dat[,2:6]))    # design matrix (with intercept)

# Simulate response probabilities and observed responses (y, 1: response, 0: no response)
prob_rspns_int <- exp(X_int %*% rbind(beta_coefs, TE_lp)) /
  (1 + exp(X_int %*% rbind(beta_coefs, TE_lp)))           # response probabilities
int_dat$y <- rbinom(n_int, 1, prob = prob_rspns_int)      # observed responses
int_binary_df <- int_dat
usethis::use_data(int_binary_df, overwrite = TRUE)
