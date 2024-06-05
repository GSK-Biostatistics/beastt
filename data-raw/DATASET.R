
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
