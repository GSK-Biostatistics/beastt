
#' Bootstrap Binary Data
#'
#' @param external_dat External Data to bootstrap
#' @param n number of rows in the output dataset
#' @param imbal_var optional variable to to create imbalance with. If left
#'   `NULL` the distribution of the covariates will match the distribution in
#'   the external dataset. The imbalance variable must be binary.
#' @param imbal_prop optional imbalance proportion. Needed if there is an
#'   imbalance variable. This controls the proportion of subjects with the
#'   reference value.
#' @param ref_val
#'
#' @returns dataframe with the same number of columns as the reference dataframe,
#'   but n number of rows
#' @export
#'
#' @examples
#' bootstrap_binary(external_dat, n = 100000)
#' @importFrom rlang quo_is_null
#' @importFrom dplyr slice row_number
bootstrap_binary <- function(external_dat, n,
                             imbal_var = NULL, imbal_prop = NULL, ref_val = 0){
  var <- enquo(imbal_var)
  if(quo_is_null(var)){
    no_imbal_ids <- sample(1:nrow(external_dat), size = pop_size, replace = TRUE)     # sample IDs
    pop <- slice(external_dat, no_imbal_ids)
  } else if(!quo_is_null(var) && !is.null(imbal_prop)){
    uni_val <- external_dat |>
      pull({{imbal_var}}) |>
      discard(is.na) |>
      unique()
    if(length(uni_val) > 2){
      cli_abort("{.arg imbal_var} must be binary")
    }
    indexed_df <- external_dat |>
      mutate(`__row__` = row_number())
    ref_inds <- indexed_df |>
      filter({{imbal_var}} == ref_val) |>
      pull(`__row__`)
    non_ref <- indexed_df |>
      filter({{imbal_var}} != ref_val) |>
      pull(`__row__`)
    ref_size = round(imbal_prop*n)
    non_ref_size = n-ref_size
    imbal_ids <- c( sample(ref_inds,
                           size = ref_size, replace = TRUE),
                    sample(non_ref,
                           size = non_ref_size, replace = TRUE) )  # sample IDs
    pop <- slice(external_dat, imbal_ids)

  } else{
    cli_abort("{.arg imbal_var} and {.arg imbal_prop} must both have values")
  }
  pop
}


#' Calculate Conditional Drift for Binary
#'
#' In order to properly simulate values to test inverse probability weighting,
#' we need to convert the marginal drift and treatment effect to effects
#' conditioned on the predictive covariates we will use in the IPW weighting.
#'
#' @param population Large population dataset, where the columns of the dataset
#'   correspond to the covariates for the logistic regression
#' @param glm logistic regression model object
#' @param marg_drift vector of marginal drift values
#' @param marg_trt_eff vector of marginal treatment effect values
#'
#' @returns tibble of all combinations of the marginal drift and treatment
#'   effect. For each row the conditional drift and treatment effect has been
#'   calculated as well as the true control response rate and true treatment
#'   effect.
#' @export
#'
#' @importFrom dplyr left_join ungroup
#' @importFrom purrr map2_dbl
calc_cond_binary <- function(population, glm, marg_drift, marg_trt_eff){
  if(!all(class(glm) == c("glm","lm"))){
    cli_abort("{.arg beta_coefs} must be a glm object")
  }

  cov_vec <-   names(glm$coefficients) |>
    discard(\(x) x == "(Intercept)")
  beta_coefs <- glm$coefficients

  if(all(cov_vec %in% colnames(population))){
    X_IC = as.matrix(cbind(int = 1, select(population,cov_vec)))
    # Find optimal delta value
    # RR from EC arm if covariates were same as IC arm
    EC_RR_star <- mean(inv_logit(X_IC %*% beta_coefs))
    # assume True external control response rate
    external_cont_RR <- mean(glm$y)
    scenarios <- crossing(
      marg_drift,marg_trt_eff
    ) |>
      filter(marg_drift + marg_trt_eff + external_cont_RR < 1)

    delta_df <- scenarios |>
      pull(marg_drift) |>
      unique() |>
      map(function(Delta){
        optim_delta <- function(x){
          abs( Delta - ( mean(inv_logit(X_IC %*% beta_coefs + x)) - EC_RR_star ) )
        }
        delta_val <- optimize( f = optim_delta, lower = -5, upper = 5 )$minimum

        # Calculate the "true" value of theta0 (control RR)
        theta0_true_val <- mean(inv_logit(X_IC %*% beta_coefs + delta_val))
        c("marg_drift" = Delta, "conditional_drift" = delta_val, "true_control_RR" = theta0_true_val)
      }) |>
      reduce(bind_rows)

    cond_df <- scenarios |>
      left_join(delta_df, by = "marg_drift") |>
      mutate(conditional_trt_eff =
               map2_dbl(marg_trt_eff, conditional_drift,
                               function(marg_trt_eff, conditional_drift){
                                 if(marg_trt_eff == 0){
                                   beta1_val <- 0
                                 } else {
                                   optim_beta1 <- function(x){
                                     abs( marg_trt_eff - ( mean(inv_logit(X_IC %*% beta_coefs + conditional_drift + x)) -
                                                             mean(inv_logit(X_IC %*% beta_coefs + conditional_drift)) ) )
                                   }
                                   beta1_val <- optimize( f = optim_beta1, lower = -5, upper = 5 )$minimum
                                 }
                                 beta1_val
                               })) |>
      rowwise() |>
      mutate(true_trt_eff =  mean(inv_logit(X_IC %*% beta_coefs + conditional_drift + conditional_trt_eff)))|>
      ungroup()
  } else {
    cli_abort("Not all covariates in {.arg beta_coefs} are in the population")
  }
  cond_df

}

#' @noRd
inv_logit <- function(x){
  exp(x)/(1+exp(x))
}
