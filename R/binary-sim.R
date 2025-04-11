
#' Bootstrap Covariate Data
#'
#' @param external_dat Data frame of the external data from which to bootstrap
#'   covariate vectors
#' @param n Number of rows in the output dataset
#' @param imbal_var Optional variable indicating which covariate's distribution
#'   should be altered to incorporate imbalance compared to the external data.
#'   If left `NULL`, the distributions of all covariates in the output dataset
#'   will match the distributions in the external dataset. The imbalance
#'   variable must be binary.
#' @param imbal_prop Optional imbalance proportion, required if an imbalance
#'   variable is specified. This defines the proportion of individuals with the
#'   reference value of the imbalance variable in the returned dataset. This
#'   can either be a single proportion or a vector of proportions, in which case
#'   a list of datasets is returned.
#' @param ref_val Optional value corresponding to the reference level of the
#'   binary imbalance variable, if specified
#'
#' @details Covariate data can be generated for `n` individuals enrolled in
#'   the internal trial by bootstrap sampling entire covariate vectors from the
#'   external data, thus preserving the correlation between the covariates. If
#'   both `imbal_var` = `NULL` and `imbal_prop` = `NULL`, the function returns
#'   a single data frame in which the distributions of each covariate align
#'   with the covariate distributions from the external data (i.e., balanced
#'   covariate distributions across the two trials). Alternatively, covariate
#'   imbalance can be incorporated into the generated sample with respect to a
#'   binary covariate (`imbal_var`) such that a specified proportion
#'   (`imbal_prop`) of individuals in the resulting sample will have the
#'   reference level (`ref_val`) of this imbalance covariate. In this case,
#'   stratified bootstrap sampling is employed with the imbalance covariate as
#'   the stratification factor.
#'
#'   Multiple samples with varying degrees of imbalance can be generated
#'   simultaneously by defining `imbal_prop` to be a vector of values. The
#'   function then returns a list of data frames with a length equal to the
#'   number of specified imbalance proportions.
#'
#' @returns Data frame with the same number of columns as the external data frame
#'   and n number of rows (if the length of `imbal_prop` is 0 or 1); otherwise,
#'   a list of data frames with a length equal to that of `imbal_prop`
#' @export
#'
#' @examples
#' # Return one data frame with covariate distributions similar to external data
#' samp_balance <- bootstrap_cov(ex_binary_df, n = 1000)
#'
#' # Return a list of two data frames that incorporate imbalance w.r.t. covariate 2
#' samp_imbalance <- bootstrap_cov(ex_binary_df, n = 1000, imbal_var = cov2,
#'                                 imbal_prop = c(0.25, 0.5), ref_val = 0)
#' @importFrom rlang quo_is_null as_name
#' @importFrom dplyr slice row_number
bootstrap_cov <- function(external_dat, n,
                             imbal_var = NULL, imbal_prop = NULL, ref_val = 0){
  var <- enquo(imbal_var)

  if(quo_is_null(var)){

    # Draw n covariate vectors from external data frame using bootstrap sampling
    no_imbal_ids <- sample(1:nrow(external_dat), size = n, replace = TRUE)  # sample IDs
    pop <- slice(external_dat, no_imbal_ids)

  } else if(!quo_is_null(var) && !is.null(imbal_prop)){
    if(!as_name(var) %in% colnames(external_dat)){
      cli_abort("{.arg imbal_var} is not present in {.arg imbal_var}")
    }
    var_val <- external_dat |>
      pull({{imbal_var}})
    uni_val <- unique(var_val)
    if(length(uni_val) > 2){
      cli_abort("{.arg imbal_var} must be binary without missing data")
    }

    # Identify individuals in external data with reference level of imbal_var
    ref_inds <- which(var_val == ref_val)

    # Identify individuals in external data with non-reference level of imbal_var
    non_ref <-  which(var_val != ref_val)

    if(length(imbal_prop) == 1) {

      # Draw n covariate vectors from external data frame using stratified
      # bootstrap sampling wrt sampling variable
      ref_size = round(imbal_prop*n)    # number of individuals with ref level
      non_ref_size = n-ref_size         # number of individuals with non-ref level
      imbal_ids <- c( sample(ref_inds,  # sample IDs for indvdls with ref level
                             size = ref_size, replace = TRUE),
                      sample(non_ref,   # sample IDs for indvdls with non-ref level
                             size = non_ref_size, replace = TRUE) )
      pop <- slice(external_dat, imbal_ids)

    } else {

      ## Create list with multiple data frames, each with a proportion of individuals
      ## with the ref level equal to the specified values of imbal_prop (if a vector)
      # Draw n covariate vectors from external data frame using stratified
      # bootstrap sampling wrt sampling variable
      pop <- map(imbal_prop, function(ip){
        ref_size = round(ip*n)        # number of individuals with ref level
        non_ref_size = n-ref_size     # number of individuals with non-ref level
        imbal_ids <- c( sample(ref_inds,   # sample IDs for indvdls with ref level
                               size = ref_size, replace = TRUE),
                        sample(non_ref,    # sample IDs for indvdls with non-ref level
                               size = non_ref_size, replace = TRUE) )
        slice(external_dat, imbal_ids)
      })
      names(pop) = imbal_prop

    }

  } else{
    cli_abort("{.arg imbal_var} and {.arg imbal_prop} must both have values")
  }
  pop
}


#' Calculate Conditional Drift and Treatment Effect for Binary Outcome Models
#'
#' In order to properly generate binary response data for the internal trial as
#' part of a simulation study that investigates inverse probability weighting,
#' we need to translate the desired marginal drift and treatment effect to the
#' corresponding conditional drift and treatment effect that can then be added
#' into a binary outcome model (e.g., logistic regression model) used to
#' simulate response data.
#'
#' @param population A very large data frame (e.g., number of rows \eqn{\ge}
#'   100,000) where the columns correspond to the covariates defined in the
#'   logistic regression model object. This data frame should be constructed to
#'   represent the population of the internal trial according to the assumed
#'   covariate distributions (possibly imbalanced from the external data).
#' @param glm Logistic regression model object fit using the external data
#' @param marg_drift Vector of marginal drift values
#' @param marg_trt_eff Vector of marginal treatment effect values
#'
#' @details In simulation studies that investigate the properties of inverse
#'   probability weighted Bayesian dynamic borrowing, scenarios should be
#'   considered in which the underlying response rates for the internal and
#'   external control populations differ by varying amounts due to unmeasured
#'   confounding (i.e., drift, where positive values indicate a higher response
#'   rate for the internal population). While values of drift and treatment
#'   effect (i.e., risk difference) can be defined on the marginal scale for
#'   simulation studies, we must first convert these values to the conditional
#'   scale and then include these terms, along with covariates, in a logistic
#'   regression outcome model when generating response data for the internal
#'   arms. Doing so allows us to assume a relationship between the covariates
#'   and the response variable while properly accounting for drift and treatment
#'   effect.
#'
#'   To identify the conditional drift and treatment effect that correspond to
#'   specified values of marginal drift and treatment effect, we first bootstrap
#'   covariate vectors from the external data (e.g., \eqn{N \ge 100,000}) to
#'   construct a "population" that represents both the internal trial
#'   (possibly incorporating intentional covariate imbalance) and the external
#'   trial \emph{after} standardizing it to match the covariate distributions
#'   of the internal trial (allowing us to control for measured confounding
#'   from potential imbalance in the covariate distributions). Measured
#'   confounding can be incorporated into the data generation by bootstrapping
#'   a very large data frame (`population`) in which the distribution of at
#'   least one covariate is intentionally varied from that of the external data;
#'   additional \emph{unmeasured} drift can be incorporated through the
#'   translation of specified marginal values (`marg_drift`) to conditional
#'   values.
#'
#'   Let \eqn{\Delta} and \eqn{\delta} denote the marginal and conditional drift,
#'   respectively. For a specified value of \eqn{\Delta}, we can identify the
#'   corresponding \eqn{\delta} as the value that, when added as an additional
#'   term in the logistic regression model (i.e., change in the intercept) for
#'   each individual in the population, increases/decreases the
#'   population-averaged conditional probabilities of response by an amount
#'   approximately equal to \eqn{\Delta}. That is, the optimal \eqn{\delta}
#'   minimizes
#'
#'   \deqn{\left| \left( \frac{1}{N} \sum_{i=1}^N \frac{\exp \left(
#'   \boldsymbol{x}_i^\prime \boldsymbol{\beta}_{EC} + \delta \right)}{1 +
#'   \exp\left(\boldsymbol{x}_i^\prime \boldsymbol{\beta}_{EC} + \delta \right)}
#'   - \frac{1}{N} \sum_{i=1}^N \frac{\exp \left(
#'   \boldsymbol{x}_i^\prime \boldsymbol{\beta}_{EC} \right)}{1 +
#'   \exp \left(\boldsymbol{x}_i^\prime \boldsymbol{\beta}_{EC} \right)}
#'   \right) - \Delta \right|,}
#'
#'   where \eqn{\boldsymbol{\beta}_{EC}} is the vector of regression coefficients
#'   from the logistic regression model (`glm`) fit to the external control data
#'   (assumed here to be the "true" covariate effects when generating response
#'   data) and \eqn{\boldsymbol{x}_i} is a vector of covariates from the
#'   bootstrapped population of size \eqn{N}. In the formula above, the first
#'   and second terms correspond to the population-averaged conditional
#'   probabilities (i.e., the marginal response rates) of the internal control
#'   population with drift and the external control population (with covariate
#'   distributions standardized to match the internal trial), respectively.
#'
#'   If we now denote the marginal and conditional treatment effect by
#'   \eqn{\Gamma} and \eqn{\gamma}, respectively, we can use a similar process
#'   to identify the optimal \eqn{\gamma} that approximately corresponds to the
#'   specified value of \eqn{\Gamma}, which is done by minimizing the following:
#'
#'   \deqn{\left| \left( \frac{1}{N} \sum_{i=1}^N \frac{\exp \left(
#'   \boldsymbol{x}_i^\prime \boldsymbol{\beta}_{EC} + \delta + \gamma \right)}{1 +
#'   \exp\left(\boldsymbol{x}_i^\prime \boldsymbol{\beta}_{EC} + \delta + \gamma \right)}
#'   - \frac{1}{N} \sum_{i=1}^N \frac{\exp \left(
#'   \boldsymbol{x}_i^\prime \boldsymbol{\beta}_{EC} + \delta \right)}{1 +
#'   \exp \left(\boldsymbol{x}_i^\prime \boldsymbol{\beta}_{EC} + \delta \right)}
#'   \right) - \Gamma \right|,}
#'
#'   where the first term is the population-averaged conditional probabilities
#'   (i.e., the marginal response rate) of the internal treated population.
#'
#' @returns tibble of all combinations of the marginal drift and treatment
#'   effect. For each row the conditional drift and treatment effect has been
#'   calculated as well as the true control response rate and true treatment
#'   effect.
#' @export
#'
#' @examples
#' library(dplyr)
#' # Model "true" regression coefficients using the external data
#' logit_mod <- glm(y ~ cov1 + cov2 + cov3 + cov4, data = ex_binary_df, family = binomial)
#'
#' # Bootstrap internal control "population" with imbalance w.r.t. covariate 2
#' pop_int_ctrl <- bootstrap_cov(ex_binary_df, n = 100000, imbal_var = cov2,
#'                               imbal_prop = 0.25, ref_val = 0) |>
#'                               select(-subjid, -y)  # keep only covariate columns
#'
#' # Convert the marginal drift and treatment effects to conditional
#' calc_cond_binary(population = pop_int_ctrl, glm = logit_mod,
#'                  marg_drift = c(-.1, 0, .1), marg_trt_eff = c(0, .15))
#'
#' @importFrom dplyr left_join ungroup all_of
#' @importFrom purrr map2_dbl
#' @importFrom stats optimize
calc_cond_binary <- function(population, glm, marg_drift, marg_trt_eff){
  if(!inherits(glm, "glm")){
    cli_abort("{.arg glm} must be a glm object")
  }

  if(!inherits(population, "data.frame")){
    cli_abort("{.arg population} must be a tibble or dataframe. If you are using lists, check you haven't converted the dataframe into a list of vectors")
  }

  cov_vec <- names(glm$coefficients) |>
    discard(\(x) x == "(Intercept)")
  beta_coefs <- glm$coefficients

  if(all(cov_vec %in% colnames(population))){

    # Construct design matrix (with intercept) for the large sample ("population")
    # corresponding to the internal control (IC) population
    X_IC = as.matrix(cbind(int = 1, select(population, all_of(cov_vec))))

    # Marginal response rate (RR) for the external control (EC) population AFTER
    # standardizing it to match the covariate distributions of the IC population
    # (possibly imbalanced)
    EC_RR_star <- mean(inv_logit(X_IC %*% beta_coefs))

    # Remove scenarios for which the EC RR + marginal drift + marginal trt effect
    # is outside the range (0,1)
    external_cont_RR <- mean(glm$y)   # assumed "true" EC response rate
    scenarios <- crossing(
      marg_drift,marg_trt_eff
    ) |>
      filter(marg_drift + marg_trt_eff + external_cont_RR > 0 &
             marg_drift + marg_trt_eff + external_cont_RR < 1)

    # Identify the optimal conditional drift value ("delta") that corresponds to
    # the defined marginal drift value ("Delta"). If we calculate the marginal RR
    # for the IC population by averaging over each individual's conditional
    # probability of response (conditional on their covariates and delta), the
    # optimal value of the conditional drift is the delta that results in a
    # difference between the IC RR and the EC RR (after standardization) that is
    # approximately equal to the marginal drift Delta.
    delta_df <- scenarios |>
      pull(marg_drift) |>
      unique() |>
      map(function(Delta){
        if(Delta == 0){
          delta_val <- 0    # set delta = 0 if Delta = 0
        } else {
          optim_delta <- function(x){
            IC_RR <- mean(inv_logit(X_IC %*% beta_coefs + x))   # IC response rate
            abs( Delta - ( IC_RR - EC_RR_star ) )
          }
          delta_val <- optimize( f = optim_delta, lower = -5, upper = 5 )$minimum
        }

        # Calculate the "true" IC RR for each defined value of drift
        IC_RR_true_val <- mean(inv_logit(X_IC %*% beta_coefs + delta_val))
        c("marg_drift" = Delta, "conditional_drift" = delta_val, "true_control_RR" = IC_RR_true_val)
      }) |>
      bind_rows()

    # Identify the optimal conditional trt effect value ("gamma") that corresponds
    # to the defined marginal trt effect value ("Gamma"). If we calculate the
    # marginal RR for the internal treated (IT) population by averaging over each
    # individual's conditional probability of response (conditional on their
    # covariates, delta, and gamma), the optimal value of the conditional trt
    # effect is the gamma that results in a difference between the IT RR and the
    # IC RR that is approximately equal to the marginal trt effect (Gamma).

    cond_df <- scenarios |>
      left_join(delta_df, by = "marg_drift") |>
      mutate(conditional_trt_eff =
               map2_dbl(marg_trt_eff, .data$conditional_drift,
                        function(Gamma, delta){
                          if(Gamma == 0){
                            gamma_val <- 0    # set gamma = 0 if Gamma = 0
                          } else {
                            optim_gamma <- function(x){
                              IT_RR <- mean(inv_logit(X_IC %*% beta_coefs + delta + x))   # IT response rate
                              abs( Gamma - ( IT_RR - mean(inv_logit(X_IC %*% beta_coefs + delta)) ) )
                            }
                            gamma_val <- optimize( f = optim_gamma, lower = -5, upper = 5 )$minimum
                          }
                          gamma_val
                        })) |>
      rowwise() |>
      # Calculate the "true" IT RR for each defined value of drift and treatment effect
      mutate(true_trt_RR =  mean(inv_logit(X_IC %*% beta_coefs + .data$conditional_drift + .data$conditional_trt_eff)))|>
      ungroup()
  } else {
    cli_abort("Not all covariates in {.arg beta_coefs} are in the population")
  }
  cond_df

}


#' Inverse Logit Function
#'
#' @param x Real number(s) to take the inverse logit of
#'
#' @returns Vector of probabilities between 0 and 1
#' @export
#'
#' @examples
#' inv_logit(0.5)
#'
#' @details
#' This function is a short hand to \eqn{\exp(x)/(1 + \exp(x))}.
#'
inv_logit <- function(x){
  exp(x)/(1+exp(x))
}
