
#' Calculate Power Prior Beta
#'
#' @description Calculate a (potentially inverse probability weighted) beta
#'   power prior for the control response rate using external control data.
#'
#' @param external_data This can either be a `prop_scr_obj` created by calling
#'   `create_prop_scr()` or a tibble of the external data. If it is just a
#'   tibble the weights will be assumed to be 1.
#' @param response Name of response variable
#' @param prior A beta distributional object that is the initial prior for the
#'   control response rate before the external control data are observed
#'
#' @details Weighted participant-level response data from an external study are
#'   incorporated into an inverse probability weighted (IPW) power prior for the
#'   control response rate \eqn{\theta_C}. When borrowing information from an
#'   external control arm of size \eqn{N_{EC}}, the components of the IPW power
#'   prior for \eqn{\theta_C} are defined as follows:
#'   \describe{
#'   \item{Initial prior:}{\eqn{\theta_C \sim \mbox{Beta}(\nu_0, \phi_0)}}
#'   \item{IPW likelihood of the external response data \eqn{\boldsymbol{y}_E} with
#'     weights \eqn{\hat{\boldsymbol{a}}_0}:}{\eqn{\mathcal{L}_E(\theta_C \mid
#'     \boldsymbol{y}_E, \hat{\boldsymbol{a}}_0) \propto \exp \left( \sum_{i=1}^{N_{EC}}
#'     \hat{a}_{0i} \left[ y_i \log(\theta_C) + (1 - y_i) \log(1 - \theta_C) \right] \right)}}
#'   \item{IPW power prior:}{\eqn{\theta_C \mid \boldsymbol{y}_E, \hat{\boldsymbol{a}}_0
#'     \sim \mbox{Beta} \left( \sum_{i=1}^{N_{EC}} \hat{a}_{0i} y_i + \nu_0,
#'     \sum_{i=1}^{N_{EC}} \hat{a}_{0i} (1 - y_i) + \phi_0 \right)}}
#'   }
#'
#'   Defining the weights \eqn{\hat{\boldsymbol{a}}_0} to equal 1 results in a
#'   conventional beta power prior.
#'
#' @return Beta power prior object
#' @export
#'
#' @importFrom rlang enquo
#' @importFrom stats dbeta rbeta
#' @importFrom distributional dist_beta
#' @importFrom dplyr summarise
#' @family power prior
#' @examples
#' library(distributional)
#' library(dplyr)
#' # This function can be used directly on the data
#' calc_power_prior_beta(external_data = ex_binary_df,
#'   response = y,
#'   prior = dist_beta(0.5, 0.5))
#'
#' # Or this function can be used with a propensity score object
#' ps_obj <- calc_prop_scr(internal_df = filter(int_binary_df, trt == 0),
#'   external_df = ex_binary_df,
#'   id_col = subjid,
#'   model = ~ cov1 + cov2 + cov3 + cov4)
#'
#' calc_power_prior_beta(ps_obj,
#'   response = y,
#'   prior = dist_beta(0.5, 0.5))
calc_power_prior_beta <- function(external_data, response, prior){
  if(is_prop_scr(external_data)){
    data <- external_data$external_df
    weights <- data$`___weight___`
  } else if(is.data.frame(external_data)){
    data <- external_data
    weights <- 1
  } else {
    cli_abort("{.agr external_data} either a dataset or `prop_scr` object type")
  }

  # Check the response
  response <- enquo(response)
  check <- safely(select)(data, !!response)
  if(!is.null(check$error)){
    cli_abort("{.agr response} was not found in {.agr external_data}")
  }

  prior_checks(prior, "beta")

  hyperparameter <- parameters(prior)

  params <- data |>
    mutate(shape1_response = weights*!!response,
           shape2_response = weights*(1-!!response)) |>
    summarise(sum(.data$shape1_response) + hyperparameter['shape1'],
              sum(.data$shape2_response) + hyperparameter['shape2']
    )

  dist_beta(shape1 = params$shape1,
            shape2 = params$shape2)
}

#' Calculate Power Prior Normal
#'
#' @description Calculate a (potentially inverse probability weighted) normal
#'   power prior using external data.
#'
#' @param external_data This can either be a `prop_scr_obj` created by calling
#'   `create_prop_scr()` or a tibble of the external data. If it is just a
#'   tibble the weights will be assumed to be 1. Only the external data for the
#'   arm(s) of interest should be included in this object (e.g., external
#'   control data if creating a power prior for the control mean)
#' @param response Name of response variable
#' @param prior Either `NULL` or a normal distributional object that is the
#'   initial prior for the parameter of interest (e.g., control mean) before the
#'   external data are observed
#' @param external_sd Standard deviation of external response data if assumed
#'   known. It can be left as `NULL` if assumed unknown
#'
#' @details Weighted participant-level response data from an external study are
#'   incorporated into an inverse probability weighted (IPW) power prior for the
#'   parameter of interest \eqn{\theta} (e.g., the control mean if borrowing
#'   from an external control arm). When borrowing information from an external
#'   dataset of size \eqn{N_{E}}, the IPW likelihood of the external response
#'   data \eqn{y_E} with weights \eqn{\hat{\boldsymbol{a}}_0} is defined as
#'
#'   \deqn{\mathcal{L}_E(\theta \mid \boldsymbol{y}_E, \hat{\boldsymbol{a}}_0,
#'   \sigma_{E}^2) \propto \exp \left( -\frac{1}{2 \sigma_{E}^2}
#'   \sum_{i=1}^{N_{E}} \hat{a}_{0i} (y_i - \theta)^2 \right).}
#'
#'   The `prior` argument should be either a distributional object with a family
#'   type of `normal` or `NULL`, corresponding to the use of a normal initial
#'   prior or an improper uniform initial prior (i.e., \eqn{\pi(\theta) \propto
#'   1}), respectively.
#'
#'   The `external_sd` argument can be a positive value if the external standard
#'   deviation is assumed known or left as `NULL` otherwise. If `external_sd =
#'   NULL`, then `prior` must be `NULL` to indicate the use of an improper
#'   uniform initial prior for \eqn{\theta}, and an improper prior is defined
#'   for the unknown external standard deviation such that \eqn{\pi(\sigma_E^2)
#'   \propto (\sigma_E^2)^{-1}}. The details of the IPW power prior for each
#'   case are as follows:
#'   \describe{
#'   \item{`external_sd = positive value` (\eqn{\sigma_E^2} known):}{With
#'   either a proper normal or an improper uniform initial prior, the IPW
#'   weighted power prior for \eqn{\theta} is a normal distribution.}
#'   \item{`external_sd = NULL` (\eqn{\sigma_E^2} unknown):}{With improper
#'   priors for both \eqn{\theta} and \eqn{\sigma_E^2}, the marginal IPW weighted
#'   power prior for \eqn{\theta} after integrating over \eqn{\sigma_E^2} is
#'   a non-standardized \eqn{t} distribution.}
#'   }
#'
#'   Defining the weights \eqn{\hat{\boldsymbol{a}}_0} to equal 1 results in a
#'   conventional normal (or \eqn{t}) power prior if the external standard
#'   deviation is known (unknown).
#'
#' @return Normal power prior object
#' @export
#' @importFrom rlang enquo is_quosure
#' @importFrom dplyr pull
#' @importFrom distributional parameters dist_normal
#' @family power prior
#' @examples
#' library(distributional)
#' library(dplyr)
#' # This function can be used directly on the data
#' calc_power_prior_norm(ex_norm_df,
#'   response = y,
#'   prior = dist_normal(50, 10),
#'   external_sd = 0.15)
#'
#' # Or this function can be used with a propensity score object
#' ps_obj <- calc_prop_scr(internal_df = filter(int_norm_df, trt == 0),
#'                        external_df = ex_norm_df,
#'                        id_col = subjid,
#'                        model = ~ cov1 + cov2 + cov3 + cov4)
#' calc_power_prior_norm(ps_obj,
#'                      response = y,
#'                      prior = dist_normal(50, 10),
#'                      external_sd = 0.15)
#'
calc_power_prior_norm <- function(external_data, response, prior = NULL, external_sd = NULL){
  if(is_prop_scr(external_data)){
    data <- external_data$external_df
    weights <- data$`___weight___`
  } else if(is.data.frame(external_data)){
    data <- external_data
    weights <- 1
  } else {
    cli_abort("{.agr external_data} either a dataset or `prop_scr` object type")
  }

  # Check the response
  response <- enquo(response)
  check <- safely(select)(data, !!response)
  if(!is.null(check$error)){
    cli_abort("{.agr response} was not found in {.agr external_data}")
  }

  # mean of IP-weighted power prior
  vars <- data |>
    summarise(
      tot_ipw = sum(weights),
      weight_resp = sum(weights*!!response))
  weight_resp <- vars |> pull(.data$weight_resp)
  tot_ipw <- vars |> pull(.data$tot_ipw)

  ec_sd_test <- !is.null(external_sd) && is.numeric(external_sd) && external_sd > 0

  if(is.null(prior)){
    if(ec_sd_test){
      # IF AN IMPROPER INITIAL PRIOR IS USED - PROPORTIONAL TO 1
      # Hyperparameters of power prior (normal distribution)
      sd2_hat <- external_sd^2 / tot_ipw  # variance of IPW power prior
      mean_hat <- weight_resp/tot_ipw # mean of IP-weighted power prior
      out_dist <- dist_normal(mu = mean_hat, sigma = sqrt(sd2_hat))
    } else {
      out_dist <- calc_t(pull(data, !!response),
                         n= nrow(data),
                         W = weights)
    }
  } else if(ec_sd_test) {
    prior_checks(prior, "normal")
    hyperparameter <- parameters(prior)

    sd2_hat <- ( tot_ipw/external_sd^2 +
                  hyperparameter$sigma^-2 )^-1           # variance of IP-weighted power prior
    mean_hat <- (weight_resp/external_sd^2 +
                   hyperparameter$mu/hyperparameter$sigma^2 ) * sd2_hat          # mean of IP-weighted power prior
    out_dist <- dist_normal(mu = mean_hat, sigma = sqrt(sd2_hat))
  } else {
    cli_abort("{.agr external_sd} must be a positive number if a prior is being supplied")
  }
  out_dist

}

#' Calculate a t distribution power prior
#'
#' @param Y Response
#' @param n Number of participants
#' @param W Optional vector of weights
#'
#' @return t distributional object
#' @importFrom distributional dist_student_t
#' @importFrom mixtools normalmixEM
#' @noRd
calc_t <- function(Y, n, W =NULL){
  # Degrees of freedom
  df <- n - 1
  # Location hyperparameter (easier to write in matrix form than scalar form)
  Y_vec <- as.matrix(Y)          # response vector
  Z <- matrix(1, nrow = n, ncol = 1)     # n x 1 vector of 1s
  if(is.null(W)){
    A <- diag(length(Y))
  } else {
    A <- diag(W)      # diagonal matrix with weights along diagonals
  }

  theta <- as.numeric(                              # location hyperparameter
    solve(t(Z) %*% A %*% Z) %*%
      t(Z) %*% A %*% Y_vec
  )

  # Dispersion hyperparameter (easier to write in matrix form than scalar form)
  V <- Z %*% solve(t(Z) %*% A %*% Z) %*%     # n x n matrix
    t(Z) %*% A
  tau2 <- as.numeric(                                    # dispersion hyperparameter
    df^-1 * solve(t(Z) %*% A %*% Z) %*%
    (t(Y_vec) %*% (A - t(V) %*% A %*% V) %*% Y_vec)
  )

  # Power prior object
  dist_student_t(df = df,              # degrees of freedom
                           mu = theta,           # location hyperparameter
                           sigma = sqrt(tau2))   # scale hyperparameter (sqrt of dispersion)

}


#' Prior checks
#'
#' @param prior Distributional object
#' @param family String of the type of dist the prior should be
#'
#' @return Errors if prior object is unsuitable
#'
#' @keywords internal
#' @importFrom stats family
#' @noRd
prior_checks <- function(prior, family){
  if(!is_distribution(prior)){
    cli_abort("Needs to be a distributional prior")
  } else if(length(prior) > 1){
    cli_abort("Needs to be a single prior, not a vector of priors")
  } else if(family(prior) != family){
    cli_abort("Prior need to be a {family} distribution")
  }
}

