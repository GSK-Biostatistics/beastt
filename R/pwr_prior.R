
#' Calculate Power Prior Beta
#'
#' @param prior a beta {distributional} object that is the prior of the external data
#' @param weighted_obj A `prop_scr_obj` created by calling `create_prop_scr()`
#' @param response Name of response variable
#'
#' @return beta power prior object
#' @export
#'
#' @importFrom rlang enquo
#' @importFrom stats dbeta rbeta
#' @importFrom distributional dist_beta
#' @importFrom dplyr summarise
#' @family power prior
calc_power_prior_beta <- function(prior, weighted_obj, response){
  test_prop_scr(weighted_obj)
  response <- enquo(response)

  prior_checks(prior, "beta")

  hyperparameter <- parameters(prior)

  params <- weighted_obj$external_df |>
    mutate(shape1_response = .data$`___weight___`*!!response,
           shape2_response = .data$`___weight___`*(1-!!response)) |>
    summarise(sum(.data$shape1_response) + hyperparameter['shape1'],
              sum(.data$shape2_response) + hyperparameter['shape2']
    )

  dist_beta(shape1 = params$shape1,
            shape2 = params$shape2)
}

#' Calculate Power Prior Normal
#'
#' @param prior either `NULL` or a normal {distributional} object that is the
#'   prior of the external data
#' @param weighted_obj A `prop_scr_obj` created by calling `create_prop_scr()`
#' @param response Name of response variable
#' @param external_control_sd Standard deviation of external control arm if
#'   assumed known.
#'
#' @return beta power prior object
#' @export
#' @importFrom rlang enquo is_quosure
#' @importFrom dplyr pull
#' @importFrom distributional parameters dist_normal
#' @family power prior
calc_power_prior_norm <- function(prior, weighted_obj, response, external_control_sd = NULL){
  test_prop_scr(weighted_obj)
  response <- enquo(response)

  # mean of IP-weighted power prior
  vars <- weighted_obj$external_df |>
    summarise(
      tot_ipw = sum(.data$`___weight___`),
      weight_resp = sum(.data$`___weight___`*!!response))
  weight_resp <- vars |> pull(.data$weight_resp)
  tot_ipw <- vars |> pull(.data$tot_ipw)

  if(is.null(external_control_sd) && !is.numeric(external_control_sd)){
    cli_abort("{.agr external_control_sd} must be a number")
  }

  if(is.null(prior)){
    # IF AN IMPROPER INITIAL PRIOR IS USED - PROPORTIONAL TO 1
    # Hyperparameters of power prior (normal distribution)
    sd2_hat <- external_control_sd^2 / tot_ipw  # variance of IPW power prior
    mean_hat <- weight_resp/tot_ipw # mean of IP-weighted power prior

  } else {
    prior_checks(prior, "normal")
    hyperparameter <- parameters(prior)

    sd2_hat <- ( tot_ipw/external_control_sd^2 +
                  hyperparameter$sigma^-2 )^-1           # variance of IP-weighted power prior
    mean_hat <- (weight_resp/external_control_sd^2 +
                   hyperparameter$mu/hyperparameter$sigma^2 ) * sd2_hat          # mean of IP-weighted power prior

  }

  dist_normal(mu = mean_hat, sigma = sqrt(sd2_hat))
}




#' Prior checks
#'
#' @param prior distributional object
#' @param family string of the type of dist the prior should be
#'
#' @return errors if prior object is unsuitable
#'
#' @keywords internal
#' @importFrom stats family
#' @noRd
prior_checks <- function(prior, family){
  if(!is_distribution(prior)){
    cli_abort("Needs to be a normal or beta prior")
  } else if(length(prior) > 1){
    cli_abort("Needs to be a single prior, not a vector of priors")
  } else if(family(prior) != family){
    cli_abort("Prior need to be a {family} distribution")
  }
}

