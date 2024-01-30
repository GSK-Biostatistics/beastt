
#' Calculate Power Prior
#'
#' Function to calculate a power prior given a distribution, hyperparameters, a
#' weighted object and a response variable. The weighted object will typically
#' be a `prop_scr_obj` created by calling `create_prop_scr()`
#'
#' @param prior a {distributional} object that is the prior of the external data
#' @param weighted_obj A `prop_scr_obj` created by calling `calc_prop_scr()`
#' @param response Name of response variable
#' @param ... additional parameters needed depending on the type of power prior
#'   you are trying to calculate TODO: Make this a bulleted list with the
#'   additional parameters for each distirbution
#'
#' @return `power_prior` object
#' @export
#' @importFrom dplyr summarise mutate
#' @importFrom rlang enquo quo_is_missing
#' @importFrom distributional is_distribution
calc_power_prior <- function(prior, weighted_obj, response, ...){
  test_prop_scr(weighted_obj)
  response <- enquo(response)

  if(quo_is_missing(response)){
    cli_abort("{.arg response} must be a valid variable name")
  }
  # TODO add a check to see if the response variable is in the dataset

  if(is.null(prior) ||family(prior) == "normal"){
    calc_power_prior_norm(prior, weighted_obj, response, ...)
  } else if(family(prior) == 'beta'){
    calc_power_prior_beta(prior, weighted_obj, response, ...)
  }

}

#' Calculate Power Prior Beta
#'
#' @param prior a beta {distributional} object that is the prior of the external data
#' @param weighted_obj A `prop_scr_obj` created by calling `create_prop_scr()`
#' @param response_var Name of response variable
#'
#' @return beta power prior object
#'
#' @noRd
#' @importFrom rlang enquo
#' @importFrom stats dbeta rbeta
#' @importFrom distributional dist_beta
calc_power_prior_beta <- function(prior, weighted_obj, response){
  test_prop_scr(weighted_obj)

  if(!is_quosure(response)){
    response <- enquo(response)
  }

  prior_checks(prior, "beta")

  hyperparameter <- parameters(prior)

  params <- weighted_obj$external_df |>
    mutate(shape1_response = .data$`___weight___`*!!response,
           shape2_response = .data$`___weight___`*(1-!!response)) |>
    summarise(sum(.data$shape1_response) + hyperparameter['shape1'],
              sum(.data$shape2_response) + hyperparameter['shape2']
    )

  out <-  structure(
    list(
      pwr_prior = dist_beta(shape1 = params$shape1,
                            shape2 = params$shape2)
    ),
    class = c("beta", "power_prior")
  )

  out
}

#' Calculate Power Prior Normal
#'
#' @param prior either `NULL` or a normal {distributional} object that is the
#'   prior of the external data
#' @param weighted_obj A `prop_scr_obj` created by calling `create_prop_scr()`
#' @param response_var Name of response variable
#' @param external_control_sd Standard deviation of external control arm if
#'   assumed known.
#'
#' @return beta power prior object
#' @noRd
#' @importFrom rlang enquo is_quosure
#' @importFrom stats dnorm rnorm
#' @importFrom dplyr pull
calc_power_prior_norm <- function(prior, weighted_obj, response, external_control_sd = NULL){
  test_prop_scr(weighted_obj)

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
    sd_hat <- external_control_sd^2 / tot_ipw  # variance of IPW power prior
    mean_hat <- weight_resp/tot_ipw # mean of IP-weighted power prior

  } else {
    prior_checks(prior, "normal")
    hyperparameter <- parameters(prior)

    sd_hat <- ( tot_ipw/external_control_sd^2 +
                  hyperparameter$sigma^-2 )^-1           # variance of IP-weighted power prior
    mean_hat <- (weight_resp/external_control_sd^2 +
                   hyperparameter$mu/hyperparameter$sigma^2 ) * sd_hat          # mean of IP-weighted power prior

  }

  out <-  structure(
    list(
      pwr_prior = dist_normal(mu = mean_hat, sigma = sd_hat)
    ),
    class = c("normal", "power_prior")
  )

  out
}

#' @export
#' @importFrom distributional parameters
print.power_prior <- function(x, ...){
  params <- parameters(x$pwr_prior)
  param_txt <- paste0(names(params), ": ", signif(params, 3))
  cli_h1("{family(x$pwr_prior)}")
  cli_bullets(c("*" = "{param_txt}"))

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


