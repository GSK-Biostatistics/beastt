
#' Calculate Power Prior
#'
#' Function to calculate a power prior given a distribution, hyperparameters, a
#' weighted object and a response variable. The weighted object will typically
#' be a `prop_scr_obj` created by calling `create_prop_scr()`
#'
#' @param dist Name of the distribution of the external data response
#' @param hyperparameter Named vector of the hyperparameters for a given
#'   distribution. If you are using an improper prior then value can be left as
#'   `NULL`
#' @param weighted_obj A `prop_scr_obj` created by calling `calc_prop_scr()`
#' @param response_var Name of response variable
#'
#' @return `power_prior` object
#' @export
#' @importFrom dplyr summarise mutate
#' @importFrom rlang enquo
calc_power_prior <- function(dist = c("binom", "norm"), hyperparameter = NULL, weighted_obj, response_var){
  match.arg(dist)
  test_prop_scr(weighted_obj)
  response <- enquo(response_var)

  if(!is.vector(hyperparameter)){
    cli_abort("hyperpatermeters must be a vector not a list")
  }

  if(dist == 'binom'){
    calc_power_prior_beta(hyperparameter, weighted_obj, response)
  } else if (dist == "norm"){
    calc_power_prior_norm(hyperparameter, weighted_obj, response)
  }

}

#' Calculate Power Prior Beta
#'
#' @param hyperparameter Vector with the names shape1, shape2
#' @param weighted_obj A `prop_scr_obj` created by calling `create_prop_scr()`
#' @param response_var Name of response variable
#'
#' @return beta power prior object
#' @noRd
#' @importFrom rlang ensym
#' @importFrom stats dbeta rbeta
calc_power_prior_beta <- function(hyperparameter, weighted_obj, response){
  test_prop_scr(weighted_obj)

  response <- ensym(response)
  diff <- setdiff(names(hyperparameter), c("shape1", "shape2"))
  if(length(diff) >  0){
    # TODO add a check to see what it is missing to make a more explicit error
    cli_abort(c("Incorrect hyperpatermeters for a binomial distibution",
                "i" = "Bionominal becomes a beta distribution so needs `shape1` and `shape2`"))
  }

  params <- weighted_obj$external_df |>
    mutate(shape1_response = .data$`___weight___`*!!response,
           shape2_response = .data$`___weight___`*(1-!!response)) |>
    summarise(shape1 = sum(.data$shape1_response) + hyperparameter['shape1'],
              shape2 = sum(.data$shape2_response) + hyperparameter['shape2']
    )

  out <-  structure(
    list(
      model = "beta",
      parameters = c(shape1 = params$shape1,
                     shape2 = params$shape2),
      density = function(x){
        dbeta(x, params$shape1, params$shape2)
      },
      random = function(x){
        rbeta(x, params$shape1, params$shape2)
      }
    ),
    class = c("beta", "power_prior")
  )

  out
}

#' Calculate Power Prior Normal
#'
#' @param hyperparameter Either Null or a vector with the names mean, sd
#' @param weighted_obj A `prop_scr_obj` created by calling `create_prop_scr()`
#' @param response_var Name of response variable
#' @param external_control_sd Standard deviation of external control arm if
#'   assumed known. If left
#'
#' @return beta power prior object
#' @noRd
#' @importFrom rlang ensym
#' @importFrom stats dnorm rnorm
#' @importFrom dplyr pull
calc_power_prior_norm <- function(hyperparameter, weighted_obj, response, external_control_sd = NULL){
  test_prop_scr(weighted_obj)
  response <- ensym(response)

  # mean of IP-weighted power prior
  vars <- weighted_obj$external_df |>
    summarise(
      tot_ipw = sum(.data$`___weight___`),
      weight_resp = sum(.data$`___weight___`*{{response}}))
  weight_resp <- vars |> pull(.data$weight_resp)
  tot_ipw <- vars |> pull(.data$tot_ipw)

  if(is.null(external_control_sd) && !is.numeric(external_control_sd)){
    cli_abort("{.agr external_control_sd} must be a number")
  }

  if(is.null(hyperparameter)){
    # IF AN IMPROPER INITIAL PRIOR IS USED - PROPORTIONAL TO 1
    # Hyperparameters of power prior (normal distribution)
    sd_hat <- external_control_sd^2 / tot_ipw  # variance of IPW power prior
    mean_hat <- weight_resp/tot_ipw # mean of IP-weighted power prior

  } else {
    diff <- setdiff(names(hyperparameter), c("mean", "sd"))
    if(length(diff) >  0){
      # TODO add a check to see what it is missing to make a more explicit error
      cli_abort(c("Incorrect hyperpatermeters for a normal distibution",
                  "i" = "Normal distributions either needs `mean` and `sd`, or to be left `NULL`"))
    }

    sd_hat <- ( tot_ipw/external_control_sd^2 +
                  hyperparameter$sd^-2 )^-1           # variance of IP-weighted power prior
    mean_hat <- (weight_resp/external_control_sd^2 +
                   hyperparameter$mean/hyperparameter$sd^2 ) * sd_hat          # mean of IP-weighted power prior

  }

  out <-  structure(
    list(
      model = "norm",
      parameters = c(mean = mean_hat,
                     sd = sd_hat),
      density = function(x){
        dnorm(x, mean_hat, sd_hat)
      },
      random = function(x){
        rnorm(x, mean_hat, sd_hat)
      }
    ),
    class = c("normal", "power_prior")
  )

  out
}

#' @export
print.power_prior <- function(x, ...){
  param_txt <- paste0(names(x$parameters), ": ", round(x$parameters, 3))
  cli_h1("Power Prior")
  cli_bullets(c("i" = "Distribution : {.field {x$model}}",
                "*" = "{param_txt}"))

}





