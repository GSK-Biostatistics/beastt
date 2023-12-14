
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
#' @param weighted_obj A `prop_scr_obj` created by calling `create_prop_scr()`
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
calc_power_prior_beta <- function(hyperparameter, weighted_obj, response){
  diff <- setdiff(names(hyperparameter), c("shape1", "shape2"))
  if(length(diff) >  0){
    # TODO add a check to see what it is missing to make a more explicit error
    cli_abort(c("Incorrect hyperpatermeters for a binomial distibution",
                "i" = "Bionominal becomes a beta distribution so needs `shape1` and `shape2`"))
  }

  params <- weighted_obj$external_df |>
    mutate(shape1_response = `___weight___`*!!response,
           shape2_response = `___weight___`*(1-!!response)) |>
    summarise(shape1 = sum(shape1_response) + hyperparameter['shape1'],
              shape2 = sum(shape2_response) + hyperparameter['shape2']
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
#'
#' @return beta power prior object
#' @noRd
calc_power_prior_norm <- function(hyperparameter, weighted_obj, response){
  if(is.null(hyperparameter)){
    # mean of IP-weighted power prior
    mean_hat <- weighted_obj$external_df |>
      summarise(mean_hat = sum(`___ipw___`*!!response)/sum(`___ipw___`)) |>
      pull(mean_hat)

    # variance of IPW power prior
    # tau2_hat_EC <- sd_EC^2 / sum(weighted_obj$external_df$`___ipw___`)

  } else {
    diff <- setdiff(names(hyperparameter), c("mean", "sd"))
    if(length(diff) >  0){
      # TODO add a check to see what it is missing to make a more explicit error
      cli_abort(c("Incorrect hyperpatermeters for a normal distibution",
                  "i" = "Normal distributions either needs `mean` and `sd`, or to be left `NULL`"))
    }
  }

  out <-  structure(
    list(
      model = "norm",
      parameters = c(mean = mean_hat,
                     sd = params$shape2),
      density = function(x){
        dnorm(x, mean_hat, params$shape2)
      },
      random = function(x){
        rnorm(x, mean_hat, params$shape2)
      }
    ),
    class = c("normal", "power_prior")
  )

  out
}


