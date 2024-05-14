#' Calculate Posterior Normal
#'
#' @param internal_data This can either be a propensity score object or a tibble
#'   of the internal data.
#' @param response Name of response variable
#' @param internal_control_sd Assumed known SD of the internal control
#' @param prior distributional object, if you would like a mixture distribution
#'
#' @return distributional object
#' @export
#'
calc_post_norm<- function(
    internal_data,
    response,
    internal_control_sd,
    prior
){
  if(is_prop_scr(internal_data)){
    data <- internal_data$internal_df
    nIC <- internal_data$internal_df |>
      dplyr::pull(!!internal_data$id_col) |>
      unique() |>
      length()
  } else if(is.data.frame(internal_data)) {
    data <- internal_data
    nIC <- nrow(internal_data)
  } else{
    cli_abort("{.agr internal_data} either a dataset or `prop_scr` object type")
  }
  response <- enquo(response)

  check <- safely(select)(data, !!response)
  if(!is.null(check$error)){
    cli_abort("{.agr response} was not found in {.agr internal_data}")
  }

  dist_ls <- parameters(prior)$dist[[1]] |>
    map(class) |>
    map(\(x) x[1]) |>
    unlist()


  if(all(dist_ls == "dist_normal") & !is.null(internal_control_sd)){
    prior_means <- parameters(prior)$dist[[1]]|>
      map(\(x) x$mu) |>
      unlist()
    prior_sds <- parameters(prior)$dist[[1]]|>
      map(\(x) x$sigma) |>
      unlist()
    prior_ws <- parameters(prior)$w[[1]]

    # Sum of responses and standard error of response in internal control arm
    sum_resp <- pull(data, !!response) |>
      sum()
    se_IC <- internal_control_sd / sqrt(nIC)
    browser()
    # K x 1 vectors of means and SDs of each component of posterior distribution for mu_C
    post_sds <- (nIC/sd_IC^2 + 1/prior_sds^2)^-.5        # vector of SDs
    post_means <- post_sds^2 * (sum_resp/sd_IC^2 + prior_means/prior_sds^2)   # vector of means

    # K x 1 vector of posterior weights (unnormalized) corresponding to each component of the
    # posterior distribution for mu_C
    log_post_w_propto <- log(prior_ws) - .5*log(2*pi * sd_IC^2 * prior_sds^2 * post_sds^-2) -
      .5 * sum_resp^2/sd_IC^2 - .5 * prior_means^2/prior_sds^2 +
      .5 * post_means^2 * post_sds^-2

    # K x 1 vector of posterior weights (normalized) corresponding to each component of the
    # posterior distribution for mu_C
    adj_log_post_w_propto <- exp(log_post_w_propto - max(log_post_w_propto))  # subtract max of log weights before exponentiating
    post_w_norm <- adj_log_post_w_propto / sum(adj_log_post_w_propto)      # normalized posterior weights



  }


}
