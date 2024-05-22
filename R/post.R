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
#' @importFrom dplyr pull
#' @importFrom purrr safely map map2
#' @importFrom distributional dist_normal dist_mixture is_distribution family parameters
#'
calc_post_norm<- function(
    internal_data,
    response,
    internal_control_sd,
    prior
){
  # Checking internal data and response variable
  if(is_prop_scr(internal_data)){
    data <- internal_data$internal_df
    nIC <- internal_data$internal_df |>
      pull(!!internal_data$id_col) |>
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

  # Checking the dirstibution and getting the family
  if(!is_distribution(prior)){
    cli_abort("{.agr prior} must be a distributional object")
  }
  prior_fam <- family(prior)

  if(prior_fam == "normal"){
    x <- parameters(prior)
    prior_means <- x$mu
    prior_sds <- x$sigma
  } else if(prior_fam == "student_t") {
    # TODO
  } else if(prior_fam == "mixture") {
    dist_ls <- get_base_families(prior)
    prior_means <- parameters(prior)$dist[[1]]|>
      map(\(x) x$mu) |>
      unlist()
    prior_sds <- parameters(prior)$dist[[1]]|>
      map(\(x) x$sigma) |>
      unlist()
    prior_ws <- parameters(prior)$w[[1]]
  } else {
    cli_abort("{.agr prior} must be either normal, t, or a mixture of normals and t")
  }


  # Sum of responses and standard error of response in internal control arm
  sum_resp <- pull(data, !!response) |>
    sum()
  # se_IC <- internal_control_sd / sqrt(nIC)

  # K x 1 vectors of means and SDs of each component of posterior distribution for mu_C
  post_sds <- (nIC/internal_control_sd^2 + 1/prior_sds^2)^-.5        # vector of SDs
  post_means <- post_sds^2 * (sum_resp/internal_control_sd^2 + prior_means/prior_sds^2)   # vector of means
  if(prior_fam == "normal"){
    final_dist <- dist_normal(post_means, post_sds)
  } else if(prior_fam %in% c("student_t", "mixture")) {

    # Create a list of normal distributions
    post_ls <- map2(post_means, post_sds, function(mu, sigma){
      dist_normal(mu, sigma)
    })


    # K x 1 vector of posterior weights (unnormalized) corresponding to each component of the
    # posterior distribution for mu_C
    log_post_w_propto <- log(prior_ws) - .5*log(2*pi * internal_control_sd^2 * prior_sds^2 * post_sds^-2) -
      .5 * sum_resp^2/internal_control_sd^2 - .5 * prior_means^2/prior_sds^2 +
      .5 * post_means^2 * post_sds^-2

    # K x 1 vector of posterior weights (normalized) corresponding to each component of the
    # posterior distribution for mu_C
    adj_log_post_w_propto <- exp(log_post_w_propto - max(log_post_w_propto))  # subtract max of log weights before exponentiating
    post_w_norm <- adj_log_post_w_propto / sum(adj_log_post_w_propto)      # normalized posterior weights
    final_dist <- dist_mixture(!!!post_ls, weights = post_w_norm)
  }
  final_dist
}

#' Internal function to approximate T distributions to a mixture of normals
#'
#' @param x distributional object
#'
#' @return String of families
#' @noRd
t_to_normal <- function(x){

}


#' Internal function to help determine the family of an object
#'
#' @param x distributional object
#'
#' @return String of families
#' @noRd
get_base_families <- function(x) {
  fam <- family(x)
  is_modified <- family(x) %in% c("mixture", "transformed", "inflated", "truncated")
  if(any(is_modified)) {
    fam[is_modified] <- lapply(parameters(x[is_modified])$dist, function(dist) lapply(dist, get_base_families))
  }
  fam
}

