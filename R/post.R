#' Calculate Posterior Normal
#'
#' @param internal_data This can either be a propensity score object or a tibble
#'   of the internal data.
#' @param response Name of response variable
#' @param prior distributional object, possibly a mixture distribution
#' @param internal_control_sd Standard deviation of internal control response data if
#'   assumed known. It can be left as `NULL` if assumed unknown
#'
#' @return distributional object
#' @export
#' @importFrom dplyr pull
#' @importFrom purrr safely map map2
#' @importFrom distributional dist_normal dist_mixture is_distribution parameters
#'
calc_post_norm<- function(
    internal_data,
    response,
    prior,
    internal_control_sd = NULL
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

  # Check resonse exsists in the data and calculate the sum
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

  # With known control SD
  if(!is.null(internal_control_sd)) {
    # Sum of responses in internal control arm
    sum_resp <- pull(data, !!response) |>
      sum()
    if(prior_fam == "normal"){
      out_dist <- calc_norm_post(prior, internal_control_sd, nIC, sum_resp)
    } else if (prior_fam == "student_t") {
      out_dist <- t_to_mixnorm(prior) |>
        calc_mixnorm_post(internal_control_sd, nIC, sum_resp)
    } else if(prior_fam == "mixture") {
      out_dist <- mix_t_to_mix(prior) |>
        calc_mixnorm_post(internal_control_sd, nIC, sum_resp)
    } else {
      cli_abort("{.agr prior} must be either normal, t, or a mixture of normals and t")
    }
  } else {
    if(prior_fam == "student_t") {
      out_dist <- t_to_mixnorm(prior) |>
        calc_t_post(nIC, pull(data, !!response))
    } else if(prior_fam == "mixture") {
      out_dist <- mix_t_to_mix(prior) |>
        calc_t_post(nIC, pull(data, !!response))
    } else {
      cli_abort("{.agr prior} must be either normal, t, or a mixture of normals and t")
    }
  }

  out_dist
}

#' Internal function to approximate T distributions to a mixture of normals
#'
#' @param x distributional object
#'
#' @return String of families
#' @noRd
#' @importFrom purrr quietly
#' @importFrom stats quantile
t_to_mixnorm <- function(x){
  # distribution is constrained to equal the mean of the t distribution
  # quantiles of t distribution
  t_qntl <- x |>
    quantile(p = seq(.001, .999, by = .001)) |>
    unlist()

  normEM_obj <- quietly(normalmixEM)(t_qntl,
                                     k = 2,                            # number of normal components
                                     mean.const = c(mean(x), mean(x)), # constraints on the means
                                     maxit = 1000)$result                     # maximum number of iterations

  norm_mix_w <- normEM_obj$lambda       # 2 x 1 vector of weights associated with each normal component
  norm_mix_mu <- normEM_obj$mu          # 2 x 1 vector of means associated with each normal component
  norm_mix_sigma <- normEM_obj$sigma    # 2 x 1 vector of SDs associated with each normal component
  norm_ls <- map2(norm_mix_mu, norm_mix_sigma, \(mu, sigma) dist_normal(mu, sigma))
  dist_mixture(!!!norm_ls, weights = c(norm_mix_w[1], 1-norm_mix_w[1]))
}


#' Takes a mixture distribution, checks if there is any t distributions and changes them to a mixture of normals
#'
#' @param x
#'
#' @return mixture distributional object
#' @importFrom distributional new_dist
#' @noRd
mix_t_to_mix <- function(x){
  dists <- get_base_families(x) |> unlist()
  if(!all(dists %in% c("normal", "student_t"))){
    cli_abort("Mixtures must be only made up of normal and t distributions")
  } else if("student_t" %in% dists){
    t_check <- TRUE
    while(t_check){
      i<-which(dists == "student_t")[1]
      new_mix <- parameters(x)$dist[[1]][[i]] |>
        t_to_mixnorm()
      old_weight <- parameters(x)$w[[1]][i]
      new_weights <- parameters(new_mix)$w[[1]]*old_weight

      #Removing the old mixture from lost of distributions and adding the new one in
      new_dist_ls <- c(parameters(x)$dist[[1]][-i],
                       parameters(new_mix)$dist[[1]]) |>
        # This converts the list back to distributional objects
        map(function(dist){
          new_dist(!!!parameters(dist), class=class(dist)[1])
        })
      # This is to make sure it adds to 1
      new_weight_vec <- c(parameters(x)$w[[1]][-i], new_weights[1],
                          1-sum(parameters(x)$w[[1]][-i], new_weights[1]))
      x <- dist_mixture(!!!new_dist_ls, weights = new_weight_vec)
      t_check <- "student_t" %in% unlist(get_base_families(x))
    }
  }
  x
}


#' @param prior mixture prior
#' @param n_ic number of internal response
#' @param sum_resp sum of the internal response
#'
#' @return mixture distributional
#' @noRd
calc_mixnorm_post <- function(prior,internal_control_sd,  n_ic, sum_resp){

  prior_means <- mix_means(prior)
  prior_sds <- mix_sigmas(prior)
  prior_ws <- parameters(prior)$w[[1]]

  # K x 1 vectors of means and SDs of each component of posterior distribution for mu_C
  post_sds <- (n_ic/internal_control_sd^2 + 1/prior_sds^2)^-.5        # vector of SDs
  post_means <- post_sds^2 * (sum_resp/internal_control_sd^2 + prior_means/prior_sds^2)   # vector of means

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
  dist_mixture(!!!post_ls, weights = post_w_norm)
}

#' @param prior normal prior
#' @param n_ic number of internal response
#' @param sum_resp sum of the internal response
#'
#' @return normal distributional
#' @noRd
calc_norm_post <- function(prior,internal_control_sd, n_ic, sum_resp){
  x <- parameters(prior)

  # K x 1 vectors of means and SDs of each component of posterior distribution for mu_C
  post_sds <- (n_ic/internal_control_sd^2 + 1/x$sigma^2)^-.5        # vector of SDs
  post_means <- post_sds^2 * (sum_resp/internal_control_sd^2 + x$mu/x$sigma^2)   # vector of means
  final_dist <- dist_normal(post_means, post_sds)
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

calc_t_post <- function(prior, nIC, response){
  # Get the integrated likelihood and convert it to a mixture of normals
  pi_tilde_IC <- calc_t(response, n= nIC)
  int_like <- t_to_mixnorm(pi_tilde_IC)
  like_mean <- mix_means(int_like)
  like_sds <- mix_sigmas(int_like)
  like_ws <- parameters(int_like)$w[[1]]

  prior_means <- mix_means(prior)
  prior_sds <- mix_sigmas(prior)
  prior_ws <- parameters(prior)$w[[1]]

  # 2*K x 1 vectors of means and SDs of each normal component of posterior distribution for mu_C
  K <- length(prior_means)
  post_sds <- c((1/like_sds[1]^2 + 1/prior_sds^2)^-.5,    # vector of SDs
                    (1/like_sds[2]^2 + 1/prior_sds^2)^-.5)
  post_means <- post_sds^2 *     # vector of means
    c((like_mean[1]/like_sds[1]^2 + prior_means/prior_sds^2),
      (like_mean[2]/like_sds[2]^2 + prior_means/prior_sds^2))

  # Create a list of normal distributions
  post_ls <- map2(post_means, post_sds, function(mu, sigma){
    dist_normal(mu, sigma)
  })

  # 2*K x 1 vector of posterior weights (unnormalized) corresponding to each component of the
  # posterior distribution for mu_C
  log_post_w_propto <-
    c(log(prior_ws) + log(like_ws[1]) - .5*log(2*pi * like_sds[1]^2 * prior_sds^2 * post_sds[1:K]^-2) -
        .5 * like_mean[1]^2/like_sds[1]^2 - .5 * prior_means^2/prior_sds^2 +
        .5 * post_means[1:K]^2 * post_sds[1:K]^-2,
      log(prior_ws) + log(like_ws[2]) - .5*log(2*pi * like_sds[2]^2 * prior_sds^2 * post_sds[(K+1):(2*K)]^-2) -
        .5 * like_mean[2]^2/like_sds[2]^2 - .5 * prior_means^2/prior_sds^2 +
        .5 * post_means[(K+1):(2*K)]^2 * post_sds[(K+1):(2*K)]^-2)

  # 2*K x 1 vector of posterior weights (normalized) corresponding to each component of the
  # posterior distribution for mu_C
  adj_log_post_w_propto <- exp(log_post_w_propto - max(log_post_w_propto))  # subtract max of log weights before exponentiating
  post_w_norm <- adj_log_post_w_propto / sum(adj_log_post_w_propto)
  sum_fixed_w <- c(1-sum(post_w_norm[-1]), post_w_norm[-1])
  dist_mixture(!!!post_ls, weights = sum_fixed_w)

}

mix_means <- function(x){
  parameters(x)$dist[[1]]|>
    map(\(dist) dist$mu) |>
    unlist()
}

mix_sigmas <- function(x){
  parameters(x)$dist[[1]]|>
    map(\(dist) dist$sigma) |>
    unlist()
}
