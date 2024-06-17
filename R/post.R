#' Calculate Posterior Normal
#'
#' @description Calculate a posterior distribution that is normal (or a mixture
#'   of normal components). Only the relevant treatment arms from the internal
#'   dataset should be read in (e.g., only the control arm if constructing a
#'   posterior distribution for the control mean).
#'
#' @param internal_data This can either be a propensity score object or a tibble
#'   of the internal data.
#' @param response Name of response variable
#' @param prior A distributional object corresponding to a normal distribution,
#'   a t distribution, or a mixture distribution of normal and/or t components
#' @param internal_sd Standard deviation of internal response data if
#'   assumed known. It can be left as `NULL` if assumed unknown
#'
#' @details For a given arm of an internal trial (e.g., the control arm or an
#'   active treatment arm) of size \eqn{N_I}, suppose the response data are normally
#'   distributed such that \eqn{y_i \sim N(\theta, \sigma_I^2)}, \eqn{i=1,\ldots,N_I}.
#'   If \eqn{\sigma_I^2} is assumed known, the posterior distribution for \eqn{\theta}
#'   is written as
#'
#'   \deqn{\pi( \theta \mid \boldsymbol{y}, \sigma_{I}^2 ) \propto \mathcal{L}(\theta \mid \boldsymbol{y}, \sigma_{I}^2) \; \pi(\theta),}
#'
#'   where \eqn{\mathcal{L}(\theta \mid \boldsymbol{y}, \sigma_{I}^2)} is the
#'   likelihood of the response data from the internal arm and \eqn{\pi(\theta)}
#'   is a prior distribution on \eqn{\theta} (either a normal distribution, a
#'   \eqn{t} distribution, or a mixture distribution with an arbitrary number of
#'   normal and/or \eqn{t} components). Any \eqn{t} components of the prior for
#'   \eqn{\theta} are approximated with a mixture of two normal distributions.
#'
#'   If \eqn{\sigma_I^2} is unknown, the marginal posterior distribution for
#'   \eqn{\theta} is instead written as
#'
#'   \deqn{\pi( \theta \mid \boldsymbol{y} ) \propto \left\{ \int_0^\infty \mathcal{L}(\theta, \sigma_{I}^2 \mid \boldsymbol{y}) \; \pi(\sigma_{I}^2) \; d\sigma_{I}^2 \right\} \times \pi(\theta).}
#'
#'   In this case, the prior for \eqn{\sigma_I^2} is chosen to be
#'   \eqn{\pi(\sigma_{I}^2) = (\sigma_I^2)^{-1}} such that
#'   \eqn{\left\{ \int_0^\infty \mathcal{L}(\theta, \sigma_{I}^2 \mid \boldsymbol{y}) \; \pi(\sigma_{I}^2) \; d\sigma_{I}^2 \right\}}
#'   becomes a non-standardized \eqn{t} distribution. This integrated likelihood
#'   is then approximated with a mixture of two normal distributions.
#'
#'   If `internal_sd` is supplied a positive value and `prior` corresponds to a
#'   single normal distribution, then the posterior distribution for \eqn{\theta}
#'   is a normal distribution. If `internal_sd = NULL` or if other types of prior
#'   distributions are specified (e.g., mixture or t distribution), then the
#'   posterior distribution is a mixture of normal distributions.
#'
#' @return distributional object
#' @export
#' @importFrom dplyr pull
#' @importFrom purrr safely map map2
#' @importFrom distributional dist_normal dist_mixture is_distribution parameters
#' @examples
#' library(distributional)
#' library(dplyr)
#' post_treated <- calc_post_norm(internal_data = filter(int_norm_df, trt == 1),
#'                                response = y,
#'                                prior = dist_normal(50, 10),
#'                                internal_sd = 0.25)
#'
calc_post_norm<- function(
    internal_data,
    response,
    prior,
    internal_sd = NULL
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

  # Check response exists in the data and calculate the sum
  response <- enquo(response)
  check <- safely(select)(data, !!response)
  if(!is.null(check$error)){
    cli_abort("{.agr response} was not found in {.agr internal_data}")
  }

  # Checking the distribution and getting the family
  if(!is_distribution(prior)){
    cli_abort("{.agr prior} must be a distributional object")
  }
  prior_fam <- family(prior)

  # With known internal SD
  if(!is.null(internal_sd)) {
    # Sum of responses in internal arm
    sum_resp <- pull(data, !!response) |>
      sum()
    if(prior_fam == "normal"){
      out_dist <- calc_norm_post(prior, internal_sd, nIC, sum_resp)
    } else if (prior_fam == "student_t") {
      out_dist <- t_to_mixnorm(prior) |>
        calc_mixnorm_post(internal_sd, nIC, sum_resp)
    } else if(prior_fam == "mixture") {
      out_dist <- mix_t_to_mix(prior) |>
        calc_mixnorm_post(internal_sd, nIC, sum_resp)
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



#' Calculate Posterior Beta
#'
#' @description Calculate a posterior distribution that is beta (or a mixture
#'   of beta components). Only the relevant treatment arms from the internal
#'   dataset should be read in (e.g., only the control arm if constructing a
#'   posterior distribution for the control response rate).
#'
#' @param internal_data This can either be a propensity score object or a tibble
#'   of the internal data.
#' @param response Name of response variable
#' @param prior A distributional object corresponding to a beta distribution
#'   or a mixture distribution of beta components
#'
#' @details For a given arm of an internal trial (e.g., the control arm or an
#' active treatment arm) of size \eqn{N_I}, suppose the response data are binary
#' such that \eqn{y_i \sim \mbox{Bernoulli}(\theta)}, \eqn{i=1,\ldots,N_I}. The
#' posterior distribution for \eqn{\theta} is written as
#'
#' \deqn{\pi( \theta \mid \boldsymbol{y} ) \propto \mathcal{L}(\theta \mid \boldsymbol{y}) \; \pi(\theta),}
#'
#' where \eqn{\mathcal{L}(\theta \mid \boldsymbol{y})} is the likelihood of the
#' response data from the internal arm and \eqn{\pi(\theta)} is a prior
#' distribution on \eqn{\theta} (either a beta distribution or a mixture
#' distribution with an arbitrary number of beta components). The posterior
#' distribution for \eqn{\theta} is either a beta distribution or a mixture of
#' beta components depending on whether the prior is a single beta
#' distribution or a mixture distribution.
#'
#' @return distributional object
#' @export
#' @importFrom dplyr pull
#' @importFrom purrr safely map map2
#' @importFrom distributional dist_beta dist_mixture is_distribution parameters
#' @examples
#' library(dplyr)
#' library(distributional)
#' calc_post_beta(internal_data = filter(int_binary_df, trt == 1),
#'                               response = y,
#'                               prior = dist_beta(0.5, 0.5))
calc_post_beta<- function(internal_data, response, prior){
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

  # Check response exists in the data and calculate the sum
  response <- enquo(response)
  check <- safely(select)(data, !!response)
  if(!is.null(check$error)){
    cli_abort("{.agr response} was not found in {.agr internal_data}")
  }

  # Checking the distribution and getting the family
  if(!is_distribution(prior)){
    cli_abort("{.agr prior} must be a distributional object")
  }
  prior_fam <- family(prior)
  all_fam <- get_base_families(prior) |> unlist()
  if(all(all_fam == "beta")){
    # Sum of responses in internal control arm
    sum_resp <- pull(data, !!response) |>
      sum()
    if(prior_fam == "beta"){
      shape1_new <- parameters(prior)$shape1 + sum_resp
      shape2_new <- parameters(prior)$shape2 + nIC - sum_resp
      out_dist <- dist_beta(shape1_new, shape2_new)
    } else if(prior_fam == "mixture"){
      shape1_vec <- parameters(prior)$dist[[1]]|>
        map(\(dist) dist$shape1) |>
        unlist()
      shape1_new <- shape1_vec + sum_resp
      shape2_vec <- parameters(prior)$dist[[1]]|>
        map(\(dist) dist$shape2) |>
        unlist()
      shape2_new <- shape2_vec + nIC - sum_resp

      weights <- parameters(prior)$w[[1]]
      log_post_w_propto <- lbeta(shape1_new, shape2_new) -lbeta(shape1_vec, shape2_vec) + log(weights)
      adj_log_post_w_propto <- exp(log_post_w_propto - max(log_post_w_propto))  # subtract max of log weights before exponentiating
      post_w_norm <- adj_log_post_w_propto / sum(adj_log_post_w_propto)      # normalized posterior weights

      dist_ls <- map2(shape1_new, shape2_new, dist_beta)
      out_dist <- dist_mixture(!!!dist_ls, weights = correct_weights(post_w_norm))
    }
  } else {
    cli_abort("{.agr prior} must be either beta or a mixture of betas")
  }

  out_dist
}



#' Internal function to approximate t distributions to a mixture of normals
#'
#' @param x A distributional object
#'
#' @return String of families
#' @noRd
#' @importFrom purrr quietly
#' @importFrom stats quantile
#' @importFrom distributional variance
t_to_mixnorm <- function(x){
  # distribution is constrained to equal the mean of the t distribution
  # quantiles of t distribution
  t_qntl <- x |>
    quantile(p = seq(.001, .999, by = .001)) |>
    unlist()

  normEM_obj <- quietly(normalmixEM)(t_qntl,
                                     k = 2,                                # number of normal components
                                     mean.const = c(mean(x), mean(x)),     # constraints on the means
                                     lambda = c(.5, .5),                   # starting values for weights (can be equal)
                                     mu = c(mean(x), mean(x)),             # starting values for mean (must be equal)
                                     sigma = c(sqrt(variance(x)) + .001,   # starting values for SD (should NOT be equal)
                                               sqrt(variance(x)) + .002),
                                     maxit = 1000)$result                  # maximum number of iterations

  norm_mix_w <- normEM_obj$lambda       # 2 x 1 vector of weights associated with each normal component
  norm_mix_mu <- normEM_obj$mu          # 2 x 1 vector of means associated with each normal component
  norm_mix_sigma <- normEM_obj$sigma    # 2 x 1 vector of SDs associated with each normal component
  norm_ls <- map2(norm_mix_mu, norm_mix_sigma, \(mu, sigma) dist_normal(mu, sigma))
  dist_mixture(!!!norm_ls, weights = correct_weights(norm_mix_w))
}


#' Takes a mixture distribution, checks if there is any t distributions and changes them to a mixture of normals
#'
#' @param x A mixture distributional object
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
      new_weight_vec <- c(parameters(x)$w[[1]][-i], new_weights) |>
        correct_weights()
      x <- dist_mixture(!!!new_dist_ls, weights = new_weight_vec)
      t_check <- "student_t" %in% unlist(get_base_families(x))
    }
  }
  x
}


#' @param prior Mixture prior
#' @param internal_sd Standard deviation of internal response data
#' @param n_ic Number of internal response
#' @param sum_resp Sum of the internal response
#'
#' @return Mixture distributional
#' @noRd
calc_mixnorm_post <- function(prior, internal_sd, n_ic, sum_resp){

  prior_means <- mix_means(prior)
  prior_sds <- mix_sigmas(prior)
  prior_ws <- parameters(prior)$w[[1]]

  # K x 1 vectors of means and SDs of each component of posterior distribution for theta
  post_sds <- (n_ic/internal_sd^2 + 1/prior_sds^2)^-.5        # vector of SDs
  post_means <- post_sds^2 * (sum_resp/internal_sd^2 + prior_means/prior_sds^2)   # vector of means

  # Create a list of normal distributions
  post_ls <- map2(post_means, post_sds, function(mu, sigma){
    dist_normal(mu, sigma)
  })

  # K x 1 vector of posterior weights (unnormalized) corresponding to each component of the
  # posterior distribution for theta
  log_post_w_propto <- log(prior_ws) - .5*log(2*pi * internal_sd^2 * prior_sds^2 * post_sds^-2) -
    .5 * sum_resp^2/internal_sd^2 - .5 * prior_means^2/prior_sds^2 +
    .5 * post_means^2 * post_sds^-2

  # K x 1 vector of posterior weights (normalized) corresponding to each component of the
  # posterior distribution for theta
  adj_log_post_w_propto <- exp(log_post_w_propto - max(log_post_w_propto))  # subtract max of log weights before exponentiating
  post_w_norm <- adj_log_post_w_propto / sum(adj_log_post_w_propto)      # normalized posterior weights
  dist_mixture(!!!post_ls, weights = correct_weights(post_w_norm))
}

#' @param prior Normal prior
#' @param internal_sd Standard deviation of internal response data
#' @param n_ic Number of internal response
#' @param sum_resp Sum of the internal response
#'
#' @return Normal distributional
#' @noRd
calc_norm_post <- function(prior, internal_sd, n_ic, sum_resp){
  x <- parameters(prior)

  # K x 1 vectors of means and SDs of each component of posterior distribution for theta
  post_sds <- (n_ic/internal_sd^2 + 1/x$sigma^2)^-.5        # vector of SDs
  post_means <- post_sds^2 * (sum_resp/internal_sd^2 + x$mu/x$sigma^2)   # vector of means
  final_dist <- dist_normal(post_means, post_sds)
}


#' Internal function to help determine the family of an object
#'
#' @param x Distributional object
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

  # 2*K x 1 vectors of means and SDs of each normal component of posterior distribution for theta
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
  # posterior distribution for theta
  log_post_w_propto <-
    c(log(prior_ws) + log(like_ws[1]) - .5*log(2*pi * like_sds[1]^2 * prior_sds^2 * post_sds[1:K]^-2) -
        .5 * like_mean[1]^2/like_sds[1]^2 - .5 * prior_means^2/prior_sds^2 +
        .5 * post_means[1:K]^2 * post_sds[1:K]^-2,
      log(prior_ws) + log(like_ws[2]) - .5*log(2*pi * like_sds[2]^2 * prior_sds^2 * post_sds[(K+1):(2*K)]^-2) -
        .5 * like_mean[2]^2/like_sds[2]^2 - .5 * prior_means^2/prior_sds^2 +
        .5 * post_means[(K+1):(2*K)]^2 * post_sds[(K+1):(2*K)]^-2)

  # 2*K x 1 vector of posterior weights (normalized) corresponding to each component of the
  # posterior distribution for theta
  adj_log_post_w_propto <- exp(log_post_w_propto - max(log_post_w_propto))  # subtract max of log weights before exponentiating
  post_w_norm <- adj_log_post_w_propto / sum(adj_log_post_w_propto)
  dist_mixture(!!!post_ls, weights = correct_weights(post_w_norm))

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
