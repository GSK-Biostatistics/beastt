#' Calculate Posterior Normal
#'
#' @description Calculate a posterior distribution that is normal (or a mixture
#'   of normal components). Only the relevant treatment arms from the internal
#'   dataset should be read in (e.g., only the control arm if constructing a
#'   posterior distribution for the control mean).
#'
#' @param internal_data A tibble of the internal data.
#' @param response Name of response variable
#' @param prior A distributional object corresponding to a normal distribution,
#'   a t distribution, or a mixture distribution of normal and/or t components
#' @param internal_sd Standard deviation of internal response data if
#'   assumed known. It can be left as `NULL` if assumed unknown
#'
#' @details For a given arm of an internal trial (e.g., the control arm or an
#'   active treatment arm) of size \eqn{N_I}, suppose the response data are normally
#'   distributed such that \eqn{Y_i \sim N(\theta, \sigma_I^2)}, \eqn{i=1,\ldots,N_I}.
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
#'                                internal_sd = 0.15)
#'
calc_post_norm<- function(
    internal_data,
    response,
    prior,
    internal_sd = NULL
){
  # Checking internal data and response variable
  if(is.data.frame(internal_data)) {
    data <- internal_data
    nIC <- nrow(internal_data)
  } else{
    cli_abort("{.agr internal_data} a dataset")
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
#' @param internal_data A tibble of the internal data.
#' @param response Name of response variable
#' @param prior A distributional object corresponding to a beta distribution
#'   or a mixture distribution of beta components
#'
#' @details For a given arm of an internal trial (e.g., the control arm or an
#' active treatment arm) of size \eqn{N_I}, suppose the response data are binary
#' such that \eqn{Y_i \sim \mbox{Bernoulli}(\theta)}, \eqn{i=1,\ldots,N_I}. The
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
#'                response = y,
#'                prior = dist_beta(0.5, 0.5))
calc_post_beta<- function(internal_data, response, prior){
  # Checking internal data and response variable
  if(is.data.frame(internal_data)) {
    data <- internal_data
    nIC <- nrow(internal_data)
  } else{
    cli_abort("{.agr internal_data} a dataset")
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



#' Calculate Posterior Weibull
#'
#' @description Calculate a posterior distribution for the probability of
#'   surviving past a given analysis time(s) for time-to-event data with a
#'   Weibull likelihood. Only the relevant treatment arms from the internal
#'   dataset should be read in (e.g., only the control arm if constructing a
#'   posterior distribution for the control survival probability).
#'
#' @param internal_data This can either be a propensity score object or a tibble
#'   of the internal data.
#' @param response Name of response variable
#' @param event Name of event indicator variable (1: event; 0: censored)
#' @param prior A distributional object corresponding to a multivariate normal
#'   distribution or a mixture of 2 multivariate normals. The first element
#'   of the mean vector and the first row/column of covariance matrix correspond
#'   to the log-shape parameter, and the second element corresponds to the intercept
#'   (i.e., log-inverse-scale) parameter.
#' @param analysis_time Vector of time(s) when survival probabilities will be
#'   calculated
#' @param ... rstan sampling option. This overrides any default options. For more
#'   information, see [rstan::sampling()]
#'
#' @details For a given arm of an internal trial (e.g., the control arm or an
#'   active treatment arm) of size \eqn{N_I}, suppose the response data are time to event
#'   such that \eqn{Y_i \sim \mbox{Weibull}(\alpha, \sigma)}, where
#'
#'   \deqn{f(y_i \mid \alpha, \sigma) = \left( \frac{\alpha}{\sigma} \right) \left( \frac{y_i}{\sigma}
#'   \right)^{\alpha - 1} \exp \left( -\left( \frac{y_i}{\sigma} \right)^\alpha
#'   \right),}
#'
#'   \eqn{i=1,\ldots,N_I}. Define \eqn{\boldsymbol{\theta} = \{\log(\alpha), \beta\}}
#'   where \eqn{\beta = -\log(\sigma)} is the intercept parameter of a Weibull
#'   proportional hazards regression model. The posterior distribution for
#'   \eqn{\boldsymbol{\theta}} is written as
#'
#'   \deqn{\pi( \boldsymbol{\theta} \mid \boldsymbol{y}, \boldsymbol{\nu} ) \propto
#'   \mathcal{L}(\boldsymbol{\theta} \mid \boldsymbol{y}, \boldsymbol{\nu}) \;
#'   \pi(\boldsymbol{\theta}),}
#'
#'   where \eqn{\mathcal{L}(\boldsymbol{\theta} \mid \boldsymbol{y}, \boldsymbol{\nu}) =
#'   \prod_{i=1}^{N_I} f(y_i \mid \boldsymbol{\theta})^{\nu_i} S(y_i \mid \boldsymbol{\theta})^{1 - \nu_i}}
#'   is the likelihood of the response data from the internal arm with event indicator
#'   \eqn{\boldsymbol{\nu}} and survival function \eqn{S(y_i \mid \boldsymbol{\theta}) =
#'   1 - F(y_i \mid \boldsymbol{\theta})}, and \eqn{\pi(\boldsymbol{\theta})} is a prior
#'   distribution on \eqn{\boldsymbol{\theta}} (either a multivariate normal distribution or a
#'   mixture of two multivariate normal distributions). Note that the posterior distribution
#'   for \eqn{\boldsymbol{\theta}} does not have a closed form, and thus MCMC samples
#'   for \eqn{\log(\alpha)} and \eqn{\beta} are drawn from the posterior. These MCMC
#'   samples are used to construct samples from the posterior distribution
#'   for the probability of surviving past the specified analysis time(s) for the
#'   given arm.
#'
#'   To construct a posterior distribution for the treatment difference (i.e., the
#'   difference in survival probabilities at the specified analysis time), obtain MCMC
#'   samples from the posterior distributions for the survival probabilities under
#'   each arm and then subtract the control-arm samples from the treatment-arm
#'   samples.
#'
#' @return stan posterior object
#' @export
#' @importFrom rstan sampling
#'
#'
calc_post_weibull <- function(internal_data,
                              response, event,
                              prior,
                              analysis_time,
                              ...){
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

  # Check response exists in the data
  response <- enquo(response)
  check <- safely(select)(data, !!response)
  if(!is.null(check$error)){
    cli_abort("{.agr response} was not found in {.agr internal_data}")
  }

  # Check event indicator exists in the data
  event <- enquo(event)
  check <- safely(select)(data, !!event)
  if(!is.null(check$error)){
    cli_abort("{.agr event} was not found in {.agr internal_data}")
  }

  # Checking the distribution and getting the family
  if(!is_distribution(prior)){
    cli_abort("{.agr prior} must be a distributional object")
  }
  prior_fam <- family(prior)

  if(prior_fam == "mvnorm"){
    mean <- parameters(prior)$mu[[1]]
    cov <- parameters(prior)$sigma[[1]]
    stan_data_post_inputs <- list(
      N = nIC,   # internal sample size of given arm
      y = pull(data, !!response),       # observed time (event or censored)
      e = pull(data, !!event),   # event indicator (1: event; 0: censored)
      K = length(analysis_time),              # number of times to compute survival probabilities
      times = as.array(analysis_time),        # time(s) when survival probability should be calculated
      mean_inf = mean,      # mean vector of informative prior component
      mean_vague = mean,  # mean vector of vague prior component
      cov_inf = cov,        # covariance matrix of informative prior component
      cov_vague = cov,    # covariance matrix of vague prior component
      w0 = 0                          # prior mixture weight associated with informative component
    )

  } else if(prior_fam == "mixture" &&  length(get_base_families(prior)[[1]]) == 2) {
    means <- mix_means(prior)
    covs <- mix_sigmas(prior)
    weight <- parameters(prior)$w[[1]][1] # Weight of robust mixture prior (RMP) associated with informative component

    stan_data_post_inputs <- list(
      N = nIC,   # internal sample size of given arm
      y = pull(data, !!response),       # observed time (event or censored)
      e = pull(data, !!event),   # event indicator (1: event; 0: censored)
      K = length(analysis_time),              # number of times to compute survival probabilities
      times = as.array(analysis_time),        # time(s) when survival probability should be calculated
      mean_inf = means[[1]],      # mean vector of informative prior component
      mean_vague = means[[2]],  # mean vector of vague prior component
      cov_inf = covs[[1]],        # covariance matrix of informative prior component
      cov_vague = covs[[2]],    # covariance matrix of vague prior component
      w0 = weight                        # prior mixture weight associated with informative component
    )



  } else {
    cli_abort("{.agr prior} must be a multivariate normal distributional object or a mixture of 2 multivariate normal objects")
    }
  post <- sampling_optional_inputs(list(stanmodels$weibullpost,
                                        data = stan_data_post_inputs,
                                        warmup = 15000,       # number of burn-in iterations to discard
                                        iter = 45000,         # total number of iterations (including burn-in)
                                        chains = 1,           # number of chains
                                        cores = 1,            # number of cores
                                        refresh = 0),
                                   ...)
  post

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



#' Extract Means of Mixture Components
#'
#' @param x A mixture distributional object
#'
#' @details If a distributional object that is a mixture of two or more normal
#'   distributions is read in, the function will return a numeric object with
#'   the means of each normal component. If the distributional object is a
#'   mixture of two or more multivariate normal distributions, the function
#'   will return a list with the mean vectors of each multivariate normal
#'   component.
#'
#' @return numeric or list object
#' @export
#'
#' @examples
#' library(distributional)
#' mix_norm <- dist_mixture(comp1 = dist_normal(1, 10),
#'                          comp2 = dist_normal(1.5, 12),
#'                          weights = c(.5, .5))
#' mix_means(mix_norm)
mix_means <- function(x){
  out <- parameters(x)$dist[[1]]|>
    map(\(dist) dist$mu)
  if(length(unlist(out)) == length(out)){
    out <- out |> unlist()
  }
  out

}



#' Extract Standard Deviations of Mixture Components
#'
#' @param x A mixture distributional object
#'
#' @details If a distributional object that is a mixture of two or more normal
#'   distributions is read in, the function will return a numeric object with
#'   the standard deviations of each normal component. If the distributional
#'   object is a mixture of two or more multivariate normal distributions, the
#'   function will return a list with the covariance matrices of each multivariate
#'   normal component.
#'
#' @return numeric or list object
#' @export
#'
#' @examples
#' library(distributional)
#' mix_norm <- dist_mixture(comp1 = dist_normal(1, 10),
#'                          comp2 = dist_normal(1.5, 12),
#'                          weights = c(.5, .5))
#' mix_sigmas(mix_norm)
mix_sigmas <- function(x){
  out <- parameters(x)$dist[[1]]|>
    map(\(dist) dist$sigma)
  if(length(unlist(out)) == length(out)){
    out <- out |> unlist()
  }
  out
}
