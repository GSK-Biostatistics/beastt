#' Plot Distribution
#'
#' @param ... Distributional object(s) to plot. When passing multiple objects
#'   naming them will change the labels in the plot, else they will use the
#'   distributional format
#'
#' @return ggplot object that is the density of the provided distribution
#' @export
#'
#' @importFrom distributional is_distribution
#' @importFrom ggplot2 ggplot theme_bw
#' @importFrom ggdist stat_slab
#' @importFrom purrr map_chr map_lgl
#' @importFrom grDevices rainbow
#' @examples
#' library(distributional)
#' plot_dist(dist_normal(0, 1))
#' plot_dist(dist_multivariate_normal(mu = list(c(1, 2)), sigma = list(matrix(c(4, 2, 2, 3), ncol=2))))
#' #Plotting Multiple
#' plot_dist(dist_normal(0, 1), dist_normal(10, 5))
#' plot_dist('Prior' = dist_normal(0, 1), 'Posterior' = dist_normal(10, 5))
plot_dist <- function(...){
  input <- list(...)
  if(!all(map_lgl(input, is_distribution))){
    cli_abort("Must be given distributional objects")
  }

  if(any(map_lgl(input, is_mvnorm))){
    plot_dist_mvnorm(input)
  } else {
    Distributions <- names(input)
    n <- length(input)
    colors <- c("#5398BE", "#FFA21F", rainbow(n-2))

    if(is.null(Distributions) & n > 0){
      Distributions <- map_chr(input, format, width = 2)
      if(length(unique(Distributions)) != n)
        Distributions <- paste0(1:n,": ", Distributions)
    }
    plot <- ggplot(data.frame(), aes(xdist = input)) +
      stat_slab(aes(color = Distributions), fill = NA) +
      labs(y = "Density", x = "") +
      scale_fill_manual(values = colors) +
      scale_color_manual(values = colors) +
      theme_bw()

    if(n==1){
      plot + stat_slab(aes(fill = Distributions), alpha = 0.5)
    } else {
      plot
    }

  }
}

#' Plot multivariate normal distribution contours
#'
#' @param dist_list list of distributional objects
#'
#' @noRd
#' @return ggplot
#'
#' @importFrom ggplot2 geom_contour scale_color_manual
#' @importFrom dplyr rowwise mutate
#' @importFrom tidyr crossing
#' @importFrom stats density
#' @importFrom purrr map_lgl map_dbl
#' @importFrom grDevices rainbow
plot_dist_mvnorm <- function(dist_list){

  if(!all(map_lgl(dist_list, is_mvnorm))){
    cli_abort("Must be given multivariate normal objects")
  }
  check_pkg_installed("mvtnorm")
  tib <- data.frame(x1=numeric(), x2=numeric(), z=numeric(), id=factor())
  for (i in 1:length(dist_list)){
    dist <- dist_list[[i]]
    mus <- mean(dist)
    vars <- variance(dist)
    if(length(mus) == 1 & length(vars) == 1){
      cli_abort("Cannot plot multivariate normals with an index of 1")
    }
    sigmas <- map_dbl(vars, sqrt)
    x1 <- seq(from=mus[1]-2.5*sigmas[1], to=mus[1]+2.5*sigmas[1], length.out=100)
    x2 <- seq(from=mus[2]-2.5*sigmas[2], to=mus[2]+2.5*sigmas[2], length.out=100)
    temp <- crossing(x1, x2) |>
      rowwise() |>
      mutate(z = unlist(density(dist, c(x1, x2))), id=as.factor(i))
    tib <- rbind(tib, temp)
  }

  n <- length(dist_list)
  if(n == 1){
    ggplot(tib, aes(x=x1, y=x2, z=.data$z)) +
      geom_contour(colour = "black") +
      labs(x=expression(paste("log(", alpha, ")")), y=expression(beta)) +
      theme_bw()
  } else {
    colors <- c("#5398BE", "#FFA21F", rainbow(n-2))
    ggplot(tib, aes(x=x1, y=x2, z=.data$z, colour=.data$id)) +
      geom_contour() +
      labs(x=expression(paste("log(", alpha, ")")), y=expression(beta)) +
      theme_bw() + scale_color_manual(values = colors)
  }
}

#' Check if a distributional object belongs to mvnorm family
#'
#' @param dist a distributional object
#'
#' @noRd
#' @return logical value
is_mvnorm <- function(dist){
  ifelse(family(dist)=="mvnorm", TRUE, FALSE)
}


#' Correct weights to always sum to 1
#'
#' @param x vector of weights
#'
#' @noRd
#' @return vector of weights that sum to 1
correct_weights <- function(x){
  n <- length(x)
  x[n] <- 1 - sum(x[1:(n-1)])
  x
}

#' Robustify Normal Distributions
#'
#' @description Adds vague normal component, where the level of vagueness is controlled by
#' the `n` parameter
#' @param prior Normal or Multivariate Normal distributional object
#' @param n Number of theoretical participants (or events, for time-to-event data)
#' @param weights Vector of weights, where the first number corresponds to the
#'   informative component and the second is the vague
#'
#' @details In cases with a normal endpoint, a robust mixture prior can be created by
#'    adding a vague normal component to any normal prior with mean \eqn{\theta}
#'    and variance \eqn{\sigma^2}.The vague component is calculated to have the
#'    same mean \eqn{\theta} and variance equal to \eqn{\sigma^2 \times n}, where
#'    `n` is the specified number of theoretical participants. If robustifying a normal
#'    power prior that was calculated from external control data and `n` is defined as
#'    the number of external control participants, and the vague component would
#'    then correspond to one external control participant's worth of data.
#'
#' @return mixture distribution
#' @export
#'
#' @examples
#' library(distributional)
#' robustify_norm(dist_normal(0,1), n = 15)
robustify_norm <- function(prior, n, weights = c(0.5, 0.5)){
  prior_fam <- family(prior)
  if(prior_fam == "normal"){
    prior_param <- parameters(prior)
    robust_prior <- dist_mixture(informative = prior,
                                 vague = dist_normal(prior_param$mu, sqrt(prior_param$sigma^2*n)),
                                 weights = weights
    )
    robust_prior
  } else if(prior_fam == "mvnorm") {
    robust_prior <- robustify_mvnorm(prior, n, weights)
  } else {
    cli_abort("{.agr prior} must be either normal or a multivariate normal")
  }
  robust_prior
}

#' Robustify Multivariate Normal Distributions
#'
#' @description Adds vague normal component, where the level of vagueness is controlled by
#' the `n` parameter
#' @param prior Multivariate Normal distributional object
#' @param n Number of theoretical participants (or events, for time-to-event data)
#' @param weights Vector of weights, where the first number corresponds to the
#'   informative component and the second is the vague
#'
#' @details In cases with a time-to-event endpoint, a robust mixture prior can be
#'    created by adding a vague multivariate normal component to any multivariate
#'    normal prior with mean vector \eqn{\boldsymbol{\mu}} and covariance matrix
#'    \eqn{\boldsymbol{\Sigma}}. The vague component is calculated to have the
#'    same mean vector \eqn{\boldsymbol{\mu}} and covariance matrix equal to
#'    \eqn{\boldsymbol{\Sigma} \times n}, where `n` is the specified number of
#'    theoretical events.
#'
#' @return mixture distribution
#' @export
#'
#' @examples
#' library(distributional)
#' robustify_mvnorm(
#'       dist_multivariate_normal(mu = list(c(1, 0)), sigma = list(c(10, 5))),
#'        n = 15)
robustify_mvnorm <- function(prior, n, weights = c(0.5, 0.5)){
  prior_checks(prior, "mvnorm")
  prior_param <- parameters(prior)
  robust_prior <- dist_mixture(informative = prior,
                               vague = dist_multivariate_normal(prior_param$mu,
                                                        list(prior_param$sigma[[1]]*n)),
                               weights = weights
  )
  robust_prior
}

#' Sampling with updated inputs
#'
#' @param defaults list of default options for the stan model
#' @param ...
#'
#' @return stan sample
#' @noRd
sampling_optional_inputs <- function(defaults, ...){
  input_vals <- list(...)
  defaults[names(input_vals)] <- input_vals
  stan_samp <- do.call(sampling, defaults)
}
