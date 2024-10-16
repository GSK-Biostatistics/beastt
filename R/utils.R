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
#' @importFrom ggdist stat_slabinterval
#' @importFrom purrr map_chr map_lgl
#' @importFrom grDevices rainbow
#' @examples
#' library(distributional)
#' plot_dist(dist_normal(0, 1))
#' #Plotting Multiple
#' plot_dist(dist_normal(0, 1), dist_normal(10, 5))
#' plot_dist('Prior' = dist_normal(0, 1), 'Posterior' = dist_normal(10, 5))
plot_dist <- function(...){
  input <- list(...)
  if(!all(map_lgl(input, is_distribution))){
    cli_abort("Must be given distributional objects")
  }

  Distributions <- names(input)
  n <- length(input)
  if(n > 2) {
    colors <- c("#5398BE", "#FFA21F", rainbow(n-2))

  } else {
    colors <- c("#5398BE", "#FFA21F")
  }

  fill_alpha <- ifelse(n > 1, 0.5, 1)
  if(is.null(Distributions) & n > 0){
    Distributions <- map_chr(input, format)
    if(length(unique(Distributions)) != n)
      Distributions <- paste0(1:n,": ", Distributions)

  }

  ggplot(data.frame(), aes(xdist = input, fill = Distributions)) +
    stat_slabinterval(alpha = fill_alpha) +
    labs(y = "Density", x = "") +
    scale_fill_manual(values = colors) +
    theme_bw()

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
#' @param prior Normal distributional object
#' @param n Number of theoretical participants
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
  prior_checks(prior, "normal")
  prior_param <- parameters(prior)
  robust_prior <- dist_mixture(prior,
                               dist_normal(prior_param$mu, sqrt(prior_param$sigma^2*n)),
                               weights = weights
  )
  robust_prior
}
