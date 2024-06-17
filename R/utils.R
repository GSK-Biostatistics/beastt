#' Title
#'
#' @param distribution Distributional object to plot
#'
#' @return ggplot object that is the density of the provided distribution
#' @export
#'
#' @importFrom distributional is_distribution
#' @importFrom ggplot2 ggplot theme_bw
#' @importFrom ggdist stat_slabinterval
plot_dist <- function(distribution){
  if(!is_distribution(distribution)){
    cli_abort("Must be given a distributional object")
  }

  ggplot(data.frame(), aes(xdist = distribution, y =1)) +
    stat_slabinterval() +
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
