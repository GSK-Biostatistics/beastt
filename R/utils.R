#' Plot Distribution
#'
#' @param ... Distributional object(s) to plot. If there are multiple objects you can name them.
#'
#' @return ggplot object that is the density of the provided distribution
#' @export
#'
#' @importFrom distributional is_distribution
#' @importFrom ggplot2 ggplot theme_bw
#' @importFrom ggdist stat_slabinterval
#' @importFrom purrr map_chr map_lgl
plot_dist <- function(...){
  input <- list(...)
  if(!all(map_lgl(input, is_distribution))){
    cli_abort("Must be given distributional objects")
  }

  Distributions <- names(input)
  n <- length(input)
  fill_alpha <- ifelse(n > 1, 0.5, 1)
  if(is.null(Distributions) & n > 1){
    Distributions <- map_chr(input, format)
    if(unique(Distributions) != n)
      Distributions <- paste0(1:n,": ", Distributions)

  }
  ggplot(data.frame(), aes(xdist = input, fill = Distributions)) +
    stat_slabinterval(alpha = fill_alpha) +
    labs(y = "Density", x = "") +
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
