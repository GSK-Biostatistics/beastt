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
