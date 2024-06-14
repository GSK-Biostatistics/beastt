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
plot_dist <- function(...){
  input <- list(...)
  if(!all(purrr::map_lgl(input, is_distribution))){
    cli_abort("Must be given distributional objects")
  }

  Distributions <- names(input)
  n <- length(input)
  fill_alpha <- ifelse(n > 1, 0.5, 1)
  if(is.null(Distributions) & n > 1){
    Distributions <- purrr::map_chr(input, format)
    if(unique(Distributions) != n)
      Distributions <- paste0(1:n,": ", Distributions)

  }
  ggplot(data.frame(), aes(xdist = input, fill = Distributions)) +
    stat_slabinterval(alpha = fill_alpha) +
    labs(y = "Density", x = "") +
    theme_bw()

}
