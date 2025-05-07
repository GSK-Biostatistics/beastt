#' Create Sweet Spot Plots for Multiple Simulation Scenarios
#'
#' Creates visualization plots to help identify the "sweet spot" in borrowing
#' strategies across different simulation scenarios. For each unique scenario
#' defined by the combination of variables in `scenario_vars`, the function
#' produces a plot showing power, Type I error, and power prior distribution.
#'
#' @param .data A data frame containing iteration-level simulation results.
#' @param scenario_vars A vector of unquoted column names that define the unique
#'   scenarios to create plots for. Each unique combination of values in these
#'   columns will generate a separate plot.
#' @param trt_diff An unquoted column name representing the treatment
#'   difference. Used to identify scenarios with null effect (trt_diff = 0) for
#'   Type I error calculation.
#' @param control_marg_param An unquoted column name to be used as the x-axis in the plots.
#'   This is typically the control end point of interest ton the marginal scale (e.g. Control Response Rate)
#' @param power_prior An unquoted column name containing distributional objects
#'   that represent the power prior distribution.
#' @param h0_prob An unquoted column name containing the probability of
#'   accepting the null hypothesis when using borrowing.
#' @param h0_prob_no_borrowing An unquoted column name containing the
#'   probability of accepting the null hypothesis without borrowing.
#'
#' @return A list of ggplot objects, one for each unique scenario defined by
#'   `scenario_vars`. Each plot shows:
#'   \itemize{
#'     \item Power curves with (and without borrowing if `h0_prob_no_borrowing` is given)
#'     \item Type I error rates with (and without borrowing if `h0_prob_no_borrowing` is given)
#'     \item Average Distribution of the power priors
#'   }
#'
#' @details The function calculates power and Type I error rates for each
#' scenario and visualizes them together with the average power prior distribution. This
#' helps identify the "sweet spot" where borrowing provides power gains while
#' maintaining acceptable Type I error rates.
#'
#' Type I error is calculated using scenarios where `trt_diff` equals 0. Power
#' is calculated across all scenarios.
#'
#'
#' @examples
#'
#' # Assuming sim_results is a data frame with simulation results in the shape of binary template code
#' plots <- sweet_spot_plot(
#'   .data = binary_sim_df,
#'   scenario_vars = c("population", "marg_trt_eff"),
#'   trt_diff = marg_trt_eff,
#'   control_marg_param = true_control_RR,
#'   power_prior = pwr_prior,
#'   h0_prob = reject_H0_yes,
#'   h0_prob_no_borrowing = no_borrowing_reject_H0_yes
#' )
#'
#' # Display the first plot
#' plots[[1]]
#'
#'
#' @export
sweet_spot_plot <- function(.data, scenario_vars,
                             trt_diff, control_marg_param,
                             power_prior,
                             h0_prob, h0_prob_no_borrowing
                             ){

  data_grouped <- .data |>
    dplyr::group_by(dplyr::across({{scenario_vars}}), {{control_marg_param}})

  # Calculate type 1 Error for all scenarios
  type1_df <- data_grouped |>
    dplyr::filter({{trt_diff}} == 0) |>
    dplyr::summarise(type1_borrowing = mean({{h0_prob}}),
              type1_no_borrowing = mean({{h0_prob_no_borrowing}}),
              .groups = "drop") |>
    dplyr::select({{scenario_vars}}, {{control_marg_param}}, .data$type1_borrowing,
                  .data$type1_no_borrowing, -{{trt_diff}}) |>
    tidyr::pivot_longer(c(.data$type1_borrowing, .data$type1_no_borrowing),
                        names_prefix = "type1_",
                        names_to = "borrowing_status", values_to = "Type 1 Error")

  # Calculate power for all scenarios
  power <- data_grouped |>
    dplyr::summarise(h0_prob_borrowing = mean({{h0_prob}}),
              h0_prob_no_borrowing = mean({{h0_prob_no_borrowing}}),
              .groups = "drop") |>
    dplyr::select({{scenario_vars}}, {{control_marg_param}}, .data$h0_prob_borrowing,
                  .data$h0_prob_no_borrowing) |>
    tidyr::pivot_longer(c(.data$h0_prob_borrowing, .data$h0_prob_no_borrowing),
                        names_prefix = "h0_prob_",
                        names_to = "borrowing_status", values_to = "Power")

  # Remove the treatment difference from the scenario vector if there
  diff_col_str <- rlang::as_string(rlang::ensym(trt_diff))
  sc_vars_no_trt_diff <- scenario_vars |>
    purrr::discard(\(x) x == diff_col_str)
  # Combine power and type 1 error, then nest each dataset down so there should
  # only be one row per scenario. Dropping where the trt difference is 0
  plot_df <- power |>
    dplyr::left_join(type1_df, by = dplyr::join_by({{control_marg_param}}, !!!sc_vars_no_trt_diff, "borrowing_status")) |>
    dplyr::filter({{trt_diff}} != 0) |>
    tidyr::pivot_longer(c("Power", .data$`Type 1 Error`)) |>
    dplyr::group_by(dplyr::across({{scenario_vars}})) |>
    tidyr::nest()

  # Get the average power prior for each scenario,
  # dropping where the trt difference is 0
  prior <- .data |>
    dplyr::group_by(dplyr::across({{scenario_vars}})) |>
    dplyr::filter({{trt_diff}} != 0) |>
    dplyr::summarise(pwr_prior = avg_dist({{power_prior}}), .groups = "drop_last") |>
    dplyr::select({{scenario_vars}}, .data$pwr_prior)

  # Create a plot fro each scenario
  plot_ls<- plot_df |>
    dplyr::left_join(prior, dplyr::join_by(!!!scenario_vars)) |>
    purrr::pmap(\(...){
      inputs <- list(...)

      scenarios_cols <- inputs[!names(inputs) %in% c("data", "pwr_prior")]
      title = purrr::map2_chr(names(scenarios_cols), scenarios_cols,
                      \(x, y) stringr::str_c(x, y, sep = ": ")) |>
        stringr::str_c(collapse = ", ")

      x_vals <- dplyr::pull(inputs$data, {{control_marg_param}})
      colors <- c("#5398BE", "#FFA21F")

      ggplot() +
        ggdist::stat_slab(aes(xdist = inputs$pwr_prior, fill = "Power Prior"),
                          alpha = 0.6,
                          color = "grey80", size = 0.25,
                          # show.legend = TRUE
                          ) +
        ggplot2::geom_line(data = inputs$data,
                           aes(x = {{control_marg_param}}, y = .data$value,
                               linetype = .data$borrowing_status, color = .data$name),
                           size = 0.75) +
        ggplot2::scale_y_continuous(name = "Power",
                                    sec.axis = ggplot2::sec_axis(~ ., name = "Type 1 Error")) +
        ggplot2::scale_x_continuous(limits = c(min(x_vals), max(x_vals))) +
        ggplot2::scale_color_manual(values = colors, name = "") +
        ggplot2::scale_linetype_manual(name = "",
                                       values = c("borrowing" = "solid",
                                                  "no_borrowing" = "dashed"),
                                       labels = c("Borrowing", "No Borrowing")) +
        ggplot2::scale_fill_manual(name = "", values = "grey85") +
        ggplot2::guides(fill = ggplot2::guide_legend(order = 1),  # Set the order of the fill legend
               color = ggplot2::guide_legend(order = 2),  # Set the order of the color legend
               linetype = ggplot2::guide_legend(order = 3)) + # Set the order of the linetype legend
        ggplot2::ggtitle(title) +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.position = "bottom")

    })
  plot_ls

}



#' Calculate Average Distribution from Multiple Distributional Objects
#'
#' Computes a single "average" distribution from a vector of distributional
#' objects. This function calculates the mean of each parameter across all
#' input distributions and returns a new distributional object with these
#' averaged parameters.
#'
#' @param x A vector of distributional objects of the same family (beta,
#'   normal, or multivariate normal).
#'
#' @return A single distributional object of the same family as the input,
#'   with parameters set to the average of all input distribution parameters.
#'
#' @details
#' The function supports three distribution families:
#' \itemize{
#'   \item Beta distributions: Averages the shape1 and shape2 parameters
#'   \item Normal distributions: Averages the mean and sd parameters
#'   \item Multivariate normal distributions: Averages the location vectors
#'     and scale matrices
#' }
#'
#' For multivariate normal distributions, both the location vector and scale
#' matrix are averaged element-wise.
#'
#'
#' @examples
#' library(distributional)
#'
#' # Beta distributions
#' beta_dists <- c(
#'   dist_beta(shape1 = 2, shape2 = 5),
#'   dist_beta(shape1 = 3, shape2 = 3),
#'   dist_beta(shape1 = 4, shape2 = 2)
#' )
#' avg_dist(beta_dists)
#'
#' # Normal distributions
#' norm_dists <- c(
#'   dist_normal(mu = 0, sigma = 1),
#'   dist_normal(mu = 2, sigma = 2),
#'   dist_normal(mu = 4, sigma = 3)
#' )
#' avg_dist(norm_dists) |> parameters()
#'
#' # Multivariate normal distributions
#' mvn_dists <- c(
#'   dist_multivariate_normal(mu = list(c(0, 0)), sigma = list(matrix(c(1, 0, 0, 1), nrow = 2))),
#'   dist_multivariate_normal(mu = list(c(1, 1)), sigma = list(matrix(c(2, 0, 0, 2), nrow = 2)))
#' )
#'
#' avg_dist(mvn_dists) |> parameters()
#' @export
avg_dist <- function(x){
  if(!distributional::is_distribution(x)){
    cli_abort("{.arg x} must be a distributional object")
  }
  all_fam <- unique(family(x))
  if(length(all_fam) > 1){
    cli_abort("All distributional objects in {.agr x} must be of the same type")
  } else if(!all_fam %in% c("beta", "mvnorm", "normal")){
    cli_abort("Only beta, normal and multivariate normal distributions are supported")
  }
  dist_params <- distributional::parameters(x)
  if(all_fam == "mvnorm"){
    avg_params <- purrr::map(dist_params, \(y){
      purrr::list_transpose(y) |> purrr::map_dbl(mean)
    })

    avg_params$sigma <- avg_params$sigma |>
      matrix(ncol = sqrt(length(avg_params$sigma)))
    avg_params <- purrr::map(avg_params, list)
  } else {
    avg_params <- colMeans(dist_params) |> as.list()
  }

  dist_fx <- switch(all_fam,
                    "beta" = distributional::dist_beta,
                    "mvnorm" = distributional::dist_multivariate_normal,
                    "normal" = distributional::dist_normal)
  do.call(dist_fx,avg_params)


}
