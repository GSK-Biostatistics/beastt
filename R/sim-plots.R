#' Create Sweet Spot Plots for Multiple Simulation Scenarios
#'
#' Create visualization plots to help identify the "sweet spot" in borrowing
#' strategies across different simulation scenarios. For each unique scenario
#' defined by the combination of variables in `scenario_vars`, the function
#' produces a plot showing power, type I error, and the distribution of the
#' design prior for the control marginal parameter for approaches with and
#' without borrowing.
#'
#' @param .data A data frame containing iteration-level simulation results.
#' @param scenario_vars A vector of quoted column names corresponding
#'   to variables used to define unique simulation scenarios. Each unique
#'   combination of values in these columns will generate a separate plot.
#' @param trt_diff An unquoted column name representing the treatment
#'   difference. Used to identify scenarios with null effect (trt_diff = 0) for
#'   type I error calculation.
#' @param control_marg_param An unquoted column name to be used as the x-axis in
#'   the plots. This is typically the control endpoint of interest on the
#'   marginal scale (e.g., control response rate).
#' @param design_prior An unquoted column name containing distributional objects
#'   that represent the design prior distribution for the control marginal
#'   parameter (e.g., posterior distribution using the external control data).
#'   Used to aid visualization of which values of the control marginal parameter
#'   are assumed to be plausible. Default is `NULL`, in which case no design
#'   prior is plotted. See Details for more information.
#' @param h0_prob An unquoted column name containing the probability of
#'   rejecting the null hypothesis when when borrowing external data.
#' @param h0_prob_no_borrowing An unquoted column name containing the
#'   probability of rejecting the null hypothesis when not borrowing
#'   external data.
#'
#' @return A list of ggplot objects, one for each unique scenario defined by
#'   `scenario_vars`. Each plot shows:
#'   \itemize{
#'     \item Power curves for the cases with and without borrowing
#'     \item Type I error rates for the cases with and without borrowing
#'     \item Distribution of the design prior (if `design_prior` is specified)
#'   }
#'
#' @details The function calculates power and type I error rates for BDB approaches
#' that borrow from external data (e.g., use of a robust mixture prior with positive
#' weight on the informative component) and an approach that does not
#' borrow from external data (e.g., use of a vague prior) under each scenario
#' and visualizes them together as a function of the underlying control marginal
#' parameter of interest (e.g., control response rate for binary outcomes) that
#' may vary as a result of drift. This helps identify the "sweet spot" where borrowing
#' results in higher power and lower type I error rates compared to not borrowing.
#' Type I error is calculated using scenarios where `trt_diff` equals 0, and power
#' is calculated for all scenarios with positive values of `trt_diff`.
#'
#' If `design_prior` is non-`NULL`, the design prior distribution is included
#' in the plot to provide insight into which values of the control marginal
#' parameter are plausible given this assumed design prior. We note that
#' `design_prior` can represent any informative prior that potentially
#' incorporates the external control data (e.g., the posterior distribution of
#' the control marginal parameter constructed using the external data and a
#' vague prior). Each element of the vector corresponding to `design_prior` must
#' be a distributional object with a family equal to "beta", "normal", or
#' "mixture" (where each component is either "beta" or "normal"). For the
#' time-to-event case in which a multivariate normal prior is assumed for the
#' control log-shape and intercept of a Weibull proportional hazards model,
#' this distribution must first be translated into a univariate beta design
#' prior for the control survival probability at some prespecified time.
#' This approximation can be done using [approx_mvn_at_time()]. If the
#' design priors in the vector indicated by `design_prior` differ across
#' iterations within a given scenario (e.g., using the IPW power prior as the
#' iteration-specific design prior), then the average distribution will be
#' plotted (i.e., a distribution of the same family with the hyperparameters
#' averaged across iterations).
#'
#' @references
#' Best, N., Ajimi, M., Neuenschwander, B., Saint-Hilary, G., & Wandel, S.
#' (2024). Beyond the Classical Type I Error: Bayesian Metrics for Bayesian
#' Designs Using Informative Priors. \emph{Statistics in Biopharmaceutical Research},
#' 17(2), 183–196. \doi{10.1080/19466315.2024.2342817}
#'
#' @examples
#'
#' # Assuming binary_sim_df is a data frame with simulation results in the shape of binary template code
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
                        names_to = "borrowing_status", values_to = "Type I Error")

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
    tidyr::pivot_longer(c("Power", .data$`Type I Error`)) |>
    dplyr::group_by(dplyr::across({{scenario_vars}})) |>
    tidyr::nest()

  # Get the average power prior for each scenario,
  # dropping where the trt difference is 0
  prior <- .data |>
    dplyr::group_by(dplyr::across({{scenario_vars}})) |>
    dplyr::filter({{trt_diff}} != 0) |>
    dplyr::summarise(pwr_prior = avg_dist({{power_prior}}), .groups = "drop_last") |>
    dplyr::select({{scenario_vars}}, .data$pwr_prior)

  # Create a plot for each scenario
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
                                    sec.axis = ggplot2::sec_axis(~ ., name = "Type I Error")) +
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
#' Compute a single "average" distribution from a vector of distributional
#' objects. This function calculates the mean of each hyperparameter across all
#' input distributions and returns a new distributional object of the same
#' family with these averaged hyperparameters.
#'
#' @param x A vector of distributional objects of the same family (beta,
#'   normal, multivariate normal, or mixture).
#'
#' @return A single distributional object of the same family as the input,
#'   with hyperparameters set equal to the average of all input distribution
#'   hyperparameters.
#'
#' @details
#' The function supports four distribution families:
#' \itemize{
#'   \item Beta distributions: Averages the shape1 and shape2 hyperparameters
#'   \item Normal distributions: Averages the mean and standard deviation
#'     hyperparameters
#'   \item Multivariate normal distributions: Averages the location vectors
#'     and covariance matrices
#'   \item Mixture distributions: Same as above for each distribution type,
#'     where averaging is done by component. Also averages the mixture weight.
#' }
#'
#' For multivariate normal distributions, both the location vector and covariance
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
#' avg_dist(beta_dists) |> parameters()
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
#' avg_dist(mvn_dists) |> parameters()
#'
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
