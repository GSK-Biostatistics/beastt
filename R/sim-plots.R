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
#'   probability of rejecting the null hypothesis when not borrowing external
#'   data.
#' @param highlight Logical value to indicate if you want sweet spot
#'   highlighting or not. If `TRUE` the sweet spot (where borrowing increase
#'   power and reduces type 1 error) will be highlighted.
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
#'
#' @examples
#' library(dplyr)
#' # Assuming binary_sim_df is a data frame with simulation results in the shape of binary template code
#' plots <- sweet_spot_plot(
#'   .data = binary_sim_df,
#'   scenario_vars = c("population", "marg_trt_eff"),
#'   trt_diff = marg_trt_eff,
#'   control_marg_param = true_control_RR,
#'   prior = pwr_prior,
#'   h0_prob = reject_H0_yes,
#'   h0_prob_no_borrowing = no_borrowing_reject_H0_yes
#' )
#'
#' # Display the first plot
#' plots[[1]]
#'
#' tte_plots <- tte_sim_df |>
#'  mutate(beta_appox = approx_mvn_at_time(mix_prior, time = 12)) |>
#'  sweet_spot_plot(
#'    scenario_vars = c("population", "marg_trt_eff"),
#'    trt_diff = marg_trt_eff,
#'    control_marg_param = true_control_surv_prob,
#'    prior = beta_appox,
#'    h0_prob = reject_H0_yes,
#'    h0_prob_no_borrowing = no_borrowing_reject_H0_yes
#'  )
#'
#' tte_plots[[1]]
#'
#' @export
sweet_spot_plot <- function(.data, scenario_vars,
                            trt_diff, control_marg_param,
                            design_prior = NULL,
                            h0_prob, h0_prob_no_borrowing,
                            highlight = TRUE
){

  if(nrow(dplyr::filter(.data, {{trt_diff}} == 0)) == 0){
    cli_abort("Unable to calculate Type 1 Error without a scenario where `trt_diff` equals 0")
  }

  data_grouped <- .data |>
    dplyr::group_by(dplyr::across({{scenario_vars}}),
                    {{control_marg_param}},
                    {{trt_diff}})

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
    dplyr::select({{scenario_vars}}, {{control_marg_param}}, {{trt_diff}},
                  .data$h0_prob_borrowing, .data$h0_prob_no_borrowing) |>
    tidyr::pivot_longer(c(.data$h0_prob_borrowing, .data$h0_prob_no_borrowing),
                        names_prefix = "h0_prob_",
                        names_to = "borrowing_status", values_to = "Power")

  # Remove the treatment difference from the scenario vector if there
  sc_vars_no_trt_diff <- .data |>
    dplyr::ungroup() |>
    dplyr::select({{scenario_vars}}, {{trt_diff}}) |>
    dplyr::select(-{{trt_diff}}) |>
    colnames()

  # Combine power and type 1 error, then nest each dataset down so there should
  # only be one row per scenario. Dropping where the trt difference is 0
  plot_df <- power |>
    dplyr::left_join(type1_df, by = dplyr::join_by({{control_marg_param}}, !!!sc_vars_no_trt_diff, "borrowing_status")) |>
    dplyr::filter({{trt_diff}} != 0) |>
    tidyr::pivot_longer(c("Power", .data$`Type I Error`)) |>
    dplyr::group_by(dplyr::across({{scenario_vars}})) |>
    tidyr::nest()

  design_prior_test <- !rlang::quo_is_null(rlang::enquo(design_prior))

  if(design_prior_test){
    prior_col <- .data |> dplyr::pull({{design_prior}})
    if(!distributional::is_distribution(prior_col)){
      cli_abort("`design_prior` must be a column of distributional objects")
    }
    all_fam <- unique(family(prior_col))
    if(all_fam %in% c("mvnorm")){
      cli_abort("Multivariate `design_prior` need to be approximated as a beta, see `approx_mvn_at_time()`")
    }

    # Get the average design prior for each scenario,
    # dropping where the trt difference is 0
    prior <- .data |>
      dplyr::group_by(dplyr::across({{scenario_vars}})) |>
      dplyr::filter({{trt_diff}} != 0) |>
      dplyr::summarise(des_prior = avg_dist({{design_prior}}), .groups = "drop_last") |>
      dplyr::select({{scenario_vars}}, .data$des_prior)

    # Create a plot for each scenario
    quite_join <- purrr::quietly(dplyr::left_join)
    plot_df <- plot_df |>
      quite_join(prior) |>
      _$result
  }

  plot_ls<-  plot_df |>
    purrr::pmap(\(...){
      inputs <- list(...)

      scenarios_cols <- inputs[!names(inputs) %in% c("data", "des_prior")]
      title = purrr::map2_chr(names(scenarios_cols), scenarios_cols,
                              \(x, y) {
                                y_one <- stringr::str_c(y, collapse = ", ")
                                stringr::str_c(x, y_one, sep = ": ")

                              }) |>
        stringr::str_c(collapse = ", ")

      x_vals <- dplyr::pull(inputs$data, {{control_marg_param}})
      colors <- c("#5398BE", "#FFA21F")

      type1_range <- inputs$data |>
        dplyr::filter(.data$name == "Type 1 Error") |>
        dplyr::pull() |>
        range()

      power_range <- inputs$data |>
        dplyr::filter(.data$name == "Power") |>
        dplyr::pull() |>
        range()

      # For when there is an inflated type 1 error
      inflate_fct = min(power_range / type1_range)/2
      inflate_fct = ifelse(all(power_range > type1_range), inflate_fct, 1)

      scaled_df <- inputs$data |>
        dplyr::mutate(
          value = dplyr::case_when(
            .data$name == "Type 1 Error" ~  .data$value*inflate_fct,
            TRUE ~ .data$value
          )
        )

      if(design_prior_test){
        plot <- ggplot() +
          ggdist::stat_slab(aes(xdist = inputs$des_prior, fill = "Design Prior"),
                            alpha = 0.5,
                            color = "grey80", linewidth = 0.25
          )
      } else {
        plot <- ggplot()
      }

      plot <- plot +
        ggplot2::geom_line(data = scaled_df,
                           aes(x = {{control_marg_param}}, y = .data$value,
                               linetype = .data$borrowing_status, color = .data$name),
                           linewidth = 0.75) +
        ggplot2::scale_y_continuous(name = "Power",
                                    sec.axis = ggplot2::sec_axis(~ ./inflate_fct, name = "Type 1 Error")) +
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

      if(highlight){
        # Reshape to make borrowing and no-borrowing seperate columns
        reshaped_df <- scaled_df |>
          tidyr::pivot_wider(id_cols = c({{control_marg_param}}, name),
                      names_from = borrowing_status)

        # For each point, calculate the slope from the previous point to that one
        # Then calculate the intercept to determine when the borrowing and no borrowing line would cross
        line_cross_df <- reshaped_df |>
          dplyr::group_by(name) |>
          dplyr::mutate(
            dplyr::across(c(borrowing, no_borrowing), \(x){
              slope = (x - dplyr::lag(x)) / ({{control_marg_param}} - dplyr::lag({{control_marg_param}}))
            }, .names = "{.col}_slope"),
            int_borrowing = borrowing - borrowing_slope*{{control_marg_param}},
            int_no_borrowing = no_borrowing - no_borrowing_slope*{{control_marg_param}},
            line_cross = (int_borrowing - int_no_borrowing) / (no_borrowing_slope - borrowing_slope))

        # Get the minimum point when borrowing is greater than no borrowing for type 1 and power
        # These values represent the min and max value
        edge_vec <- line_cross_df |>
          dplyr::filter(name == "Power",
                        {{control_marg_param}} %in% c(min({{control_marg_param}}), max({{control_marg_param}}))
                        ) |>
          dplyr::pull(borrowing)
        dirction_test <- ifelse(edge_vec[1] < edge_vec[2], "positive", "negative")

        check_fx <- switch(dirction_test,
                           "positive" = min,
                           "negative" = max)
        highlight_range <- line_cross_df |>
          dplyr::group_by(name) |>
          dplyr::filter(borrowing > no_borrowing) |>
          dplyr::filter({{control_marg_param}} == check_fx({{control_marg_param}})) |>
          dplyr::select(name, line_cross) |>
          tidyr::pivot_wider(names_from = name, values_from = line_cross)

        sweet_spot_check <- ifelse(dirction_test == "positive",
                      highlight_range$Power > highlight_range$`Type 1 Error`,
                      highlight_range$Power < highlight_range$`Type 1 Error`)
        if((is.na(highlight_range$Power) | is.na(highlight_range$`Type 1 Error`)) ||
          sweet_spot_check
          ) {
          cli::cli_warn("No sweet spot avaliable to highlight")
        } else {
          plot <- plot +
            ggplot2::geom_rect(aes(xmin = highlight_range$Power,
                               xmax = highlight_range$`Type 1 Error`,
                               ymin = 0, ymax = 1),
                               alpha = 0.25, fill = "#93A646"
                               )

          }
      }
      plot
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
#'   For multivariate normal distributions, both the location vector and
#'   covariance matrix are averaged element-wise.
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
  } else if(!all_fam %in% c("beta", "mvnorm", "normal", "mixture")){
    cli_abort("Only beta, normal, multivariate normal distributions and mixtures of those are supported")
  }
  dist_params <- distributional::parameters(x)
  if(all_fam == "mvnorm"){
    avg_params <- purrr::map(dist_params, \(y){
      purrr::list_transpose(y) |> purrr::map_dbl(mean)
    })

    avg_params$sigma <- avg_params$sigma |>
      matrix(ncol = sqrt(length(avg_params$sigma)))
    avg_params <- purrr::map(avg_params, list)
  } else if(all_fam == "mixture"){
    avg_weights <-  do.call(rbind, dist_params$w) |>
      colMeans()
    comp_dists <- c(1:length(avg_weights))|>
      map(\(index){
        dist_list <- dist_params$dist|>
          purrr::map(\(dist){
            dist[[index]]
          })
        class(dist_list) <- c("distribution", "vctrs_vctr", class(dist_list))
        avg_dist(dist_list)
      })
    avg_params <- c(comp_dists, weights = list(correct_weights(avg_weights)))

  } else {
    avg_params <- colMeans(dist_params) |> as.list()
  }
  dist_fx <- switch(all_fam,
                    "beta" = distributional::dist_beta,
                    "mvnorm" = distributional::dist_multivariate_normal,
                    "normal" = distributional::dist_normal,
                    "mixture" = distributional::dist_mixture)
  do.call(dist_fx,avg_params)

}


#' Approximate Multivariate Normal Distribution as Beta at a Specific Time
#'
#' Converts a multivariate normal distribution for Weibull parameters (or a mixture
#' of these distributions) into an approximate beta distribution for the survival
#' probability at a specific time point. This is particularly useful for visualizing
#' survival probabilities in sweet spot plots
#'
#' @param x A vector of distributional objects that must be either multivariate normal
#'   distributions or mixtures of multivariate normal distributions. For Weibull models,
#'   these represent distributions of the log(shape) and log(scale) parameters.
#' @param time A numeric value specifying survival time at which to calculate the
#'   survival probability.
#'
#' @return A vector of beta distributional (or mixture of beta distributional)
#'   objects approximating the survival probabilities at the specified time
#'   point. If the input is a mixture distribution, the output will be a mixture
#'   of beta distributions with the same weights.
#'
#' @seealso [sweet_spot_plot()]
#'
#' @details
#' The function performs the following steps:
#' \itemize{
#'   \item For each multivariate normal distribution, it generates 10,000 samples of the
#'     Weibull parameters
#'   \item Calculates the corresponding survival probabilities at the specified time
#'     using the Weibull survival function
#'   \item Fits a beta distribution to match the mean and variance of these survival
#'     probabilities
#'   \item For mixture distributions, it performs this approximation for each component
#'     and creates a new mixture with the same weights
#' }
#'
#' The conversion uses the relationship between Weibull parameters and survival
#' probability: S(t) = exp(-(t*exp(log_scale))^exp(log_shape)).
#'
#'
#' @examples
#'
#' library(distributional)
#'
#' # Create a multivariate normal distribution for Weibull parameters
#' # (log(shape), log(scale))
#' mvn_dist <- dist_multivariate_normal(
#'   mu = list(c(0, -1)),  # log(shape) = 0, log(scale) = -1
#'   sigma = list(matrix(c(0.1, 0, 0, 0.1), nrow = 2))
#' )
#'
#' # Approximate as beta distribution for survival at time t=12
#' beta_approx <- approx_mvn_at_time(mvn_dist, time = 12)
#'
#'
#' @export
approx_mvn_at_time <- function(x, time){
  if(!distributional::is_distribution(x)){
    cli_abort("{.arg x} must be a distributional object")
  }
  all_fam <- unique(family(x))
  if(!all_fam %in% c("mvnorm", "mixture")){
    cli_abort("Only multivariate normal distributions and mixtures of those are supported")
  }
  dist_list <- purrr::map(x, \(prior){
    if(family(prior) == "mvnorm"){
      samples <- distributional::generate(prior, times = 10000)[[1]]
      # Sample from IPW power prior for control survival probability at time
      surv_prob_pp <- exp(-(time * exp(samples[,2]))^exp(samples[,1]) )
      # Approximate the IPW power prior for the control survival probability at time
      # t with a beta distribution
      IPW_pp_mean <- mean(surv_prob_pp)
      IPW_pp_var <- var(surv_prob_pp)
      IPW_pp_shape1 <- ((1 - IPW_pp_mean) / IPW_pp_var - 1 / IPW_pp_mean) * IPW_pp_mean^2
      IPW_pp_shape2 <- IPW_pp_shape1 * (1 / IPW_pp_mean - 1)
      dist_beta(shape1 = IPW_pp_shape1, shape2 = IPW_pp_shape2)
    } else {
      prior_params <- distributional::parameters(prior)
      dists <- prior_params$dist[[1]]
      class(dists) <- c("distribution", "vctrs_vctr", "list")
      appox_beta_list <-  approx_mvn_at_time(dists, time = time)
      dist_mixture(list(appox_beta_list), weights = prior_params$w[[1]])
    }
  })
  class(dist_list) <- c("distribution", "vctrs_vctr", class(dist_list))
  dist_list
}
