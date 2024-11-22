#' Write Code function
#'
#' @param simulation A purpose, either "Simulation" or "Analysis"
#' @param endpoint An endpoint type, one of "Binary", "Normal" or "Survival"
#' @param selections input from UI
#'
#' @importFrom rstudioapi documentNew
#' @importFrom dplyr case_match
#' @importFrom stringr str_glue
write_code <- function(simulation, endpoint, input_list){

  doc <- list(
    header = "########################################################
#################  {beastt} Template  #################\n
########################################################\n
library(beastt)\nlibrary(distributional)\nlibrary(dplyr)
    ",

    data_setup = "",
    priors = "",
    post = ""
  )

  if(simulation == "Simulation"){
    doc <- write_simulation_sect(doc, endpoint, input_list)
  } else {
    doc$data_setup <- "################# Read in Data ################# \n
internal_df <- read.csv('DATA LOCATION')\nexternal_df <- read.csv('DATA LOCATION')\n
    "
  }

  out <- doc |>
    write_prior_sect(endpoint, input_list) |>
    write_post_sect(endpoint, input_list) |>
    paste0(collapse = "\n")

  documentNew(out, type = "r")


}


write_simulation_sect <- function(doc, endpoint, input_list){
  doc
}

write_prior_sect <- function(doc, endpoint, input_list){
  if(input_list$analysis_inputs$borrType == "On control arm"){
    int_dat <- "filter(internal_df, trt == 0), #EDIT 'trt == 0' TO SELECT ONLY THE CONTROL ARM"
  } else {
    int_dat <- "internal_df"
  }
  if(length(input_list$analysis_inputs$plots$plotProp) > 0){
    prop_plots <- case_match(input_list$analysis_inputs$plots$plotProp,
                             "Histogram" ~ 'prop_scr_hist(ps_obj)',
                             "Histogram - IPW" ~ 'prop_scr_hist(ps_obj, variable = "ipw")',
                             "Density" ~ 'prop_scr_dens(ps_obj)',
                             "Density - IPW" ~ 'prop_scr_dens(ps_obj, variable = "ipw")',
                             "Love" ~ 'prop_scr_love(ps_obj, reference_line = 0.1)'
    ) |>
      paste0(collapse = "\n")


  } else {
    prop_plots <- ""
  }

  prior <- switch (endpoint,
                   "Binary" = "dist_beta(0.5, 0.5)",
                   "Normal" = "dist_normal(0, 10)"
  )

  if(input_list$analysis_inputs$robustify){
    robust_prior <- str_glue("mix_prior <- dist_mixture(pwr_prior, {prior}, weights = c(0.5, 0.5))")
  } else {
    robust_prior <- ""
  }

  if(length(input_list$analysis_inputs$plots$plotPrior) > 0){
    priors_to_plot <- case_match(input_list$analysis_inputs$plots$plotPrior,
                                 "Vague" ~ str_glue('"Vague Prior" = {prior}'),
                                 "Power Prior" ~ '"Power Prior" = pwr_prior',
                                 "Robust Mixture" ~ str_glue('"Robustified Power Prior" = mix_prior')
    ) |>
      paste0(collapse = ", \n")
    prior_plots <- str_glue('plot_dist({priors_to_plot})')

  } else {
    prior_plots <- ""
  }

  doc$priors <- str_glue(
    "# Calculate Propensity Scores by creating a prop_scr object
    ps_obj <- calc_prop_scr(internal_df = {int_dat},
                        external_df = external_df,
                        id_col = subjid, # EDIT 'subjid' IF COLUMN NAME IS DIFFERENT
                        model = ~ cov1 + cov2 + cov3 + cov4) #EDIT TO THE CORRECT COVARIATES
    ps_obj

    # Create plots related to Propensity Scores
    {prop_plots}

    # Calculate power prior
    pwr_prior <- calc_power_prior_beta(ps_obj,
                                   response = y,
                                   prior = {prior})
    {robust_prior}

    # Create prior plots
    {prior_plots}
    "
  )
  doc

}

write_post_sect <- function(doc, endpoint, input_list){
  if(input_list$analysis_inputs$robustify){
    post <- "post_mixed <- calc_post_beta(ps_obj, response = y, prior = mix_prior)"
    post_to_plot <- '"Robust Mixture Posterior" = post_mixed'
    prior_to_plot <- '"Robust Mixture Prior" = mix_prior'
  } else {
    post <- "post <- calc_post_beta(ps_obj, response = y, prior = pwr_prior)"
    post_to_plot <- '"Posterior" = post'
    prior_to_plot <- '"Power Prior" = pwr_prior'
  }

  if(length(input_list$analysis_inputs$plots$plotPost) > 0){
    posts_to_plot <- case_match(input_list$analysis_inputs$plots$plotPost,
                                "Prior" ~ str_glue('plot_dist({prior_to_plot}, {post_to_plot})'),
                                "Posterior" ~ str_glue('plot_dist({post_to_plot})')
    ) |>
      paste0(collapse = "\n")
    post_plots <- str_glue('{posts_to_plot}')

  } else {
    post_plots <- ""
  }
  doc$posts <- str_glue(
    "# Calculate posterior distribution
  {post}

  # Create posterior plots
  {post_plots}
  "
  )
  doc
}
