
write_code <- function(simulation, endpoint, selections){
  doc <- list(
    header = "########################################################
################# {beastt} Template ################# \n
########################################################\n

library(beastt)\nlibrary(distributional)\nlibrary(dplyr)
    ",

    data_setup = "",
    priors = "",
    post = ""
  )

  if(simulation == "Simulation"){
    doc <- write_simulation_sect(doc, endpoint, selection)
  } else {
    doc$data_setup <- "################# Read in Data ################# \n
internal_df <- read.csv('DATA LOCATION')\nexternal_df <- read.csv('DATA LOCATION')\n
    "
  }

  out <- doc |>
    write_prior_sect(endpoint, selections) |>
    write_post_sect(endpoint, selections) |>
    paste0(collapse = "\n")

  rstudioapi::documentNew(out, type = "r")


}


write_simulation_sect <- function(doc, endpoint, selections){
  doc
}

write_prior_sect <- function(doc, endpoint, selections){
  if(selections$borrType == "On control arm"){
    int_dat <- "filter(internal_df, trt == 0), #EDIT 'trt == 0' TO SELECT ONLY THE CONTROL ARM"
  } else {
    int_dat <- "internal_df"
  }
  if(length(selections$plots$plotProp) > 0){
    prop_plots <- dplyr::case_match(selections$plots$plotProp,
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

  pwer_fx <- switch(endpoint,
                    "Binary" = "calc_power_prior_beta",
                    "Normal" = "calc_power_prior_norm"
                    )


  doc$priors <- stringr::str_glue(
    "# Calculate Propensity Scores by creating a prop_scr object
    ps_obj <- calc_prop_scr(internal_df = {int_dat},
                        external_df = external_df,
                        id_col = subjid, # EDIT 'subjid' IF COLUMN NAME IS DIFFERENT
                        model = ~ cov1 + cov2 + cov3 + cov4) #EDIT TO THE CORRECT COVARIATES
    ps_obj
    {prop_plots}
    pwr_prior <- {pwer_fx}(ps_obj,
                                   response = y,
                                   prior = {prior})
    "
  )
  doc

}

write_post_sect <- function(doc, endpoint, selection){
  doc
}
