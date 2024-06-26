#' Create a Propensity Score Object
#'
#' @description Calculate the propensity scores and ATT inverse probability
#'   weights for participants from internal and external datasets. Only the
#'   relevant treatment arms from each dataset should be read in (e.g., only
#'   the control arm from each dataset if creating a hybrid control arm).
#'
#' @param internal_df Internal dataset with one row per subject and all the
#'   variables needed to run the model
#' @param external_df External dataset with one row per subject and all the
#'   variables needed to run the model
#' @param id_col Name of the column in both datasets used to identify each
#'   subject. It must be the same across datasets
#' @param model Model used to calculate propensity scores
#' @param ... Optional arguments
#'
#' @details For the subset of participants in both the external and internal
#'    studies for which we want to balance the covariate distributions (e.g.,
#'    external control and internal control participants if constructing a
#'    hybrid control arm), we define a study-inclusion propensity score for
#'    each participant as
#'
#'    \deqn{e(x_i) = P(S_i = 1 \mid x_i),}
#'
#'    where \eqn{x_i} denotes a vector of baseline covariates for the \eqn{i}th
#'    participant and \eqn{S_i} denotes the indicator that the participant is
#'    enrolled in the internal trial (\eqn{S_i = 1} if internal, \eqn{S_i = 0}
#'    if external). The estimated propensity score \eqn{\hat{e}(x_i)} is obtained
#'    using logistic regression.
#'
#'    An ATT inverse probability weight is calculated for each individual as
#'
#'    \deqn{\hat{a}_{0i} = \frac{\hat{e}(x_i)}{\hat{P}(S_i = s_i | x_i)} = s_i + (1 - s_i ) \frac{\hat{e}(x_i)}{1 - \hat{e}(x_i)}.}
#'
#'    In a weighted estimator, data from participants in the external study
#'    are given a weight of \eqn{\hat{e}(x_i)⁄(1 - \hat{e}(x_i))} whereas data
#'    from participants in the internal trial are given a weight of 1.
#'
#' @return `prop_scr_obj` object, with the internal and the external data and
#'   the propensity score and inverse probability weight calculated for each
#'   subject.
#' @export
#' @importFrom rlang f_rhs f_lhs `f_lhs<-` is_formula enquo as_name
#' @importFrom stringr str_split_1 str_detect
#' @importFrom cli cli_abort
#' @importFrom purrr reduce discard
#' @importFrom rlang sym
#' @importFrom dplyr mutate filter tibble as_tibble
#' @importFrom stats glm
#' @importFrom cobalt bal.tab
#' @examples
#' # This can be used for both continuous and binary data
#' library(dplyr)
#' # Continuous
#' calc_prop_scr(internal_df = filter(int_norm_df, trt == 0),
#'                        external_df = ex_norm_df,
#'                        id_col = subjid,
#'                        model = ~ cov1 + cov2 + cov3 + cov4)
#' # Binary
#' calc_prop_scr(internal_df = filter(int_binary_df, trt == 0),
#'                        external_df = ex_binary_df,
#'                        id_col = subjid,
#'                        model = ~ cov1 + cov2 + cov3 + cov4)
#'
calc_prop_scr <- function(internal_df, external_df,
                            id_col, model, ...){
  if(!is_formula(model)){
    cli_abort(c("{.arg model} parameter must be a formula",
            "*" = "Make sure all covariates are left unquoted",
            "i" = "Formula should look like `~cov1 + cov2 + cov3...`"))
  } else {
    # Given model is a formula, check variable names match names in datasets

    if(!is.null(f_lhs(model))){
      cli_abort(c("The left hand side of the formula must be left blank",
           "i" = "Formula should look like `~cov1 + cov2 + cov3...`"))
    }

    covriates <- all.vars(model)

    int_dat <-inherits(internal_df, c("matrix", "data.frame"))
    ex_dat <-inherits(external_df, c("matrix", "data.frame"))
    if(!(int_dat & ex_dat)){
      cli_abort("{.agr int_dat} and {.arg ex_dat} both must either be tibbles, data.frames, or matrices")
    }

    inter <- reduce(list(colnames(internal_df),
                       colnames(external_df),
                       covriates),intersect)
    # Writing out error messages for if the variables are missing
    if(!all(covriates %in% inter)){
      in_miss <- setdiff(covriates, colnames(internal_df))
      ex_miss <- setdiff(covriates, colnames(external_df))
      if(length(in_miss) > 0 & length(ex_miss) > 0){
        cli_abort(c(
          "{.arg model} has covariates that are not present in the provided datasets",
          "x" = "The following variables are not in {.arg internal_df}: {in_miss}.",
          "x" = "The following variables are not in {.arg external_df}: {ex_miss}.",
          "i" = "Given both datasets are missing the covriate consider checking the spelling of the covariates"
        ))
      } else if(length(in_miss) > 0 ){
        cli_abort(c(
          "{.arg model} has covariates that are not present in {.arg internal_df}",
          "x" = "The following variables are missing: {in_miss}."))
      } else {
        cli_abort(c(
          "{.arg model} has covariates that are not present in {.arg external_df}",
          "x" = "The following variables are missing: {ex_miss}."))
      }
    }
  }

  id_col_en <- enquo(id_col)
  id_str <- id_col_en |> as_name()
  if(!id_str %in% colnames(internal_df) & !id_str %in% colnames(external_df)){
    cli_abort(c(
      "{.arg id_col} is in not the provided datasets",
      "i" = "Given both datasets are missing the variable consider checking the spelling of the column name"
    ))
  } else if(!id_str %in% colnames(internal_df)){
    cli_abort("{.arg id_col} is not in {.arg internal_df}")
  } else if( !id_str %in% colnames(external_df)){
    cli_abort("{.arg id_col} is not in {.arg external_df}")
  }

  # Adding internal as the response variable in the model
  f_lhs(model) <- sym('___internal___')

  i_df <- as_tibble(internal_df) |>
    mutate(`___internal___` = TRUE)
  e_df <- as_tibble(external_df) |>
    mutate(`___internal___` = FALSE)
  all_df <- bind_rows(i_df, e_df)


  all_df_ps <- all_df |>
    mutate(`___ps___` = glm(model, data = all_df, family = binomial)$fitted,
           `___weight___` = .data$`___internal___` + (1 - .data$`___internal___`) * .data$`___ps___` / (1 - .data$`___ps___`)
    )

  # Calculating the absolute standardized mean difference
  asmd_adj <- bal.tab(select(all_df_ps, !!covriates), # df of covariates (internal and external)
                  treat = all_df_ps$`___internal___`,   # internal indicator
                  binary = "std",         # use standardized version of mean differences for binary covariates
                  continuous = "std",     # use standardized version of mean differences for continuous covariates
                  s.d.denom = "pooled",   # calculation of the denominator of SMD
                  weights = all_df_ps$`___weight___`,
                  abs = TRUE)$Balance

  asmd_unadj <- bal.tab(select(all_df_ps, !!covriates), # df of covariates (internal and external)
                      treat = all_df_ps$`___internal___`,   # internal indicator
                      binary = "std",         # use standardized version of mean differences for binary covariates
                      continuous = "std",     # use standardized version of mean differences for continuous covariates
                      s.d.denom = "pooled",   # calculation of the denominator of SMD
                      abs = TRUE)$Balance

  asmd_clean <- tibble(
    covariate = rownames(asmd_adj),
    diff_unadj = asmd_unadj[,2],
    diff_adj = asmd_adj[,3],
  )


  structure(
    list(
         internal_df = filter(all_df_ps, .data$`___internal___` == TRUE),
         external_df = filter(all_df_ps, .data$`___internal___` == FALSE),
         id_col=id_col_en, model = model,
         abs_std_mean_diff = asmd_clean
         ),
    class = c("prop_scr")
  )
}

#' @export
#' @importFrom cli cli_h1 cli_text cli_bullets
#' @importFrom dplyr select
print.prop_scr <- function(x, ..., n = 10){
  # cat the model
  cli_h1("Model")
  cli_bullets(c("*" = f_rhs(x$model)))

  cli_h1("Propensity Scores and Weights")
  x$external_df |>
    select(!!x$id_col,
           Internal = .data$`___internal___`,
           `Propensity Score` = .data$`___ps___`,
           `Inverse Probability Weight` = .data$`___weight___`) |>
    print(n = n)

  cli_h1("Absolute Standardized Mean Difference")
  print(x$abs_std_mean_diff)
}

#' @export
tidy.prop_scr <- function(x, ...){
  x$external_df |>
    select(!!x$id_col,.data$`___internal___`, .data$`___ps___`, .data$`___weight___`)
}

#' @export
glance.prop_scr <- function(x, ...){
  x$external_df |>
    bind_rows(x$internal_df)
}


#' Test If Propensity Score Object
#'
#' @param x Object to test
#'
#' @return Boolean
#' @export
#' @examples
#' library(dplyr)
#' x <- calc_prop_scr(internal_df = filter(int_norm_df, trt == 0),
#'                        external_df = ex_norm_df,
#'                        id_col = subjid,
#'                        model = ~ cov1 + cov2 + cov3 + cov4)
#' is_prop_scr(x)
#'
is_prop_scr <- function(x){
  inherits(x, "prop_scr")
}

#' Testing prop_scr
#'
#' Internal function what will return an error if the object is not a propensity
#' score
#'
#' @param x Object to test
#'
#' @return NULL
#' @noRd
#' @keywords internal
#' @importFrom cli cli_abort
test_prop_scr <- function(x){
  if(!is_prop_scr(x)){
    cli_abort(c("{.var {substitute(x)}} is not a `prop_scr` object",
                "i" = "Please use {.fun calc_prop_scr} to make a `prop_scr` object"))
  }
}


#' Histogram of the Propensity Score Object
#'
#' @description Plot overlapping histograms of the propensity scores for both
#'   the internal and external participants, or plot external IPWs.
#'
#' @param x Propensity score object
#' @param variable Variable to plot. It must be either a propensity score
#'   ("ps" or "propensity score") or inverse probability weight ("ipw" or
#'   "inverse probability weight")
#' @param ... Optional arguments for `geom_histogram`
#'
#' @return ggplot object
#' @export
#' @importFrom ggplot2 ggplot aes geom_histogram labs scale_fill_manual ggtitle
#'    theme_bw
#' @importFrom dplyr bind_rows
#' @importFrom stringr str_glue
#' @examples
#' library(dplyr)
#' ps_obj <- calc_prop_scr(internal_df = filter(int_norm_df, trt == 0),
#'                        external_df = ex_norm_df,
#'                        id_col = subjid,
#'                        model = ~ cov1 + cov2 + cov3 + cov4)
#' # Plotting the Propensity Scores
#' prop_scr_hist(ps_obj)
#' # Or plotting the inverse probability weights
#' prop_scr_hist(ps_obj, variable = "ipw")
prop_scr_hist <- function(x, variable = c("propensity score", "ps", "inverse probability weight", "ipw"),
                          ...){
  test_prop_scr(x)

  plot_var <- match.arg(variable)
  x_var <- ifelse(plot_var %in% c("propensity score", "ps"),
                          sym('___ps___'),
                          sym('___weight___'))
  x_label <- ifelse(plot_var %in% c("propensity score", "ps"),
                     "Propensity Score",
                     "Inverse Probability Weight")

  if(plot_var %in% c("propensity score", "ps")) {
    .data <- bind_rows(x$internal_df, x$external_df)
  } else {
    .data <- x$external_df
  }

  plot <-   .data|>
    ggplot(aes(x = !!x_var, fill = .data$`___internal___`)) +
    labs(y = "", x = x_label, fill = "Dataset") +
    scale_fill_manual(values = c("#FFA21F", "#5398BE"),
                      labels = c("TRUE" =  "Internal", "FALSE" = "External")) +
    ggtitle(str_glue("Histogram of {x_label}s")) +
    theme_bw()

  if(length(list(...)) == 0) {
    plot <- plot +
      geom_histogram(position = "identity", binwidth = .05)
  } else {
    plot <- plot +
      geom_histogram(position = "identity", ...)
  }

  plot

}


#' Density of the Propensity Score Object
#'
#' @description Plot overlapping density curves of the propensity scores for
#'   both the internal and external participants, or plot external IPWs.
#'
#' @param x Propensity score object
#' @param variable Variable to plot. It must be either a propensity score
#'   ("ps" or "propensity score") or inverse probability weight ("ipw" or
#'   "inverse probability weight")
#' @param ... Optional arguments for `geom_density`
#'
#' @return ggplot object
#' @export
#' @importFrom ggplot2 ggplot aes geom_density labs scale_fill_manual ggtitle
#'   theme_bw
#' @importFrom dplyr bind_rows
#' @importFrom stringr str_glue
#' @examples
#' library(dplyr)
#' ps_obj <- calc_prop_scr(internal_df = filter(int_norm_df, trt == 0),
#'                        external_df = ex_norm_df,
#'                        id_col = subjid,
#'                        model = ~ cov1 + cov2 + cov3 + cov4)
#' # Plotting the Propensity Scores
#' prop_scr_dens(ps_obj)
#' # Or plotting the inverse probability weights
#' prop_scr_dens(ps_obj, variable = "ipw")
#'
prop_scr_dens <- function(x, variable = c("propensity score", "ps", "inverse probability weight", "ipw"),
                          ...){
  test_prop_scr(x)

  plot_var <- match.arg(variable)
  x_var <- ifelse(plot_var %in% c("propensity score", "ps"),
                  sym('___ps___'),
                  sym('___weight___'))
  x_label <- ifelse(plot_var %in% c("propensity score", "ps"),
                    "Propensity Score",
                    "Inverse Probability Weight")

  if(plot_var %in% c("propensity score", "ps")) {
    .data <- bind_rows(x$internal_df, x$external_df)
  } else {
    .data <- x$external_df
  }

  plot <-   .data |>
    ggplot(aes(x = !!x_var, fill = .data$`___internal___`)) +
    labs(y = "Density", x = x_label, fill = "Dataset") +
    scale_fill_manual(values = c("#FFA21F", "#5398BE"),
                      labels = c("TRUE" =  "Internal", "FALSE" = "External")) +
    ggtitle(str_glue("Density of {x_label}s")) +
    theme_bw()

  if(length(list(...)) == 0) {
    plot <- plot +
      geom_density(alpha = 0.5)
  } else {
    plot <- plot +
      geom_density(...)
  }

  plot

}


#' Love Plot of the Absolute Standardized Mean Differences
#'
#' @description Plot the unadjusted and IPW-adjusted absolute standardized mean
#'   differences for each covariate.
#'
#' @param x Propensity score object
#' @param reference_line Numeric value of where along the x-axis the vertical
#'   reference line should be placed
#' @param ... Optional options for `geom_point`
#'
#'
#' @return ggplot object
#' @export
#' @importFrom ggplot2 ggplot aes geom_point labs scale_color_manual ggtitle
#'   theme_bw geom_vline
#' @importFrom tidyr pivot_longer
#' @examples
#' library(dplyr)
#' ps_obj <- calc_prop_scr(internal_df = filter(int_norm_df, trt == 0),
#'                        external_df = ex_norm_df,
#'                        id_col = subjid,
#'                        model = ~ cov1 + cov2 + cov3 + cov4)
#' # Plotting the Propensity Scores
#' prop_scr_love(ps_obj, reference_line = 0.1)
#'
prop_scr_love <- function(x, reference_line = NULL, ...){
  test_prop_scr(x)

  .data <- x$abs_std_mean_diff |>
    pivot_longer(-"covariate")

  plot <- .data |>
    ggplot(aes(x = .data$value, color = .data$name,
               y = .data$covariate)) +
    labs(y = "Covariates", x = "Absolute Standardized Mean Difference", color = "Sample") +
    scale_color_manual(values = c("#FFA21F", "#5398BE"),
                      labels = c("diff_unadj" =  "Unadjusted", "diff_adj" = "Adjusted")) +
    ggtitle("Covariates Balance") +
    geom_point() +
    theme_bw()

  if(!is.null(reference_line)){
    plot <- plot + geom_vline(xintercept = reference_line)
  }
  plot

}



