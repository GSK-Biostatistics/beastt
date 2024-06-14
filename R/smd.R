#' Calculate Absolute Standardized Mean Differences
#'
#' @description Calculate the unadjusted and IPW-adjusted absolute standardized
#'   mean differences for each covariate.
#'
#'
#' @param prop_scr_obj Propensity score object
#' @param treatment_var Treatment variable
#' @param ... Any additional variables
#'
#' @return Tibble of the calculated standardized mean differences for the
#'   additional variables
#' @export
#' @importFrom rlang enquo enquos as_label
calc_std_mean_diff <- function(prop_scr_obj, treatment_var, ...){
  test_prop_scr(prop_scr_obj)
  dots <- enquos(...)
  dot_str <- dots |> map(rlang::as_label) |> unlist()
  trt_var <- enquo(treatment_var)

  df <- bind_rows(prop_scr_obj$internal_df,prop_scr_obj$external_df)

  trt_check <- safely(select)(df, !!trt_var)
  if(!is.null(trt_check$error)){
    cli_abort(
      "{.arg treatment_var} variable is not present in the provided dataset.
      Please ensure that {.arg prop_scr_obj} has {!!trt_var} ")
  }
  dot_check <- safely(select)(df, !!!dots)
  if(!is.null(dot_check$error)){
    missing<- setdiff(dot_str, colnames(df))  |>
      paste(collapse = ",")
    cli_abort(
      "{.arg ...} has variable(s) that are not present in the provided dataset.
      Please ensure that {.arg prop_scr_obj} has the following variables {missing} ")
  }

  covs <- all.vars(f_rhs(prop_scr_obj$model))
  new_vars <- setdiff(dot_str, covs) |>
    unique()


  asmd_adj <- bal.tab(select(df, !!!new_vars), # df of covariates (internal and external)
                      treat = df$`___internal___`,   # internal indicator
                      binary = "std",         # use standardized version of mean differences for binary covariates
                      continuous = "std",     # use standardized version of mean differences for continuous covariates
                      s.d.denom = "pooled",   # calculation of the denominator of SMD
                      weights = df$`___weight___`,
                      abs = TRUE)$Balance

  asmd_unadj <- bal.tab(select(df, !!!new_vars), # df of covariates (internal and external)
                        treat = df$`___internal___`,   # internal indicator
                        binary = "std",         # use standardized version of mean differences for binary covariates
                        continuous = "std",     # use standardized version of mean differences for continuous covariates
                        s.d.denom = "pooled",   # calculation of the denominator of SMD
                        abs = TRUE)$Balance
  asmd_clean <- tibble(
    covariate = rownames(asmd_adj),
    diff_unadj = asmd_unadj[,2],
    diff_adj = asmd_adj[,3],
  )
  browser()
  needed_cols <- c(trt_var, sym(`___internal___`), prop_scr_obj$id_col)
  grped_df <- select(df, !!!c(new_vars, covs, trt_var), `___internal___`) |>
    dplyr::group_by(!!trt_var, `___internal___`)
  cat_vars <- grped_df |>
    select(!!trt_var, `___internal___`, where(\(x) is.character(x)|is.factor(x))) |>
    mutatae(across(where(\(x) is.character(x)|is.factor(x)), ))

  cont_vars <- grped_df |>
    select(!!trt_var, `___internal___`, where(is.numeric))

  prop_scr_obj$abs_std_mean_diff |>
    bind_rows(asmd_clean)



}
