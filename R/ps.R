#' Create a Proprensity Score Object
#'
#' @param internal_df Internal dataset with one row per subject and all the
#'   variables needed to run the model
#' @param external_df External dataset with one row per subject and all the
#'   variables needed to run the model
#' @param id_col Name of the column in both datasets used to identify each
#'   subject. It must be the same across datasets
#' @param model Model used to calculate propensity scores
#' @param ...
#'
#' @return `prop_scr_obj` object, with the internal and the external data and
#'   the propensity score and inverse probability weight calculated for each
#'   subject.
#' @export
#' @importFrom rlang f_rhs f_lhs `f_lhs<-` is_formula enquo as_name
#' @importFrom stringr str_split_1
#' @importFrom cli cli_abort
#' @importFrom purrr reduce
create_prop_scr <- function(internal_df, external_df,
                            id_col, model, ...){
  if(!is_formula(model)){
    cli_abort(c("{.arg model} parameter must be a fomula",
            "*" = "Make sure all covariates are left unquoted",
            "i" = "Formula should look like '~cov1 + cov2 + cov3...`"))
  } else {
    # Given model is a fomula, check variable names match names in datasets

    if(!is.null(f_lhs(model))){
      cli_abort(c("The left hand side of the formula must be left blank",
           "i" = "Formula should look like '~cov1 + cov2 + cov3...`"))
    }

    covriates <- f_rhs(model) |>
      format() |>
      str_split_1("\\s?\\+\\s?")

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

  i_df <- internal_df |> mutate(`___internal___` = TRUE)
  e_df <- external_df |> mutate(`___internal___` = FALSE)
  all_df <- bind_rows(i_df, e_df)


  all_df_ps <- all_df |>
    mutate(`___ps___` = glm(model, data = all_df, family = binomial)$fitted,
           `___ipw___` = `___internal___` + (1 - `___internal___`) * `___ps___` / (1 - `___ps___`)
    )

  structure(
    list(
         internal_df = filter(all_df_ps, `___internal___` == TRUE),
         external_df = filter(all_df_ps, `___internal___` == FALSE),
         id_col=id_col_en, model = model
         ),
    class = c("prop_scr_obj")
  )
}

#' @export
#' @importFrom cli cli_h1 cli_text
print.prop_scr_obj <- function(x, n = 10){
  # cat the model
  cli_h1("Model")
  cli_bullets(c("*" = f_rhs(x$model)))
  cli_h1("Propensoity Scores and Weights")
  x$external_df |>
    select(!!x$id_col, `Propensity Score` = `___ps___`,
           `Inverse Probablity Weight` = `___ipw___`) |>
    print(n = n)
}

#' @export
tidy.prop_scr_obj <- function(x){
  x$external_df |>
    select(!!x$id_col, `___ps___`, `___ipw___`)
}

#' @export
glance.prop_scr_obj <- function(x){
  x$external_df |>
    bind_rows(x$internal_df)
}

# TODO add plots

