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
#'
create_prop_scr <- function(internal_df, external_df,
                            id_col, model, ...){

  # TODO Add tests for checking data and model are consistent
  # TODO makes it so the model is put in with just ~ x + bx

  i_df <- internal_df |> mutate(`___internal___` = TRUE)
  e_df <- external_df |> mutate(`___internal___` = FALSE)
  all_df <- bind_rows(i_df, e_df)


  all_df_ps <- all_df |>
    mutate(`___ps___` = glm(model, data = all_df, family = binomial)$fitted,
           `___ipw___` = `___internal___` + (1 - `___internal___`) * `___ps___` / (1 - `___ps___`)
    )

  structure(
    list(id_col=id_col, model = model,
         internal_df = filter(all_df_ps, `___internal___` == TRUE),
         external_df = filter(all_df_ps, `___internal___` == FALSE)
         ),
    class = c("prop_scr_obj")
  )
}

print.prop_scr_obj <- function(x){
  # cat the model
  # x$external_df |>
  #   select(id, ps, ipw)
}
