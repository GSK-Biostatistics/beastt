#' Create a Proprensity Score Object
#'
#' @param internal_df Internal dataset with one row per subject and all the
#'   variables needed to run the model
#' @param external_df External dataset with one row per subject and all the
#'   variables needed to run the model
#' @param model Model used to calculate propensity scores
#' @param ...
#'
#' @return `prop_scr_obj` object, with the internal and the external data and
#'   the propensity score and inverse probability weight calculated for each subject.
#' @export
#'
create_prop_scr <- function(internal_df, external_df,
               model, ...){

  i_df <- internal_df |> mutate(`___internal__` = TRUE)
  e_df <- external_df |> mutate(`___internal__` = FALSE)
  all_df <- bind_rows(i_df, e_df)


  all_df_ps <- all_df |>
    mutate(`__ps__` = glm(model, data = all_df, family = binomial)$fitted,
           ipw = `___internal__` + (1 - `___internal__`) * `__ps__` / (1 - `__ps__`)
    )

  structure(
    list(model = model,
         internal_df = filter(all_df_ps, `___internal__` == TRUE),
         external_df = filter(all_df_ps, `___internal__` == FALSE)
         ),
    class = c("prop_scr_obj")
  )
}
