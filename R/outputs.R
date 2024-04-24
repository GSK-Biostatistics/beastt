#' Create a Pre-and-Post-Weighting Summary Table
#'
#' @param df Combined dataset with one row per subject and all the
#'   variables (subgroups) to be displayed in the summary table
#' @param vars Variables (subgroups) to be displayed in the summary table
#' @param internal Indicator variable identifying internal vs external subjects
#' @param trt Treatment variable, should be set to "Control" if borrowing from
#'   external controls
#' @param internal_df Internal dataset with one row per subject and all the
#'   variables needed to run the model
#' @param external_df External dataset with one row per subject and all the
#'   variables needed to run the model
#' @param id_col Name of the column in both datasets used to identify each
#'   subject. It must be the same across datasets
#' @param model Model used to calculate propensity scores
#' @param ... Optional arguments
#'
#' @return `tibble` object, with the adjusted and unadjusted
#' @export
#' @importFrom rlang f_rhs f_lhs `f_lhs<-` is_formula enquo as_name
#' @importFrom stringr str_split_1 str_detect str_split_i
#' @importFrom cli cli_abort
#' @importFrom purrr reduce discard
#' @importFrom rlang sym
#' @importFrom dplyr mutate filter tibble as_tibble
#' @importFrom stats glm
#' @importFrom cobalt bal.tab
#' @examples
#' require(extraDistr)
#' df <- data.frame(
#'   id_col = 1:40,
#'   cov1 = rep(c("a", "b", "c", "d"), 10),
#'   cov2 = rep(c("e", "f", "g", "h"), each = 10),
#'   internal = rdunif(40, 0, 1),
#'   trt = rep(c("Active", "Control"), each = 20)
#' )
#' model <- as.formula("~cov1 + cov2")
#' internal_df <- df %>%
#' filter(internal == 1)
#' external_df <- df %>%
#' filter(internal == 0)
#' pre_post_summary(df, vars = c("cov1", "cov2"), df$internal, trt = df$trt,
#' internal_df, external_df, model = model, id_col = id_col)
pre_post_summary <- function(df, vars, internal, trt, internal_df,
                                external_df, id_col, model, ...){

  freqs <- list()

  for(i in 1:length(vars)){

    var_name <- vars[i]

    # Calculate frequencies
    ftable <- table(df[[var_name]], internal, trt)
    freqdf <- as.data.frame(ftable)

    # Calculate percentages
    ptable <- prop.table(ftable, margin = 1) * 100
    pdf <- as.data.frame(ptable)
    names(pdf)[names(pdf) == "Freq"] <- "Perc"

    # Tidy up summaries
    fptable <- merge(freqdf, pdf)
    fptable$variable <- var_name
    names(fptable)[names(fptable) == "Var1"] <- "value"
    freqs[[var_name]] <- fptable

  }

  # Combine all summaries to a single data.frame
  freqs_out <- do.call(rbind, freqs)
  freqs_out_fmt <- freqs_out |>
    select(variable, value, internal, trt, Freq, Perc)
  rownames(freqs_out) <- 1:nrow(freqs_out)

  # Calculate propensity scores
  pscore <- calc_prop_scr(internal_df, external_df,
                          id_col, model)$abs_std_mean_diff
  pscore$variable <- str_split_i(pscore$covariate, "_", 1)
  pscore$value <- str_split_i(pscore$covariate, "_", -1)

  # Combine outputs
  out <- merge(freqs_out_fmt, pscore) |>
    select(-covariate) |>
    as_tibble()

  names(out) <- tolower(names(out))

  return(out)
}
