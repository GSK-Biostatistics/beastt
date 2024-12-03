
#' Inputs for plots
#'
#' @param id mod id
#'
#' @noRd
#' @noMd
#'
#' @importFrom shiny h4 hr fluidRow column checkboxGroupInput
plotUI <- function(id) {
  ns <- NS(id)
  tagList(
    h4("Plots"),
    hr(),
    fluidRow(
      column(4,
             shiny::checkboxGroupInput(ns("plotProp"),
                                       "Propensity Score Plot(s)",
                                       choices =
                                         c("Histogram", "Histogram - IPW",
                                           "Density", "Density - IPW",
                                           "Love"))),

      column(3,
             shiny::checkboxGroupInput(ns("plotPrior"),
                                       "Prior Plot",
                                       choices = c("Vague", "Power Prior"))),
      column(3,
             shiny::checkboxGroupInput(ns("plotPost"),
                                       "Posterior Plot",
                                       choices = c("Prior", "Posterior")))
    )
  )
}

#' Plot Inputs server
#'
#' @param id mod id
#' @param input_list list of inputs
#' @param robust robustify
#'
#' @noRd
#' @noMd
#' @importFrom shiny moduleServer moduleServer observeEvent updateCheckboxGroupInput
plotServer <- function(id, input_list, robust) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    observeEvent(robust, {
      choices <- if (is.null(robust)) {
        c("Vague", "Power Prior")
      } else if (robust) {
        c("Vague", "Power Prior", "Robust Mixture")
      } else {
        c("Vague", "Power Prior")
      }
      updateCheckboxGroupInput(session = session,
                               inputId = "plotPrior",
                               choices = choices
      )
    })

    plot_choices <- list(plotProp = input$plotProp,
           plotPrior = input$plotPrior,
           plotPost = input$plotPost)

    return(plot_choices)

  })
}

