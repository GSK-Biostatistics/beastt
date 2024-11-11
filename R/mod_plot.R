
#' Inputs for plots
#'
#' @param id mod id
#' @param robust Boolean if the
#'
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
#'
#' @noMd
#' @importFrom shiny moduleServer moduleServer
plotServer <- function(id, selections) {
  moduleServer(id, function(input, output, session) {
    observeEvent(selections()$robustify,{
      choices <- if (is.null(selections()$robustify)) {
        c("Vague", "Power Prior")
      } else if (selections()$robustify) {
        c("Vague", "Power Prior", "Robust Mixture")
      } else {
        c("Vague", "Power Prior")
      }
      updateCheckboxGroupInput(session = session,
                               inputId = "plotPrior",
                               choices = choices
      )
    })

    struct <- reactiveVal(NULL)
    observeEvent(selections(),{
      struct(selections())
    })

    observeEvent(input$plotProp,{
      struct(input$plotProp)
    })

    observeEvent(input$plotPrior,{
      struct(input$plotPrior)
    })

    observeEvent(input$plotPost,{
      struct(input$plotPost)
    })

    # old_selections <- reactiveVal(selections())
    # observeEvent(input$plotProp, {
    #   temp <- old_selections()
    #   browser()
    #   temp$plotProp <- input$plotProp
    #   old_selections(temp)
    # })

    return(struct)


    # return(reactive({list(
    #   plotProp = input$plotProp,
    #   plotPrior = input$plotPrior,
    #   plotPost = input$plotPost
    # )}))

  })
}

