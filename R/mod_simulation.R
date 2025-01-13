#' Simulation UI
#'
#' @param id mod ID
#'
#' @noRd
#' @noMd
#' @importFrom shiny NS tagList h4 radioButtons checkboxGroupInput checkboxInput
#' @importFrom bslib tooltip
simulationUI <- function(id) {
  ns <- NS(id)
  tagList(
    h4("Simulation"),
    hr(),
    fluidRow(
      column(6,
             radioButtons(ns("internal"), strong("Internal Covariate Data"),
                          choices = c("Simulate from scratch", "Bootstrap"),
                          selected = "Simulate from scratch"),
             checkboxInput(ns("covImb"), strong("Covariate Imbalance"), value = FALSE)),
      column(6,
             radioButtons(ns("parallel"), strong("Parallelisation"),
                          choices = c("by scenario", "by iteration"),
                          selected = "by scenario"))
    ),
    tooltip(
      checkboxGroupInput(ns("scenarioOptions"), strong("Options to vary"),
                         choices=c("Internal control sample size"="samplesize",
                                   "Weight of informative component of RMP"="weight",
                                   "Treatment effect"="trteffect", "Drift"="drift")),
      "Information about options"
    ),
    checkboxGroupInput(ns("opChar"), strong("Operating Characteristics"),
                       choices = c("Iteration level results"))
  )
}

#' Simulation Server
#'
#' @param id mod ID
#' @param input_list input from UI
#'
#' @noRd
#' @noMd
#'
#' @importFrom shiny renderUI h3 reactive observeEvent
simulationServer <- function(id, input_list) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    observeEvent(input$parallel, {
      if (input$parallel == "by iteration") {
        shinyjs::hide("scenarioOptions")
      } else {
        shinyjs::show("scenarioOptions")
      }
    })
    observeEvent(input$internal, {
      if (input$internal == "Bootstrap") {
        shinyjs::show("covImb")
      } else {
        shinyjs::hide("covImb")
      }
    })

    simulation_selections <- reactive({
      list(internal=input$internal,
           parallel=input$parallel,
           scenarioOptions=input$scenarioOptions,
           covImb=input$covImb,
           opChar=input$opChar)
    })

    return(simulation_selections)
  })
}
