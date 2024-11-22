#' Simulation UI
#'
#' @param id mod ID
#'
#' @noMd
#' @importFrom shiny NS tagList h4 radioButtons checkboxGroupInput checkboxInput
simulationUI <- function(id) {
  ns <- NS(id)
  tagList(
    h4("Simulation"),
    radioButtons(ns("internal"), "Internal Data",
                 choices = c("Bootstrap", "Simulate from scratch"),
                 selected = "Bootstrap"),
    radioButtons(ns("parallel"), "Parallelisation",
                 choices = c("by iteration", "by scenario"),
                 selected = "by iteration"),
    checkboxGroupInput(ns("scenarioOptions"), "",
                       choices=c("Vary internal control sample size",
                                 "Weight of informative component of RMP",
                                 "Treatment effect", "Drift")),
    checkboxInput(ns("covImb"), "Covariate Imbalance", value = FALSE),
    checkboxGroupInput(ns("opChar"), "Operating Characteristics",
                       choices = c("Iteration level results"))
  )
}

#' Simulation Server
#'
#' @param id mod ID
#'
#' @noMd
#'
#' @importFrom shiny renderUI h3 reactive observeEvent
#' @importFrom shinyjs hide show
simulationServer <- function(id, input_list) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    observeEvent(input$parallel, {
      if (input$parallel == "by iteration") {
        hide("scenarioOptions")
      } else {
        show("scenarioOptions")
      }
    })
    observeEvent(input$internal, {
      if (input$internal == "Bootstrap") {
        show("covImb")
      } else {
        hide("covImb")
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
