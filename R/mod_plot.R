plotUI <- function(id, robust) {
  ns <- NS(id)
  if(robust){
    prior_choices <- c("Vague",
                       "Power Prior",
                       "Robust Mixture")
  } else {
    prior_choices <- c("Vague",
                       "Power Prior")
  }
  tagList(
    h4("Plots"),
    hr(),
    fluidRow(
      column(4,
             shiny::checkboxGroupInput(ns("plotProp"),
                                       "Propensity Score Plot(s)",
                                       choices =
                                         c("Historgram", "Historgram - IPW",
                                           "Density", "Density - IPW",
                                           "Love"))),

      column(3,
             shiny::checkboxGroupInput(ns("plotPrior"),
                                       "Prior Plot",
                                       choices = prior_choices)),
      column(3,
             shiny::checkboxGroupInput(ns("plotPost"),
                                       "Posterior Plot",
                                       choices = c("Prior", "Posterior")))
    )
  )
}

plotServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    return(reactive({list(
      plotProp = input$plotProp,
      plotPrior = input$plotPrior,
      plotPost = input$plotPost
    )}))

  })
}

