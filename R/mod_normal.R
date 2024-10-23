normalanalysisUI <- function(id) {
  ns <- NS(id)
  tagList(
    h4("Normal Analysis"),
    shiny::selectInput(ns("borrType"), "Type of Borrowing",
                       choices=c("On control arm",
                                 "On treatment arm",
                                 "No borrowing")),
    shiny::checkboxInput(ns("robustify"), "Robustify Power Prior"),
    shiny::radioButtons(ns("stddev"), "Standard Deviation",
                        choices = c("Known", "Unknown (Student's *t* approximation)"),
                        inline = TRUE),
  )
}

normalServer <- function(id) {
  moduleServer(id, function(input, output, session) {

    input$borrType

  })
}

