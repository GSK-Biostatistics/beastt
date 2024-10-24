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
                        width = '100%',
                        choices = c("Known", "Unknown (Student's t approximation)")
                        ),
  )
}

normalServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    return(reactive({list(
      borrType = input$borrType,
      robustify = input$robustify,
      stddev = input$stddev

    )}))

  })
}

