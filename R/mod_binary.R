binaryanalysisUI <- function(id) {
  ns <- NS(id)
  tagList(
    h4("Binary Analysis"),
    shiny::selectInput(ns("borrType"), "Type of Borrowing",
                       choices=c("On control arm",
                                 "On treatment arm",
                                 "No borrowing")),
    shiny::checkboxInput(ns("robustify"), "Robustify Power Prior")
  )
}

binaryServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    reactive({list(borrType = input$borrType,
                   robustify = input$robustify)})
  })
}

