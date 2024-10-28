binaryanalysisUI <- function(id) {
  ns <- NS(id)
  tagList(
    h4("Binary Analysis"),
    hr(),
    shiny::selectInput(ns("borrType"), "Type of Borrowing",
                       choices=c("On control arm",
                                 "On treatment arm",
                                 "No borrowing")),
    shiny::checkboxInput(ns("robustify"), "Robustify Power Prior"),
    uiOutput(ns("plots"))
  )
}

binaryServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    plot_select <- plotServer("plot-select")

    output$plots <- renderUI({
      plotUI(ns("plot-select"), input$robustify)
    })

    reactive({list(borrType = input$borrType,
                   robustify = input$robustify,
                   plots = plot_select())})
  })

}

