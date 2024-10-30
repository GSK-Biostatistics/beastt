#' Binary Input UI
#'
#' @param id mod ID
#'
#' @noMd
#' @importFrom shiny NS tagList h4 hr selectInput checkboxInput uiOutput
binaryanalysisUI <- function(id) {
  ns <- NS(id)
  tagList(
    h4("Binary Analysis"),
    hr(),
    selectInput(ns("borrType"), "Type of Borrowing",
                choices=c("On control arm",
                          "On treatment arm",
                          "No borrowing")),
    checkboxInput(ns("robustify"), "Robustify Power Prior"),
    uiOutput(ns("plots"))
  )
}

#' Binary Server UI
#'
#' @param id mod it
#' @noMd
#' @importFrom shiny renderUI reactive
#' @importFrom shinyjs hide show
binaryServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    observeEvent(input$borrType, {
      if (input$borrType == "No borrowing") {
        hide("robustify")
      } else {
        show("robustify")
      }
    })

    plot_select <- plotServer("plot-select")

    output$plots <- renderUI({
      plotUI(ns("plot-select"), input$robustify)
    })

    reactive({list(borrType = input$borrType,
                   robustify = input$robustify,
                   plots = plot_select()
    )})
  })
}
