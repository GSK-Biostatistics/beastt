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
                          "No borrowing"), selected="On control arm"),
    checkboxInput(ns("robustify"), "Robustify Power Prior", value=FALSE),
    uiOutput(ns("plots"))
  )
}

#' Binary Server UI
#'
#' @param id mod it
#' @noMd
#' @importFrom shiny renderUI reactive
#' @importFrom shinyjs hide show
binaryServer <- function(id, input_list) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    observeEvent(input$borrType, {
      if (input$borrType == "No borrowing") {
        hide("robustify")
      } else {
        show("robustify")
      }
    })

    output$plots <- renderUI({
      plotUI(ns("plot-select"))
    })

    plot_select <- reactiveVal(list())

    observeEvent(input$robustify, {
      plot_list <- plotServer("plot-select", base_input, input$robustify)
      plot_select(plot_list)
    })

    binary_selections <- reactive({list(
      borrType = input$borrType,
      robustify = input$robustify,
      plots = plot_select()
    )})

    return(binary_selections)
  })
}
