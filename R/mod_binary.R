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
#' @param id mod id
#' @param input_list input from UI
#' @noMd
#' @importFrom shiny renderUI reactive reactiveVal observeEvent
binaryServer <- function(id, input_list) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    observeEvent(input$borrType, {
      if (input$borrType == "No borrowing") {
        shinyjs::hide("robustify")
      } else {
        shinyjs::show("robustify")
      }
    })

    plot_select <- reactiveVal(list())

    observeEvent(input_list()$purpose, {
      if (input_list()$purpose == "Analysis") {
        output$plots <- renderUI({
          plotUI(ns("plot-select"))
        })

        observeEvent(input$robustify, {
          plot_list <- plotServer("plot-select", input_list, input$robustify)
          plot_select(plot_list)
        })
      }
    })

    binary_selections <- reactive({list(
      borrType = input$borrType,
      robustify = input$robustify,
      plots = plot_select()
    )})

    return(binary_selections)
  })
}
