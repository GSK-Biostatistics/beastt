#' norm inputs UI
#'
#' @param id mod id
#'
#' @noMd
#' @importFrom shiny NS tagList h4 hr selectInput checkboxInput radioButtons uiOutput
normalanalysisUI <- function(id) {
  ns <- NS(id)
  tagList(
    h4("Normal Analysis"),
    hr(),
    selectInput(ns("borrType"), "Type of Borrowing",
                choices=c("On control arm",
                          "On treatment arm",
                          "No borrowing"), selected="On control arm"),
    checkboxInput(ns("robustify"), "Robustify Power Prior", value=FALSE),
    radioButtons(ns("stddev"), "Standard Deviation",
                 width = '100%',
                 choices = c("Known", "Unknown (Student's t approximation)")
    ),
    uiOutput(ns("plots"))
  )
}

#' norm input server
#'
#' @param id mod id
#' @param input_list input from UI
#'
#' @noMd
#' @importFrom shiny renderUI reactive observeEvent reactiveVal
#' @importFrom shinyjs hide show
normalServer <- function(id, input_list) {
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
      plot_list <- plotServer("plot-select", input_list, input$robustify)
      plot_select(plot_list)
    })

    normal_selections <- reactive({list(
      borrType = input$borrType,
      robustify = input$robustify,
      stddev = input$stddev
    )})

    return(normal_selections)
  })
}

