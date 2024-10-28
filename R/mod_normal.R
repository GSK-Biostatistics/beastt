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
                                 "No borrowing")),
    checkboxInput(ns("robustify"), "Robustify Power Prior"),
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
#'
#' @noMd
#' @importFrom shiny renderUI reactive
normalServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    plot_select <- plotServer("plot-select")

    output$plots <- renderUI({
      plotUI(ns("plot-select"), input$robustify)
    })


    return(reactive({list(
      borrType = input$borrType,
      robustify = input$robustify,
      stddev = input$stddev,
      plots = plot_select()

    )}))

  })
}

