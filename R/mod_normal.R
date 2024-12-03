#' norm inputs UI
#'
#' @param id mod id
#'
#' @noRd
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
#' @noRd
#' @noMd
#' @importFrom shiny renderUI reactive observeEvent reactiveVal
normalServer <- function(id, input_list) {
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

    normal_selections <- reactive({list(
      borrType = input$borrType,
      robustify = input$robustify,
      stddev = input$stddev
    )})

    return(normal_selections)
  })
}

