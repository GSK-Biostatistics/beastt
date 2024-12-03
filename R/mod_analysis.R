#' Analysis UI
#'
#' @param id mod id
#'
#' @noRd
#' @noMd
#' @importFrom shiny NS tagList uiOutput
analysisUI <- function(id) {
  ns <- NS(id)
  tagList(
    uiOutput(ns("analysisUI"))
  )
}

#' Analysis Server
#'
#' @param id mod ID
#' @param input_list input from UI
#'
#' @noRd
#' @noMd
#' @importFrom shiny renderUI h3 reactive
analysisServer <- function(id, input_list) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    bin_svr <- binaryServer("bin", input_list)
    norm_svr <- normalServer("norm", input_list)

    observeEvent(input_list()$endPoint, {
      output$analysisUI <- renderUI({
        if (input_list()$endPoint=="Binary") {
          binaryanalysisUI(ns("bin"))
        } else if (input_list()$endPoint=="Normal") {
          normalanalysisUI(ns("norm"))
        } else {
          h3("bad")
        }
      })
    })

    selected_inputs <- reactive({
      if (input_list()$endPoint == "Binary") {
        bin_svr()
      } else if (input_list()$endPoint == "Normal") {
        norm_svr()
      }
    })


    return(selected_inputs)

  })
}
