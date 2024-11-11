#' Analysis UI
#'
#' @param id mod id
#'
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
#' @param reactiveEndpoint reactive element with the type of endpoint
#'
#' @noMd
#' @importFrom shiny renderUI h3 reactive
analysisServer <- function(id, selections) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    bin_svr <- binaryServer("bin", selections)
    norm_svr <- normalServer("norm", selections)

    output$analysisUI <- renderUI({
      if (selections()$endPoint=="Binary") {
        binaryanalysisUI(ns("bin"))
      } else if (selections()$endPoint=="Normal") {
        normalanalysisUI(ns("norm"))
      } else {
        h3("bad")
      }
    })
    selected_inputs <- reactive({
      if (selections()$endPoint=="Binary") {
        bin_svr()
      } else if (selections()$endPoint=="Normal") {
        norm_svr()
      }
    })


    return(selected_inputs)

  })
}
