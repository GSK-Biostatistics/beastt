analysisUI <- function(id) {
  ns <- NS(id)
  tagList(
    uiOutput(ns("analysisUI"))
  )
}

analysisServer <- function(id, reactiveEndpoint) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    bin_svr <- binaryServer("bin")
    norm_svr <- normalServer("norm")

    output$analysisUI <- shiny::renderUI({
      if (reactiveEndpoint()=="Binary") {
        binaryanalysisUI(ns("bin"))
      } else if (reactiveEndpoint()=="Normal") {
        normalanalysisUI(ns("norm"))
      } else {
        h3("bad")
      }
    })
    selected_inputs <- reactive({
      if (reactiveEndpoint()=="Binary") {
        bin_svr()
      } else if (reactiveEndpoint()=="Normal") {
        norm_svr()
      }
    })


    return(selected_inputs)

  })
}
