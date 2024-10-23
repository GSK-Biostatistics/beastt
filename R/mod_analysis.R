analysisUI <- function(id) {
  ns <- NS(id)
  tagList(
    uiOutput("analysisUI")
  )
}

analysisServer <- function(id, reactiveEndpoint) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    observeEvent(reactiveEndpoint, {
      if (reactiveEndpoint()=="Binary") {
        # browser works here
        output$analysisUI <- shiny::renderUI({binaryanalysisUI("bin")})
        binaryServer("bin")
      } else if (reactiveEndpoint()=="Normal") {
        output$analysisUI <- renderUI(normalanalysisUI("norm"))
        normalServer("norm")
      } else {
        output$analysisUI <- renderUI({h3("bad")})
      }
    })
  })
}
