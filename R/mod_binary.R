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
binaryServer <- function(id, selections) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # want to take the reactive list selections (created as all_inputs in app.R),
    # add the inputs borrType and robustify, and create a new reactive list

    observeEvent(input$borrType, {
      toUpdate = selections()
      toUpdate$borrType <- input$borrType
      selections = reactive(toUpdate)
    })

    observeEvent(input$robustify, {
      toUpdate = selections()
      toUpdate$robustify <- input$robustify
      selections = reactive(toUpdate)
    })

    #selections <- reactiveValues(items=list())
    struct <- reactiveVal(NULL)
    # observeEvent(selections(),{
    #   struct(selections())
    # })
    #
    # observeEvent(input$borrType,{
    #   struct(input$borrType)
    # })
    #
    # observeEvent(input$robustify,{
    #   struct(input$robustify)
    # })

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

    plot_select <- plotServer("plot-select", selections)

    return(selections)

    # reactive({list(borrType = input$borrType,
    #                robustify = input$robustify,
    #                plots = plot_select()
    # )})
  })
}
