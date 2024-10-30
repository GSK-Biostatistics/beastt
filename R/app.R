#' BDB Code Template Maker
#'
#' RStudio add-in to create template BDB code
#'
#' @export
#'
#'@import shiny
#'@import miniUI
#'@importFrom shinyjs useShinyjs
#'@importFrom bslib bs_theme page_sidebar sidebar
#'@importFrom rstudioapi documentNew
bdb_code_template_maker <- function(){
  # Define UI for application
  ui <- page_sidebar(
    shinyjs::useShinyjs(),
    id = "bdb_template_app",
    theme = bs_theme(bootswatch = "united",
                     primary = "#F36633",
                     secondary = "#f39633",
                     "navbar-bg" = "#F36633"
    ),
    title ="BDB Code Template Maker",
    sidebar=sidebar(
      h3("Study Design"),
      radioButtons("purpose", "Purpose",
                          choices = c("Analysis", "Simulation")),
      selectInput("endPoint", "Endpoint Type",
                         choices=c("Binary", "Normal", "Survival")),
      actionButton(inputId="submit", label= "Submit")
    ),
    analysisUI("analysis")
  )

  # Define server logic
  server <- function(input, output, session) {
    bs_theme()
    reactiveEndpoint <- reactive(input$endPoint)
    all_inputs <- analysisServer("analysis", reactiveEndpoint)
    observeEvent(input$submit,{
      write_code(input$purpose, input$endPoint, all_inputs())
      stopApp()

    })

  }
  # # Run the application
  runGadget(ui, server
              , viewer =dialogViewer("", width = 1000, height = 800),
            )
}
