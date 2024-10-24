#' BDB Code Template Maker
#'
#' RStudio add-in to create template BDB code
#'
#' @export
#'
#'@import shiny
#'@import miniUI
#'@import bslib
#'@importFrom rstudioapi documentNew
bdb_code_template_maker <- function(){
  # Define UI for application
  ui <- page_sidebar(
    id = "bdb_template_app",
    theme = bs_theme(bootswatch = "united",
                     primary = "#F36633",
                     secondary = "#f39633",
                     "navbar-bg" = "#F36633"
    ),
    title ="BDB Code Template Maker",
    sidebar=sidebar(
      h3("Study Design"),
      shiny::radioButtons("Purpose", "Purpose",
                          choices = c("Analysis", "Simulation")),
      shiny::selectInput("endPoint", "Endpoint Type",
                         choices=c("Binary", "Normal", "Survival")),
      shiny::actionButton(inputId="submit", label= "Submit")
    ),
    analysisUI("analysis")
  )

  # Define server logic
  server <- function(input, output, session) {
    bslib::bs_theme()
    reactiveEndpoint <- reactive(input$endPoint)
    test <- analysisServer("analysis", reactiveEndpoint)
    observeEvent(test(),{
      print(test())

    })
    observeEvent(input$submit,{
    })

  }
  # # Run the application
  runGadget(ui, server
              # , viewer =dialogViewer("", width = 1000, height = 800),
            )
}
