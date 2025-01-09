#' BDB Code Template Maker
#'
#' RStudio add-in to create template BDB code
#'
#' @export
#'
#'@import shiny
#'@import miniUI
#'@importFrom shinyjs useShinyjs
#'@importFrom bslib bs_theme page_sidebar sidebar navset_card_pill nav_panel nav_show nav_hide nav_select
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
    tags$style(HTML("
      hr {
          margin: var(--bs-space-0) !important;
          padding: var(--bs-space-0) !important;
      }
      .bslib-card .card-body {
          overflow: hidden;
      }
    ")),
    title ="BDB Code Template Maker",
    sidebar=sidebar(
      h3("Study Design"),
      radioButtons("purpose", "Purpose",
                   choices = c("Simulation", "Analysis")),
      selectInput("endPoint", "Endpoint Type",
                  choices=c("Binary", "Normal", "Time to Event")),
      actionButton(inputId="submit", label= "Submit")
    ),
    navset_card_pill(
      placement = "above",
      id = "tabs",
      nav_panel(title = "Simulation",
                simulationUI("simulation")),
      nav_panel(title = "Analysis",
                analysisUI("analysis"))
    )
  )

  # Define server logic
  server <- function(input, output, session) {
    bs_theme()

    observeEvent(input$purpose, {
      if (input$purpose=="Simulation") {
        nav_show("tabs", "Simulation", select=TRUE, session=session)
      } else {
        nav_hide("tabs", "Simulation", session=session)
        nav_select("tabs", "Analysis", session=session)
      }
    })

    base_input <- reactive({
      list(purpose=input$purpose, endPoint=input$endPoint)
    })

    simulation_out <- simulationServer("simulation", base_input)
    analysis_out <- analysisServer("analysis", base_input)

    observeEvent(input$submit, {
      req(simulation_out)
      req(analysis_out)

      final_inputs <- list(simulation_inputs = simulation_out(),
                           analysis_inputs = analysis_out())

      write_code(input$purpose, input$endPoint, final_inputs)
      stopApp()

    })

  }
  # # Run the application
  runGadget(ui, server,
            viewer = dialogViewer("", width = 1000, height = 800),
  )
}
