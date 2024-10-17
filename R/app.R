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
      shiny::selectInput("borrType", "Type of Borrowing",
                         choices=c("On control arm",
                                   "On treatment arm",
                                   "No borrowing")),
      shiny::checkboxInput("robustify", "Robustify Power Prior"),
      shiny::actionButton(inputId="submit", label= "Submit")
    ),
    card(
      h3("Outputs"),
      shiny::checkboxGroupInput(
        inputId = "Plots",
        label = "",
        choices=NULL
      )
    )
  )

  # Define server logic
  server <- function(input, output, session) {
    bslib::bs_theme()
    shiny::observeEvent(input$endPoint, {
      endpoint <- input$endPoint
      if (endpoint=="Binary"){
        choices <- c("Histogram of Propensity Scores",
                     "Density of Inverse Probability Weights",
                     "Covariates Balance",
                     "Prior Distribution",
                     "Posterior Distributions")
        shiny::updateCheckboxGroupInput(session,
                                        inputId="Plots",
                                        label = "Plots",
                                        choices = choices)
      }
      else if (endpoint=="Normal"){
        choices <- character(0)
        shiny::updateCheckboxGroupInput(session,
                                        inputId="Plots",
                                        label = "",
                                        choices = choices)
      }
      else if (endpoint=="Survival"){
        choices <- character(0)
        shiny::updateCheckboxGroupInput(session,
                                        inputId="Plots",
                                        label = "",
                                        choices = choices)
      }
    })
    shiny::observeEvent(input$submit, {
      rstudioapi::documentNew(
        "############## Simulation Template ################",
        type = "r"
      )
      invisible(stopApp())
    })
  }
  # Run the application
  runGadget(ui, server, viewer =
              dialogViewer("", width = 1000, height = 800))
}
