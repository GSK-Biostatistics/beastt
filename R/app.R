#' BDB Template Maker
#'
#' RStudio add-in to create template BDB code
#'
#' @export
#'
#'@import shiny
#'@import miniUI
#'@import bslib
#'@importFrom rstudioapi documentNew
bdb_template_maker <- function(){
  # Define UI for application
  ui <- fluidPage(
    id = "bdb_template_app",
    theme = bs_theme(bootswatch = "united",
                     primary = "#F36633",
                     secondary = "#f39633",
                     "navbar-bg" = "#F36633"
    ),
    titlePanel("BDB Template Maker"),
    sidebarLayout(
      sidebarPanel(
        h3("Inputs"),
        shiny::radioButtons("Purpose", "Purpose",
                            choices = c("Analysis", "Simulation")),
        h3("Study"),
        shiny::selectInput("endPoint", "Endpoint Type",
                           choices=c("Binary", "Normal", "Survival")),
        shiny::numericInput("ssIntArm", "Sample Size Internal Arm:", value = NULL),
        shiny::numericInput("RMPWeights", "Weights for RMP:", value = NULL),
        shiny::actionButton(inputId="submit", label= "Submit")
      ),
      mainPanel(
        h3("Data"),
        shiny::checkboxGroupInput(
          inputId = "dataChks",
          label = "",
          choices = c("External Data",
                      "Differing Sample Size per Arm",
                      "Estimate 'true' Covariate Effects",
                      "Underlying SD (Normal Case)",
                      "Robustify Power Prior")

        ),
        shiny::radioButtons("covDatGen",
                            label = "Covariate Data Generation",
                            choices = c("External Data Based (bootstrap)",
                                        "Assumed known covariance",
                                        "Other")),
        shiny::selectInput("borrType", "Type of Borrowing",
                           choices=c("On control arm",
                                     "On treatment arm",
                                     "No borrowing")),
        shiny::checkboxGroupInput(
          inputId = "Plots",
          label = "Output plots",
          choices=c("Histogram of Propensity Scores",
                    "Density of Inverse Probability Weights",
                    "Covariates Balance",
                    "Prior Distribution",
                    "Posterior Distributions")
        )
      )
    )
  )

  # Define server logic
  server <- function(input, output, session) {
    bslib::bs_theme()
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
