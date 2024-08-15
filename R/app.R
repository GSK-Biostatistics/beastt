#' BDB Simulator Template Maker
#'
#' RStudio add-in to create template BDB simulation code
#'
#' @export
#'
#'@import shiny
#'@import miniUI
#'@import bslib
#'@importFrom rstudioapi documentNew
bdb_simulator <- function(){
  # Define UI for application
  ui <- page_sidebar(
    id = "bdb_sim_app",
    theme = bs_theme(bootswatch = "united",
                     primary = "#F36633",
                     secondary = "#f39633",
                     "navbar-bgkkjj" = "#F36633"
    ),
    title ="BDB Simulator Template Maker",
    sidebar = sidebar(
      h3("Simulation"),
      shiny::numericInput("seed", "Seed:", value = 1234),
      shiny::radioButtons("simType", "Simulation Structure",
                          choices = c("By scenario", "By iteration")),
      h3("Study"),
      shiny::numericInput("ssIntArm", "Sample Size Internal Arm:", value = NULL),
      shiny::numericInput("RMPWeights", "Weights for RMP:", value = NULL),
      shiny::actionButton(inputId="submit", label= "Submit")
    ),
    navset_card_pill(
      placement = "above",
      nav_panel(title = "Data",
                shiny::checkboxGroupInput(
                  inputId = "dataChks",
                  label = "",
                  choices = c("External Data",
                              "Differing Sample Size per Arm",
                            "Estimate 'true' Covariate Effects",
                            "Underlying SD (Normal Case)")

                ),
                shiny::radioButtons("covDatGen",
                                    label = "Coraviate Data Generation",
                                    choices = c("External Data Based (bootstrap)",
                                                "Assumed known covariance",
                                                "Other"))
                ),
      nav_panel(title = "Priors", p("Second tab content.")),
      nav_panel(title = "Posteriors and Inference", p("Third tab content"))
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
