library(shiny)
library(shinyjs)
library(DT)
library(sf)
library(tidyverse)
library(bslib)
library(readxl)
library(ggthemes)
library(patchwork)

source("VisNet.R")

# Important functions:
# getPropData(data, rootID, depth): returns table of prop data by root node subgraph
# visVHFA(data, rootID, depth): returns visNet by root node subgraph.

ui <- fluidPage(
  shinyjs::useShinyjs(),
  title = NULL,
  theme = bslib::bs_theme(
    version = 5, bootswatch = "flatly",
    primary = "#378134"
  ),
  fluidRow(
    column(
      width = 6,
      fileInput(
        inputId = "data_input",
        label = "Excel file input (only .xls supported):",
        accept = c(".xls"),
        width = "100%"
      )
    ),
    column(
      width = 3,
      uiOutput("rootID_select")
    ),
    column(
      width = 3,
      selectInput(
        inputId = "depth_input",
        label = "Search depth:",
        choices = c(1, 2, 3),
        selected = 3,
        multiple = FALSE,
        selectize = FALSE,
        width = "100%"
      )
    )
  ),
  fluidRow(
    column(
      width = 10,
      visNetworkOutput(
        outputId = "plot_out"
      ) %>% 
        shinycssloaders::withSpinner(caption = "Loading", type = 3, color = "#378134", color.background = "white")
    ),
    column(
      width = 2,
      HTML(
        paste0(
          "<p>",
          "Selected ID shape: ●<br>",
          "Matching properties from the same dataset: ■<br>",
          "Matching properties from other datasets: ▲", "</p>",
          "<p>",
          "Selected ID and matching properties color: Green<br>",
          "Grand List color: Gold<br>",
          "HDS color: Orange<br>",
          "JFO color: Orchid<br>",
          "IA color: Slate blue<br>",
          "DoARH color: Red", "</p>"
        )
      )
    )
  ),
  fluidRow(
    column(
      width = 12,
      DTOutput(
        outputId = "table_out"
      ) %>%
        shinycssloaders::withSpinner(caption = "Loading", type = 3, color = "#378134", color.background = "white")
    )
  )
)

server <- (function(input, output, session) {
  shinyjs::hide("hide_this")
  dataPath <- reactive({
    shiny::req(input$data_input)
    inFile <- (input$data_input)
    inFile$datapath
  })
  data <- reactive({
    read_xls(dataPath(), sheet = 1)
  })
  propData <- reactive({
    getPropData(data(), dataPath(), input$rootID_input, input$depth_input) %>% 
      select(where(~!all(is.na(.x))))
  })
  visNet <- reactive({
    visVHFA(data(), dataPath(), input$rootID_input, input$depth_input)
  })
  output$rootID_select <- renderUI ({
    data_rootID <- data()
    shiny::req(data_rootID)
    if ("ID" %in% names(data_rootID)) {
      rootID_options <- sort(unique(data_rootID$ID))
      selectizeInput(
        inputId = "rootID_input",
        label = "Select an ID to view:",
        choices = rootID_options,
        selected = rootID_options[100]
      )
    } else {
      stop("No IDs found!")
    }
  })
  output$plot_out <- renderVisNetwork({visNet()})
  output$table_out <- renderDT({
    datatable(
      propData(),
      options = list(
        autoWidth = TRUE,
        scrollX = TRUE
      ),
      class = "compact cell-border"
    )
  })
})

# Create Shiny app ----
shinyApp(ui = ui, server = server)