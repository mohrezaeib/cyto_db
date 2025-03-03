# ui.R
library(shiny)

ui <- fluidPage(
  titlePanel("Single-Molecule Viewer with RDKit + reticulate"),
  fluidRow(
    column(4, textInput("searchTerm", "Search in all SDF fields:", value = "")),
    column(2, actionButton("searchBtn", "Search"))
  ),
  fluidRow(
    column(
      12,
      wellPanel(
        h4("Search Results:"),
        textOutput("searchCount"),
        actionButton("prevBtn", "Previous"),
        actionButton("nextBtn", "Next")
      )
    )
  ),
  fluidRow(
    column(6, plotOutput("molPlot")),
    column(6, tableOutput("molFields"))
  )
)
