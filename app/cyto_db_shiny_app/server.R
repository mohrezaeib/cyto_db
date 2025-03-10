# server.R
library(shiny)
library(base64enc)
library(dplyr)
library(png)
library(grid)
server <- function(input, output, session) {
  rv <- reactiveValues(matchedIdx = df$mol_idx, selected = 1)
  
  observeEvent(input$searchBtn, {
    term <- tolower(input$searchTerm)
    if (term == "") {
      rv$matchedIdx <- df$mol_idx
    } else {
      matched <- df$mol_idx[grepl(term, df$searchString, fixed = TRUE)]
      rv$matchedIdx <- matched
    }
    rv$selected <- 1
  })
  
  observeEvent(input$prevBtn, {
    if (rv$selected > 1) {
      rv$selected <- rv$selected - 1
    }
  })
  
  observeEvent(input$nextBtn, {
    if (rv$selected < length(rv$matchedIdx)) {
      rv$selected <- rv$selected + 1
    }
  })
  
  output$searchCount <- renderText({
    paste("Found", length(rv$matchedIdx), "matching molecules.")
  })
  
  output$molPlot <- renderPlot({
    if (length(rv$matchedIdx) == 0) {
      plot.new()
      text(0.5, 0.5, "No matches found.")
      return()
    }
    actualIdx <- rv$matchedIdx[rv$selected]
    molData <- mols_data_r[[actualIdx + 1]]
    b <- base64decode(molData$base64_png)
    tf <- tempfile(fileext = ".png")
    writeBin(b, tf)
    img <- readPNG(tf)
    grid.raster(img)
  })
  
  output$molFields <- renderTable({
    if (length(rv$matchedIdx) == 0) return(NULL)
    actualIdx <- rv$matchedIdx[rv$selected]
    molData <- mols_data_r[[actualIdx + 1]]
    data.frame(Field = names(molData$fields), Value = unlist(molData$fields), stringsAsFactors = FALSE)
  })
}