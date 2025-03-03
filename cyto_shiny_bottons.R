library(shiny)
library(reticulate)
library(dplyr)
library(base64enc)
library(png)
library(grid)

# ----------------- Configure reticulate -----------------
use_virtualenv("/home/mda25/myenv", required = TRUE)

# ----------------- Define Python helper code -----------------
py_run_string("
import base64
import io
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import rdDepictor
rdDepictor.SetPreferCoordGen(True)

def parse_sdf_to_data(sdf_path):
    suppl = Chem.SDMolSupplier(sdf_path, sanitize=True, removeHs=True, strictParsing=True)
    results = []
    mol_idx = 0
    
    for mol in suppl:
        if mol is None:
            continue
        Chem.rdDepictor.Compute2DCoords(mol)
        img = Draw.MolToImage(mol, size=(1000,1000))
        buf = io.BytesIO()
        img.save(buf, format='PNG')
        b64 = base64.b64encode(buf.getvalue()).decode('utf-8')
        
        fields_dict = {}
        for prop_name in mol.GetPropNames():
            fields_dict[prop_name] = mol.GetProp(prop_name)
        
        results.append({'mol_idx': mol_idx, 'base64_png': b64, 'fields': fields_dict})
        mol_idx += 1
    
    return results
")

# ----------------- Global: Parse SDF once -----------------
sdf_path <- "/home/mda25/cyto_db/CytoLabs_Database.sdf"
cat("Parsing SDF file at startup...\n")
mols_data <- py$parse_sdf_to_data(sdf_path)
cat("Finished parsing. Found", length(mols_data), "molecule(s).\n")

mols_data_r <- lapply(mols_data, function(x) {
  list(mol_idx = x$mol_idx, base64_png = x$base64_png, fields = as.list(x$fields))
})

df <- do.call(rbind, lapply(mols_data_r, function(x) {
  fieldText <- paste(unlist(x$fields), collapse = " ")
  data.frame(mol_idx = x$mol_idx, searchString = tolower(fieldText), stringsAsFactors = FALSE)
}))

# ----------------- Shiny UI -----------------
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

# ----------------- Shiny Server -----------------
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

# ----------------- Run App -----------------
shinyApp(ui, server)