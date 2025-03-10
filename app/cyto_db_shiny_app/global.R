library(shiny)
library(reticulate)

# ----------------- Configure reticulate -----------------
use_virtualenv("/venv", required = TRUE)

# ----------------- Source Python script -----------------
source_python("parse_sdf.py")

# ----------------- Global: Parse SDF once -----------------
sdf_path <- "./CytoLabs_Database.sdf"
cat("Parsing SDF file at startup...\n")
mols_data <- parse_sdf_to_data(sdf_path)
cat("Finished parsing. Found", length(mols_data), "molecule(s).\n")

mols_data_r <- lapply(mols_data, function(x) {
  list(mol_idx = x$mol_idx, base64_png = x$base64_png, fields = as.list(x$fields))
})

df <- do.call(rbind, lapply(mols_data_r, function(x) {
  fieldText <- paste(unlist(x$fields), collapse = " ")
  data.frame(mol_idx = x$mol_idx, searchString = tolower(fieldText), stringsAsFactors = FALSE)
}))
