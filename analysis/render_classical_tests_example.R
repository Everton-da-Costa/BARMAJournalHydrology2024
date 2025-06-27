# Render the R Markdown document to PDF
rmarkdown::render(
  input = "classical_tests_example.Rmd", # Source Rmd file
  output_format = "pdf_document",             # Output format: PDF
  output_dir = "../output",                   # Directory for generated PDF
  output_file = "classical_tests_example.pdf" # Name of output PDF file
)
