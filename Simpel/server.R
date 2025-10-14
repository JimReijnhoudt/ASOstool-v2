# Imports
library(shinythemes)
library(shiny)
library(shinydashboard)
library(GenomicFeatures)
library(AnnotationDbi)
library(BSgenome.Hsapiens.NCBI.GRCh38)
library(biomaRt)
library(Biostrings)
library(tidyverse)
library(cluster)
library(rlang)
library(dplyr)
library(DT)
library(shinyBS)

# Dont forget the source. 
source("../tools/RNaseH_script.R")

server <- function(input, output, session) {
  
  target_region_select <- read_excel("../Files/SRGAP2A_filtered.xls")
  nucleobase_select <- read_excel("../Files/SRGAP2A_filtered.xls")
  
  output$results1 <- renderDataTable({
    datatable(target_region_select, rownames = FALSE) %>%
      formatStyle(names(target_region_select), color = "black")
  })
  
  output$results2 <- renderDataTable({
    datatable(nucleobase_select, rownames = FALSE) %>%
      formatStyle(names(nucleobase_select), color = "black")
  })
  
  # This is the added code for the functionality of RNase H script.
  observeEvent(input$results1_rows_selected, {
    row_number <- input$results1_rows_selected
    
    if (length(row_number) > 0) {
    row_data <- target_region_select[row_number, ]
    
    output$rnaseh_title <- renderText({
      paste0("RNase H results for: ", row_data$name)
    })
    
    output$rnaseh_info <- renderText({
      HTML(paste0("length of sequence: ", row_data$length, "<br>",
             "Oligo sequence: ", row_data$oligo_seq))
    })
    
    output$rnaseh_results <- renderDataTable({
      datatable(rnaseh_results(row_data$name))
    })
    
    updateTabsetPanel(session, "tabs_main", selected = "RNase H cleavage results")
    proxy1 <- dataTableProxy("results1")
    selectRows(proxy1, NULL)
    }
  })
}