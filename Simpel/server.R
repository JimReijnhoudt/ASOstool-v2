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
  
  # Reading in the csv data from a run to test the implementation of the RNaseH_script.R
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
  
  # A stored value for use in the second table. 
  selected_target <- reactiveVal(NULL)
  
  # This is the added code for the functionality of RNase H script.
  # The first observer onclick calls the RNaseH_script and return the values in a table on the next tabpanel.
  rnaseh_table_input <- function(row_number, data, table_id) {
    if (length(row_number) > 0) {
      
    row_data <- target_region_select[row_number, ]
    selected_target(row_data)
    
    output$rnaseh_title <- renderText({
      paste0("RNase H results for: ", row_data$name)
    })
    
    output$rnaseh_info <- renderText({
      HTML(paste0("length of sequence: ", row_data$length, "<br>",
             "Oligo sequence: ", row_data$oligo_seq))
    })
    
    output$rnaseh_results <- renderDataTable({
      datatable(rnaseh_results(row_data$name), selection = list(mode = 'single', selected = 1))
    })
    
    updateTabsetPanel(session, "tabs_main", selected = "RNase H cleavage results")
    
    proxy1 <- dataTableProxy("results1")
    selectRows(proxy1, NULL)
    }
  }
  
  observeEvent(input$results1_rows_selected, {
    rnaseh_table_input(input$results1_rows_selected, target_region_select, "results1")
  })
  
  observeEvent(input$results2_rows_selected, {
    rnaseh_table_input(input$results2_rows_selected, nucleobase_select, "results2")
  })
  
  # The second oberserver gives a visual of the cleavage site on the target sequence. 
  observeEvent(input$rnaseh_results_rows_selected, {
    row_number <- input$rnaseh_results_rows_selected
    if (length(row_number) == 0) return()
    
    row_data <- selected_target()
    if (is.null(row_data)) return()
    
      rnaseh_data <- rnaseh_results(row_data$name)
      selected_row <- rnaseh_data[row_number, ]
      
      position_string <- str_split_1(selected_row$position, " - ")
      start_pos <- as.numeric(position_string[1])
      end_pos <- as.numeric(position_string[2])
      
      target_seq <- row_data$name
      
      seq_visual <- paste0(
        "5' ",
        substr(target_seq, 1, start_pos - 1),
        "<span style='background-color: lightblue; color: red; font-weight: bold;'>",
        substr(target_seq, start_pos, start_pos + 6),
        "|",
        substr(target_seq, start_pos + 7, end_pos),
        "</span>",
        substr(target_seq, end_pos + 1, nchar(target_seq)),
        " 3'"
      )
      
      output$cleavage_visual <- renderUI({
        HTML(paste0("<div style='font-family: monospace; white-space: pre; font-size: 18px;'>", seq_visual, "</div>"))
      })
    })
}