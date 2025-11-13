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
library(openxlsx)

# Dont forget the source. 
source("../tools/RNaseH_script.R")

server <- function(input, output, session) {
  
  # Reading in the csv data from a run to test the implementation of the RNaseH_script.R
  target_region_select <- read_excel("../Files/SRGAP2A_filtered.xls")
  nucleobase_select <- read_excel("../Files/SRGAP2A_filtered.xls")
  
  # Render the tables.
  output$results1 <- renderDataTable({
    datatable(target_region_select, rownames = FALSE) %>%
      formatStyle(names(target_region_select), color = "black")
  })
  
  output$results2 <- renderDataTable({
    datatable(nucleobase_select, rownames = FALSE) %>%
      formatStyle(names(nucleobase_select), color = "black")
  })
  
  # A stored value for use in the second table and download. 
  selected_target <- reactiveVal(NULL)
  rnaseh_stored <- reactiveVal(NULL)
  
  # This is the added code for the functionality of RNase H script.
  # The first function calls the RNaseH_script and return the values in a table on the next tabpanel.
  rnaseh_table_input <- function(row_number, data, table_id) {
    if (length(row_number) > 0) {
      
    row_data <- data[row_number, ]
    selected_target(row_data)
    
    rnaseh_data <- rnaseh_results(
      selected_row_name = row_data$name,
      mod_5prime = 0,
      mod_3prime = 0
    )
    rnaseh_stored(rnaseh_data)
    
    output$rnaseh_title <- renderText({
      paste0("RNase H results for: ", row_data$name)
    })
    
    output$rnaseh_info <- renderText({
      HTML(paste0("length of sequence: ", row_data$length, "<br>",
             "Oligo sequence: ", row_data$oligo_seq))
    })
    
    output$rnaseh_results <- renderDataTable({
      datatable(rnaseh_data, selection = list(mode = 'single', selected = 1))
    })
    
    updateTabsetPanel(session, "tabs_main", selected = "RNase H cleavage results")
    
    proxy1 <- dataTableProxy("results1")
    selectRows(proxy1, NULL)
    }
  }
  
  # Download handler 
  output$download_rnaseh <- downloadHandler(
    filename = function() {
      row_data <- selected_target()
      if (is.null(row_data)){
        "RNaseH_results.xlsx"
      } else {
        paste0("RNaseH_results_", row_data$name, ".xlsx")
      }
    },
    content = function(file) {
      data <- rnaseh_stored()
      
      if (is.null(data)) data <- data.frame(Message = "No avalable data")
    
      data[] <- lapply(data, as.character)
      
      wb <- createWorkbook()
      addWorksheet(wb, "RNaseH_results")
      writeData(wb, "RNaseH_results", data)
      saveWorkbook(wb, file, overwrite = TRUE)
    }
  )
  
  # Table observers
  observeEvent(input$results1_rows_selected, {
    rnaseh_table_input(input$results1_rows_selected, target_region_select, "results1")
  })
  
  observeEvent(input$results2_rows_selected, {
    rnaseh_table_input(input$results2_rows_selected, nucleobase_select, "results2")
  })
  
  # Apply end modifications. 
  observeEvent(input$add_mods, {
    row_data <- selected_target()
    if (is.null(row_data)) return()
      
    rnaseh_data <- rnaseh_results(
      selected_row_name = row_data$name,
      mod_5prime = input$mod_5prime,
      mod_3prime = input$mod_3prime
    )
    rnaseh_stored(rnaseh_data)
    
    output$rnaseh_results <- renderDataTable({
      datatable(rnaseh_data, selection = list(mode = 'single', selected = 1))
    })
  })
  
  # The second observer object gives a visual of the cleavage site on the target sequence. 
  observeEvent(input$rnaseh_results_rows_selected, {
    row_number <- input$rnaseh_results_rows_selected
    if (length(row_number) == 0) return()
    
    row_data <- selected_target()
    if (is.null(row_data)) return()
    
    rnaseh_data <- rnaseh_stored()
    if (is.null(rnaseh_data)) return()
      
      selected_row <- rnaseh_data[row_number, ]
      
      position_string <- str_split_1(selected_row$position, " - ")
      start_pos <- as.numeric(position_string[1])
      end_pos <- as.numeric(position_string[2])
      
      target_seq <- row_data$name
      mod5 <- input$mod_5prime
      mod3 <- input$mod_3prime
      aso_len <- nchar(target_seq)
      
      mod5_region <- substr(target_seq, 1, mod5)
      mod3_region <- substr(target_seq, aso_len - mod3 + 1, aso_len)
      
      mid_region <- substr(target_seq, mod5 + 1, aso_len - mod3)
      
      cleavage_start <- start_pos
      cleavage_end <- start_pos + 6 
      
      seq_visual <- paste0(
        "<b style='color:darkorange;'>5'</b> ",
        "<span style='font-weight:bold; color:#d95f02;'>", mod5_region, "</span>",
        substr(target_seq, mod5 + 1, cleavage_start - 1),
        "<span style='background-color: lightblue; color: red; font-weight: bold;'>",
        substr(target_seq, cleavage_start, cleavage_end),
        "|",
        substr(target_seq, cleavage_end + 1, end_pos),
        "</span>",
        substr(target_seq, end_pos + 1, aso_len - mod3),
        "<span style='font-weight:bold; color:#1b9e77;'>", mod3_region, "</span>",
        " <b style='color:darkorange;'>3'</b>"
      )
      
      output$cleavage_visual <- renderUI({
        HTML(paste0("<div style='font-family: monospace; white-space: pre; font-size: 18px;'>", seq_visual, "</div>"))
      })
    })
}