getwd()
setwd("C:/Users/fayef/Documents/BI/BI3/periode_1/XEXT/ASOstool-v2")

# install.packages("tidyverse")
library(tidyverse)

# install.packages("httr")
library(httr)


generate_urls <- function(aso_sequence, mismatch_conditions) {
  base_url <- "https://gggenome.dbcls.jp/hg38"
  conditions <- c("_RefSeqCurated_prespliced_d3g2202/",
                "_RefSeqCurated_spliced_d3g2202/", ".p13_d3g2202/")
  
  mismatch_col <- character()
  condition_col <- character()
  url_col <- character()
  
  for (mc in mismatch_conditions) {
    for (cond in conditions) {
      url <- paste0(base_url, cond, mc, aso_sequence, ".txt")

      mismatch_col <- c(mismatch_col, mc)
      condition_col <- c(condition_col, cond)
      url_col <- c(url_col, url)
    }
  }
  
  urls_df <- data.frame(
    mismatch_condition = mismatch_col, 
    condition = condition_col,
    url = url_col,
    stringsAsFactors = FALSE
  )
  return(urls_df)
}

filter_data <- function(data, criterion) {
  start <- c('NM_', 'chr', 'NR_')
  lines <- str_split(data, "\n")[[1]]
  
  line_col <- character()
  equal_col <- integer()
  mismatch_col <- integer()
  
  for (line in lines) {
    if (any(str_starts(line, start))) {
      aantal_equal <- str_count(line, "=")
      if (aantal_equal >= criterion){
        line_col <- c(line_col, line)
        
        equal_col <- c(equal_col, aantal_equal)
        
        parts <- strsplit(line, "\t")[[1]]
        mm <- as.integer(parts[13])
        mismatch_col <- c(mismatch_col, mm)
      }
    }
  }
  
  df <- data.frame(
    line = line_col,
    equal_signs = equal_col,
    mismatches = mismatch_col,
    stringsAsFactors = FALSE
  )
  return(df)
}
# 
# fetch_protein_expression <- function(protein_name) {
#   base_url <- "https://www.proteinatlas.org/api/search_download.php?search="
#   columns <- "&columns=g,gd,brain_RNA_amygdala,brain_RNA_basal_ganglia,brain_RNA_cerebellum,brain_RNA_cerebral_cortex,brain_RNA_choroid_plexus,brain_RNA_hippocampal_formation,brain_RNA_hypothalamus,brain_RNA_medulla_oblongata,brain_RNA_midbrain,brain_RNA_pons,brain_RNA_spinal_cord,brain_RNA_thalamus&compress=no&format=tsv"
#   full_url <-paste0(base_url, protein_name, columns)
#   response <- GET(full_url)
#   
#   if (http_status(response) == 200) {
#     df <- content(r)
#   }
# }

main <- function() {
  # Define the file name at the top
  output_file_name <- "20250801_EH_UNC13A_SNP151&41_ASO146_R.xlsx"
  
  # Define the ASO sequence
  aso_sequence <- "CACACCATGCACATTCAA"  # Replace with your actual sequence
  
  # Calculate the number of mismatches based on the sequence length
  sequence_length = nchar(aso_sequence)
  mismatches_allowed = 2  # int(sequence_length * 0.1)  # 20% of the sequence length
  
  mismatch_conditions <- c("2/")
  urls_df <- generate_urls(aso_sequence, mismatch_conditions)
  
  summary_df <- data.frame(
    'Protein Hit' = character(), 
    'Source Sheets' = character(), 
    'Equal Signs' = integer(),
    'Total Mismatches' = integer(),
    'Gene description' = character(), 
    'Brain RNA - amygdala [nTPM]' = numeric(),
    'Brain RNA - basal ganglia [nTPM]' = numeric(),
    'Brain RNA - cerebellum [nTPM]' = numeric(),
    'Brain RNA - cerebral cortex [nTPM]' = numeric(),
    'Brain RNA - choroid plexus [nTPM]' = numeric(),
    'Brain RNA - hippocampal formation [nTPM]' = numeric(),
    'Brain RNA - hypothalamus [nTPM]' = numeric(),
    'Brain RNA - medulla oblongata [nTPM]' = numeric(),
    'Brain RNA - midbrain [nTPM]' = numeric(), 
    'Brain RNA - pons [nTPM]' = numeric(),
    'Brain RNA - spinal cord [nTPM]' = numeric(),
    'Brain RNA - thalamus [nTPM]' = numeric()
  )
  
  all_df <- data.frame(
    line = character(),
    equal_signs = integer(),
    mismatches = integer(),
    stringsAsFactors = FALSE
  )
  
  for (i in 1:nrow(urls_df)) { 
    mismatch_condition <- urls_df$mismatch_condition[i] 
    url <- urls_df$url[i] 
    condition <- urls_df$condition[i]
    
    print_statement <- paste("Fetching data from:", url)
    print(print_statement)
    
    filtered_results = list()
    response = GET(url)
    if (status_code(response) == 200) {
      print("Successfully fetched data.")
      text_data <- content(response, "text", encoding = "UTF-8")
      filtered_data_df <- filter_data(text_data, mismatches_allowed)
      
      all_df <- rbind(all_df, filtered_data_df)
    }
    else {
      print("Failed to fetch data.")
    }
  }
  return(list(urls = urls_df, summary = summary_df, df = all_df))
}
  
dataframes <- main()

urls <- dataframes$urls
summary_df <- dataframes$summary
df <- dataframes$df
