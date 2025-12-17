library(tidyverse)
library(httr)
# library(readr)
# library(stringr)
# library("xlsx")

# Generate URLs for the calculated number of mismatches and all condition
generate_urls <- function(sequence, mismatch_conditions, strands) {
  base_url <- "https://gggenome.dbcls.jp/hg38"
  conditions <- c("_RefSeqCurated_prespliced_d3g2202/",
                  "_RefSeqCurated_spliced_d3g2202/")
  
  mismatch_col <- character()
  strand_col <- character()
  condition_col <- character()
  url_col <- character()

  for (cond in conditions) {
    url <- paste0(base_url, cond, mismatch_conditions, strands, sequence, ".json")
    
    mismatch_col <- c(mismatch_col, mismatch_conditions)
    strand_col <- c(strand_col, strands)
    condition_col <- c(condition_col, cond)
    url_col <- c(url_col, url)
  }
  
  # save url, mismatch and condition to a dataframe
  urls_df <- data.frame(
    mismatch_condition = mismatch_col, 
    strand = strand_col,
    condition = condition_col,
    url = url_col,
    stringsAsFactors = FALSE
  )
  return(urls_df)
}

# Function to filter data based on the number of equal signs
filter_data <- function(data, criterion) {
  data_results <- data$results
  
  df_filtered <- data_results[
    grepl("^(NM_|chr|NR_)", data_results$name) &
      str_count(data_results$edit, "=") >= criterion, ]
  
  if (nrow(df_filtered) == 0) return(data.frame())
  
  n <- nrow(df_filtered)
  
  df <- data.frame(
    line = df_filtered$name,
    equal_signs = str_count(df_filtered$edit, "="),
    mismatches = df_filtered$mis,
    deletions = df_filtered$del,
    insertions = df_filtered$ins,
    subject_seq = df_filtered$sbjct,
    query_seq = df_filtered$query,
    start_target = df_filtered$position,
    end_target = df_filtered$position_end,
    snippet = df_filtered$snippet,
    snippet_start = df_filtered$snippet_pos,
    snippet_end = df_filtered$snippet_end,
    stringsAsFactors = FALSE
  )
  
  return(df)
}

# Function to fetch protein expression data from Protein Atlas
fetch_protein_expression <- function(protein_name) {
  base_url <- "https://www.proteinatlas.org/api/search_download.php?search="
  columns <- "&columns=g,gd,brain_RNA_amygdala,brain_RNA_basal_ganglia,brain_RNA_cerebellum,brain_RNA_cerebral_cortex,brain_RNA_choroid_plexus,brain_RNA_hippocampal_formation,brain_RNA_hypothalamus,brain_RNA_medulla_oblongata,brain_RNA_midbrain,brain_RNA_pons,brain_RNA_spinal_cord,brain_RNA_thalamus&compress=no&format=tsv"
  full_url <- paste0(base_url, protein_name, columns)
  response <- GET(full_url)
  
  if (status_code(response) == 200) {
    text_data <- content(response, "text", encoding = "UTF-8")
    df <- read_tsv(text_data)
    return(df)
  } else {
    print(paste("Failed to fetch data for protein", protein_name))
    return(NULL)
  }
}


all_offt <- function(sequence, mismatches_allowed) {
  # Calculate the number of mismatches based on the sequence length
  sequence_length = nchar(sequence)
  mismatch_conditions <- paste0(mismatches_allowed, "/")
  
  strands <- c("+/")
  
  urls_df <- generate_urls(sequence, mismatch_conditions, strands)
  
  # Initialize an empty data frame that will be filled with summary data
  summary_df <- data.frame()
  
  # Initialize an empty data frame that will be filled with gggenome data
  all_df <- data.frame(
    line = character(),
    equal_signs = integer(),
    mismatches = integer(),
    deletions = integer(), 
    insertions = integer(),
    subject_seq = character(),
    query_seq = character(),
    start_target = integer(),
    end_target = integer(),
    snippet = character(),
    snippet_start = integer(),
    snippet_end = integer(),
    stringsAsFactors = FALSE
  )
  
  
  # Loop over each URL and retrieve the data from it
  for (i in 1:nrow(urls_df)) {
    url <- urls_df$url[i]
    condition <- urls_df$condition[i]
    mismatch_condition <- urls_df$mismatch_condition[i]
    
    print(paste("Fetching data from:", url))
    response <- GET(url)
    
    if (status_code(response) == 200) {
      df_json <- fromJSON(
        content(response, as = "text", encoding = "UTF-8"),
        flatten = TRUE)
      filtered_data_df <- filter_data(df_json, mismatches_allowed)
      all_df <- bind_rows(all_df, filtered_data_df)
    } else {
      print("Failed to fetch data.")
      next
    }
    
    # if (nrow(all_df) == 0) {
    #   all_df <- data.frame(
    #     line = NULL,
    #     equal_signs = NULL,
    #     mismatches = NULL,
    #     deletions = NULL, 
    #     insertions = NULL,
    #     subject_seq = NULL,
    #     query_seq = NULL,
    #     start_target = NULL,
    #     end_target = NULL,
    #     snippet = NULL,
    #     snippet_start = NULL,
    #     snippet_end = NULL,
    #     stringsAsFactors = FALSE
    #   )
    # }
  }
  
  # Divide the data frame by condition
  spliced_df <- all_df[str_detect(str_to_lower(all_df$line), "\\|spliced"), ]
  prespliced_df <- all_df[str_detect(str_to_lower(all_df$line), "\\|prespliced"), ]
  other_df <- setdiff(all_df, rbind(spliced_df, prespliced_df))
  
  return(all_df)
}


# # ---------- Run ----------
# dataframes <- all_offt("TTTTTGCCATCCTGGGCGCT", 2)
# urls <- dataframes$urls
# summary_df <- dataframes$summary
# df <- dataframes$df
# spliced <- dataframes$spliced
# prespliced <- dataframes$prespliced
# other <- dataframes$other