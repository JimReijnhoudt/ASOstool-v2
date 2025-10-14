library(tidyverse)
library(httr)
# library(readr)
# library(stringr)
library("xlsx")

# Generate URLs for the calculated number of mismatches and all condition
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
  
  # save url, mismatch and condition to a dataframe
  urls_df <- data.frame(
    mismatch_condition = mismatch_col, 
    condition = condition_col,
    url = url_col,
    stringsAsFactors = FALSE
  )
  return(urls_df)
}

# Function to filter data based on the number of equal signs
filter_data <- function(data, criterion) {
  start <- c('NM_', 'chr', 'NR_')
  lines <- str_split(data, "\n")[[1]]
  
  line_col <- character()
  equal_col <- integer()
  mismatch_col <- integer()
  ins_col <- integer()
  del_col <- integer()
  
  for (line in lines) {
    if (any(str_starts(line, start))) {
      aantal_equal <- str_count(line, "=")
      if (aantal_equal >= criterion){
        line_col <- c(line_col, line)
        
        equal_col <- c(equal_col, aantal_equal)
        
        parts <- strsplit(line, "\t")[[1]]
        mm <- as.integer(parts[13])
        # mm <- sum(as.integer(tail(parts, 2)))
        mismatch_col <- c(mismatch_col, mm)
        
        del <- as.integer(parts[14])
        ins <- as.integer(parts[15])
        
        del_col <- c(del_col, del)
        ins_col <- c(ins_col, ins)
      }
    }
  }
  
  df <- data.frame(
    line = line_col,
    equal_signs = equal_col,
    mismatches = mismatch_col,
    deletions = del_col, 
    insertions = ins_col,
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


main <- function() {
  # Define the file name at the top
  output_file_name <- "20250801_EH_UNC13A_SNP151&41_ASO146.xlsx"

  # Define the ASO sequence
  aso_sequence <- "CACACCATGCACATTCAA"  # Replace with your actual sequence
  
  # Calculate the number of mismatches based on the sequence length
  sequence_length = nchar(aso_sequence)
  mismatches_allowed = 2  # int(sequence_length * 0.1)  # 20% of the sequence length
  
  mismatch_conditions <- c("2/")
  
  urls_df <- generate_urls(aso_sequence, mismatch_conditions)
  
  # Initialize an empty data frame that will be filled with summary data
  summary_df <- data.frame()
  
  # Initialize an empty data frame that will be filled with gggenome data
  all_df <- data.frame(
    line = character(),
    equal_signs = integer(),
    mismatches = integer(),
    deletions = integer(), 
    insertions = integer(),
    stringsAsFactors = FALSE
  )
  
  # Initialize an empty vector that will be filled with excel sheetnames
  sheetnames <- c()
  
  # Loop over each URL and retrieve the data from it
  for (i in 1:nrow(urls_df)) {
    url <- urls_df$url[i]
    condition <- urls_df$condition[i]
    mismatch_condition <- urls_df$mismatch_condition[i]
    
    print(paste("Fetching data from:", url))
    response <- GET(url)
    
    if (status_code(response) == 200) {
      text_data <- content(response, "text", encoding = "UTF-8")
      filtered_data_df <- filter_data(text_data, mismatches_allowed)
      all_df <- bind_rows(all_df, filtered_data_df)
    } else {
      print("Failed to fetch data.")
      next
    }
    
    # Extract all unique protein names into a vector
    protein_hits <- str_match(filtered_data_df$line, "NM_\\d+\\.\\d+\\|([^;]+)")[,2]
    protein_hits <- unique(na.omit(protein_hits))
    
    # Loop over all unique proteins and fill the summary_df data frame with the corresponding values
    for (protein_hit in protein_hits) {
      expression_data <- fetch_protein_expression(protein_hit)
      if (!is.null(expression_data) && nrow(expression_data) > 0) {

        protein_row <- filtered_data_df %>% filter(str_detect(line, protein_hit))
        max_equal <- max(protein_row$equal_signs, na.rm = TRUE)
        min_mismatch <- min(protein_row$mismatches, na.rm = TRUE)

        source_sheet <- paste0(str_replace_all(condition, "\\W+", "_"), "_", str_replace_all(mismatch_condition, "\\W+", "_"), "_mismatch")
        
        summary_row <- data.frame(
          'Protein Hit' = protein_hit,
          'Source Sheets' = source_sheet,
          'Equal Signs' = max_equal,
          'Total Mismatches' = min_mismatch,
          stringsAsFactors = FALSE
        )
        
        summary_row <- bind_cols(summary_row, expression_data[1, 2:ncol(expression_data)])
        summary_df <- bind_rows(summary_df, summary_row)
      }
      else {
        print(paste("No expression data found for protein", protein_hit))
      }
    }
    if (nrow(all_df) == 0) {
      all_df <- data.frame(
        line = NULL,
        equal_signs = NULL,
        mismatches = NULL,
        deletions = NULL, 
        insertions = NULL,
        stringsAsFactors = FALSE
      )
    }
    
    # Generate sheetnames
    if (str_detect(str_to_lower(condition), "prespliced")) {
      cond_short <- "_RefSeqCurated_p_d3g2202/"
    } 
    else if (str_detect(str_to_lower(condition), "spliced")) {
      cond_short <- "_RefSeqCurated_s_d3g2202/"
    } 
    else {
      cond_short <- ".p13_d3g2202/"
    }

    clean_cond_short <- str_replace_all(cond_short, "\\W+", "_")
    clean_mismatch <- str_replace_all(mismatch_condition, "\\W+", "_")

    sheetnames <- c(sheetnames, paste0(clean_cond_short, "_", clean_mismatch, "_mm"))
    
  }
  
  # Divide the data frame by condition
  spliced_df <- all_df[str_detect(str_to_lower(all_df$line), "\\|spliced"), ]
  prespliced_df <- all_df[str_detect(str_to_lower(all_df$line), "\\|prespliced"), ]
  other_df <- setdiff(all_df, rbind(spliced_df, prespliced_df))
  
  # Define new column names for the summary data frame
  new_colnames <- c("Protein Hit", "Source Sheets", "Equal Signs", "Total Mismatches",
                    "Gene description", "Brain RNA - amygdala [nTPM]",
                    "Brain RNA - basal ganglia [nTPM]", "Brain RNA - cerebellum [nTPM]",
                    "Brain RNA - cerebral cortex [nTPM]", "Brain RNA - choroid plexus [nTPM]",
                    "Brain RNA - hippocampal formation [nTPM]", "Brain RNA - hypothalamus [nTPM]",
                    "Brain RNA - medulla oblongata [nTPM]", "Brain RNA - midbrain [nTPM]",
                    "Brain RNA - pons [nTPM]", "Brain RNA - spinal cord [nTPM]",
                    "Brain RNA - thalamus [nTPM]")

  colnames(summary_df) <- new_colnames
  
  # Write each filtered data frame to separate sheets in the Excel output file:
  write.xlsx(prespliced_df, file = output_file_name, sheetName=sheetnames[str_detect(sheetnames, "_p_")], row.names = FALSE)
  write.xlsx(spliced_df, file = output_file_name, sheetName=sheetnames[str_detect(sheetnames, "_s_")], append=TRUE, row.names = FALSE)
  write.xlsx(other_df, file = output_file_name, sheetName=sheetnames[!str_detect(sheetnames, "_(p|s)_")], append=TRUE, row.names = FALSE)
  write.xlsx(summary_df, file = output_file_name, sheetName="Summary", append=TRUE, row.names = FALSE)
  
  return(list(urls = urls_df, summary = summary_df, df = all_df, spliced = spliced_df, prespliced = prespliced_df, other = other_df))
}

# ---------- Run ----------
dataframes <- main()
urls <- dataframes$urls
summary_df <- dataframes$summary
df <- dataframes$df
spliced <- dataframes$spliced
prespliced <- dataframes$prespliced
other <- dataframes$other
