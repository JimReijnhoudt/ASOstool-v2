library(tidyverse)
library(httr)
library(future)
library(future.apply)
# library(readr)
# library(stringr)
# library("xlsx")

# Generate URLs for the calculated number of mismatches and all condition
generate_urls <- function(sequence, mismatch_conditions, strands) {
  base_url <- "https://gggenome.dbcls.jp/hg38"
  conditions <- c("_RefSeqCurated_prespliced_d3g2202/",
                  "_RefSeqCurated_spliced_d3g2202/")

  paste0(base_url, conditions, mismatch_conditions, strands, sequence, ".json")
}

# Function to filter data based on the number of equal signs
filter_data <- function(data, criterion) {
  df_filtered <- data$results
  
  if (is.null(df_filtered) || !is.data.frame(df_filtered) || nrow(df_filtered) == 0) {
    return(data.frame())
  }
  
  df <- data.frame(
    line = df_filtered$name,
    match_string = df_filtered$edit,
    matches = str_count(df_filtered$edit, "="),
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
  df <- df %>%
    mutate(gene_name = str_extract(line, "(?<=\\|)[^;]+")) %>%
    distinct(gene_name, match_string, query_seq, .keep_all = TRUE)
  return(df)
}


all_offt <- function(sequence, mismatches_allowed) {
  # Calculate the number of mismatches based on the sequence length
  sequence_length = nchar(sequence)
  mismatch_conditions <- paste0(mismatches_allowed, "/")
  
  strands <- "+/"
  urls_vec <- generate_urls(sequence, mismatch_conditions, strands)
  
  # Initialize an empty data frame that will be filled with summary data
  out_list <- vector("list", length(urls_vec))
  
  # Loop over each URL and retrieve the data from it
  for (i in seq_along(urls_vec)) {
    url <- urls_vec[[i]]
    # condition <- urls_df$condition[i]
    # mismatch_condition <- urls_df$mismatch_condition[i]

    response <- tryCatch(
      RETRY(
        "GET",
        url,
        times = 3,
        pause_base = 2,
        pause_cap = 10,
        terminate_on = c(400, 404),
        timeout(120)
     ),
     error = function(e) NULL
    )
      
    if (is.null(response) || status_code(response) != 200) {
      out_list[[i]] <- NULL
      next
    }
    

    df_json <- fromJSON(
      content(response, as = "text", encoding = "UTF-8"),
      flatten = TRUE)
    
    
    filtered_data_df <- filter_data(df_json, mismatches_allowed)
    
    
    if (!is.data.frame(filtered_data_df) || nrow(filtered_data_df) == 0) {
      out_list[[i]] <- NULL
    } else {
      out_list[[i]] <- filtered_data_df
    }
  }
  
  bind_rows(Filter(Negate(is.null), out_list))
}
