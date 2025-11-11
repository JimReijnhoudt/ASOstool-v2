library(tidyverse)
library(httr)

# Generate URLs for the calculated number of mismatches and all condition
generate_urls <- function(sequence, mismatch_conditions, strands) {
  base_url <- "https://gggenome.dbcls.jp/hg38"
  conditions <- c("_RefSeqCurated_prespliced_d3g2202/",
                  "_RefSeqCurated_spliced_d3g2202/", ".p13_d3g2202/")
  
  mismatch_col <- character()
  strand_col <- character()
  condition_col <- character()
  url_col <- character()
  
  for (st in strands) {
    for (mc in mismatch_conditions) {
      for (cond in conditions) {
        url <- paste0(base_url, cond, mc, st, sequence, ".txt")
        
        mismatch_col <- c(mismatch_col, mc)
        strand_col <- c(strand_col, st)
        condition_col <- c(condition_col, cond)
        url_col <- c(url_col, url)
      }
    }
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


all <- function(sequence) {
  # Calculate the number of mismatches based on the sequence length
  sequence_length = nchar(sequence)
  mismatches_allowed = 2  # int(sequence_length * 0.1)  # 20% of the sequence length
  
  mismatch_conditions <- c("2/")
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
      text_data <- content(response, "text", encoding = "UTF-8")
      filtered_data_df <- filter_data(text_data, mismatches_allowed)
      all_df <- bind_rows(all_df, filtered_data_df)
    } else {
      print("Failed to fetch data.")
      next
    }
  }
  return(all_df)
}

sequences <- readLines("sequenties.txt")

run_all_on_sequence <- function(seq) {
  df <- all(seq)
  tibble(
    sequentie = seq,
    hits = nrow(df)
  )
}

results <- sequences %>%
  lapply(run_all_on_sequence) %>%
  bind_rows()