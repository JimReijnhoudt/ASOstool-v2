library(tidyverse)
library(httr)
# library(readr)
# library(stringr)
# library("xlsx")

# Generate URLs for the calculated number of mismatches and all condition
generate_urls <- function(sequence, mismatch_conditions, strands) {
  base_url <- "https://gggenome.dbcls.jp/hg38"
  conditions <- c("_RefSeqCurated_prespliced_d3g2202/",
                  "_RefSeqCurated_spliced_d3g2202/", ".p13_d3g2202/")
  
  mismatch_col <- character()
  strand_col <- character()
  condition_col <- character()
  url_col <- character()
  
  for (cond in conditions) {
    url <- paste0(base_url, cond, mismatch_conditions, strands, sequence, ".txt")
    
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
  start <- c('NM_', 'chr', 'NR_')
  lines <- str_split(data, "\n")[[1]]
  
  line_col <- character()
  equal_col <- integer()
  mismatch_col <- integer()
  ins_col <- integer()
  del_col <- integer()
  sub_col <- character()
  query_col <- character()
  
  start_target_col <- integer()
  end_target_col <- integer()
  snippet_col <- character()
  snippet_start_col <- integer()
  snippet_end_col <- integer()
  
  
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
        
        subs <- as.character(parts[9])
        quers <- as.character(parts[8])
        
        sub_col <- c(sub_col, subs)
        query_col <- c(query_col, quers)
        
        start_target_col <- c(start_target_col, as.integer(parts[3]))
        end_target_col <- c(end_target_col, as.integer(parts[4]))
        snippet_col <- c(snippet_col, as.character(parts[5]))
        snippet_start_col <- c(snippet_start_col, as.integer(parts[6]))
        snippet_end_col <- c(snippet_end_col, as.integer(parts[7]))
      }
    }
  }
  
  df <- data.frame(
    line = line_col,
    equal_signs = equal_col,
    mismatches = mismatch_col,
    deletions = del_col, 
    insertions = ins_col,
    subject_seq = sub_col,
    query_seq = query_col,
    start_target = start_target_col,
    end_target = end_target_col,
    snippet = snippet_col,
    snippet_start = snippet_start_col,
    snippet_end = snippet_end_col,
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
        
        min_deletions <- min(protein_row$deletions, na.rm = TRUE)
        min_insertions <- min(protein_row$insertions, na.rm = TRUE)
        
        chosen_row <- protein_row %>%
          filter(
            equal_signs == max_equal,
            mismatches == min_mismatch,
            deletions == min_deletions,
            insertions == min_insertions
          ) %>%
          slice(1)
        
        subject_seq <- chosen_row$subject_seq
        query_seq   <- chosen_row$query_seq
        
        start_target <- chosen_row$start_target
        end_target <- chosen_row$end_target
        snippet <- chosen_row$snippet
        snippet_start <- chosen_row$snippet_start
        snippet_end <- chosen_row$snippet_end
        
        source_sheet <- paste0(str_replace_all(condition, "\\W+", "_"), "_", str_replace_all(mismatch_condition, "\\W+", "_"), "_mismatch")
        
        summary_row <- data.frame(
          'Protein Hit' = protein_hit,
          'Source Sheets' = source_sheet,
          'Equal Signs' = max_equal,
          'Total Mismatches' = min_mismatch,
          'Total Deletions' = min_deletions, 
          'Total Insertions' = min_insertions,
          'Subject sequence' = subject_seq,
          'Query sequence' = query_seq,
          'start_target' = start_target,
          'end_target' = end_target,
          'snippet' = snippet,
          'snippet_start' = snippet_start,
          'snippet_end' = snippet_end,
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
        subject_seq = NULL,
        query_seq = NULL,
        start_target = NULL,
        end_target = NULL,
        snippet = NULL,
        snippet_start = NULL,
        snippet_end = NULL,
        stringsAsFactors = FALSE
      )
    }
  }
  
  # Divide the data frame by condition
  spliced_df <- all_df[str_detect(str_to_lower(all_df$line), "\\|spliced"), ]
  prespliced_df <- all_df[str_detect(str_to_lower(all_df$line), "\\|prespliced"), ]
  other_df <- setdiff(all_df, rbind(spliced_df, prespliced_df))
  
  # Define new column names for the summary data frame
  new_colnames <- c("Protein Hit", "Source Sheets", "Equal Signs", "Total Mismatches", "Total Deletions", "Total Insertions", "Subject sequence", "Query sequence",
                    "start_target", "end_target", "snippet", "snippet_start", "snippet_end",
                    "Gene description", "Brain RNA - amygdala [nTPM]",
                    "Brain RNA - basal ganglia [nTPM]", "Brain RNA - cerebellum [nTPM]",
                    "Brain RNA - cerebral cortex [nTPM]", "Brain RNA - choroid plexus [nTPM]",
                    "Brain RNA - hippocampal formation [nTPM]", "Brain RNA - hypothalamus [nTPM]",
                    "Brain RNA - medulla oblongata [nTPM]", "Brain RNA - midbrain [nTPM]",
                    "Brain RNA - pons [nTPM]", "Brain RNA - spinal cord [nTPM]",
                    "Brain RNA - thalamus [nTPM]")
  
  colnames(summary_df) <- new_colnames
  
  return(summary_df)
}

acc_snippet <- function(
    begin_target,
    end_target,
    begin_snippet,
    end_snippet,
    full_snippet,
    total_length = 80
) {
  # 1. Targetlengte
  target_length <- end_target - begin_target + 1
  if (target_length > total_length) {
    stop("Target is langer dan snippet.")
  }
  
  # 2. Flanks
  flank_total <- total_length - target_length
  flank_left  <- floor(flank_total / 2)
  flank_right <- ceiling(flank_total / 2)
  
  # 3. Eerste snippet (extern)
  snippet_start <- begin_target - flank_left
  snippet_end   <- end_target   + flank_right
  
  # 4. Linkergrens
  if (snippet_start < begin_snippet) {
    shift <- begin_snippet - snippet_start
    snippet_start <- begin_snippet
    snippet_end   <- snippet_end + shift
  }
  
  # 5. Rechtergrens
  if (snippet_end > end_snippet) {
    shift <- snippet_end - end_snippet
    snippet_end   <- end_snippet
    snippet_start <- snippet_start - shift
  }
  
  # 6. Sanity check lengte
  if ((snippet_end - snippet_start + 1) != total_length) {
    stop("Kon geen snippet van 80 nt maken binnen grenzen.")
  }
  
  # 7. Targetposities intern (1–80)
  target_start_internal <- begin_target - snippet_start + 1
  target_end_internal   <- end_target   - snippet_start + 1
  
  # 8. Substring indices in full_snippet
  start_i <- snippet_start - begin_snippet + 1
  end_i   <- snippet_end   - begin_snippet + 1
  
  # 9. Extract sequentie
  snippet_seq <- substr(full_snippet, start_i, end_i)
  
  # 10. Return alles
  return(
    list(
      snippet_seq             = snippet_seq,
      snippet_start_external  = snippet_start,
      snippet_end_external    = snippet_end,
      snippet_length          = total_length,
      target_start_internal   = target_start_internal,
      target_end_internal     = target_end_internal
    )
  )
}

# # ---------- Run ----------
# dataframes <- all_offt("TTTTTGCCATCCTGGGCGCT", 2)
# urls <- dataframes$urls
# summary_df <- dataframes$summary
# df <- dataframes$df
# spliced <- dataframes$spliced
# prespliced <- dataframes$prespliced
# other <- dataframes$other