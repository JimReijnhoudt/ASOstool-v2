library(readxl)
library(tidyr)
library(glue)
library(dplyr)

rnaseh_results <- function(selected_row_name, mod_5prime, mod_3prime){
  # The standard score matrix for the nucleotide positions for cleavage. 
  nucleotide_model <- read_excel("../Files/nucleotide_model.xlsx")
  
  # Target sequence. 
  aso_seq = selected_row_name
  window_size = 9
  aso_len = nchar(aso_seq)
  
  # Check if windows are possible by size.
  start_min  <- mod_5prime + 1
  start_max <- aso_len - mod_3prime - window_size + 1
  
  if (start_min > start_max) {
    showNotification(
      "Modifications to 5' and 3' ends overlap to much for sequence length.",
      type = "error",
      closeButton = TRUE
    )
    return(NULL)
  }
  
  # All windows of ASO sequence.
  starts_pos <- 1:(aso_len - window_size + 1)
  windows <- substring(aso_seq, starts_pos, starts_pos + (window_size - 1))
  
  # Remove all overlapping windows.
  valid_windows <- which(starts_pos >= start_min & starts_pos <= start_max)
  starts_pos <- starts_pos[valid_windows]
  windows <- windows[valid_windows]
  
  # To get the dinucleotides from each of the windows.  
  get_pairs <- function(window) {
    sapply(1:(nchar(window) - 1), function(i) substring(window, i, i + 1))
  }
  
  # Builds the ASO-sequence data-frame. This includes all the windows and dinucleotide pairs. 
  seq_df <- t(sapply(windows, function(w) c(window = w, get_pairs(w))))
  seq_df <- as.data.frame(seq_df, stringsAsFactors = FALSE)
  
  # Add pairs to dataframe.
  seq_df <- seq_df %>% mutate(pairs = glue("{starts_pos} - {starts_pos + (window_size - 1)}"))
  seq_df <- seq_df %>% relocate(pairs, .after = window)
  
  # To rename the columns for ease of use.
  colnames(seq_df) <- c("window", "pairs", paste0("p", 1:(window_size - 1), 2:window_size))
  
  # Mapping the values from the seq_df with nucleotide_model, first make a long dataframe of both. 
  nucleotide_long <- nucleotide_model %>% 
    pivot_longer(-dinuc_pairs, names_to = "position", values_to = "value")
  
  seq_long <- seq_df %>% 
    pivot_longer(-(window:pairs), names_to = "position", values_to = "dinuc_pairs")
  
  # Joined in long format.
  joined_long <- seq_long %>% left_join(nucleotide_long, by=c("dinuc_pairs", "position"))
  
  # Joined in short format. This is used for the final evaluation. 
  joined_short <- joined_long %>% 
    dplyr::select(window, pairs, position, value) %>% 
    pivot_wider(names_from = position, values_from = value) %>%
    arrange(window)
  
  # Final calculation for each row of the average and the multiplied score. 
  joined_short$average = rowMeans(joined_short[, !(names(joined_short) %in% c("window", "pairs"))])
  joined_short$multiply <- apply (joined_short[, !(names(joined_short) %in% c("window", "pairs", "average"))], 1, prod)
  
  # Adding the final two colomns from the joined_short to seq_df for a extra viewing option. 
  seq_scored <- seq_df %>% 
    left_join(joined_short %>% dplyr::select (window, pairs, average, multiply), by = c("window", "pairs"))
  
  # Order the dataframe by average.
  joined_short <- joined_short[order(joined_short$average, decreasing = TRUE),]
  seq_scored <- seq_scored[order(seq_scored$average, decreasing = TRUE),]
  
  # Top 10 results. 
  top_results <- head(joined_short, 10)
  
  # Simplified results. 
  simp_results <- joined_short[, c("window", "pairs", "average")]
  colnames(simp_results)[colnames(simp_results) == "pairs"] <- "position"
  
  return(simp_results)
}