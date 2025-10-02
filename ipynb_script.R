getwd()
setwd("C:/Users/fayef/Documents/BI/BI3/periode_1/XEXT/ASOstool-v2")

install.packages("tidyverse")
library(tidyverse)

# Define the file name at the top
output_file_name <- "20250801_EH_UNC13A_SNP151&41_ASO146.xlsx"

# Define the ASO sequence
aso_sequence <-  "CACACCATGCACATTCAA"  # Replace with your actual sequence

# Calculate the number of mismatches based on the sequence length
sequence_length = nchar(aso_sequence)
mismatches_allowed = 2  # int(sequence_length * 0.1)  # 20% of the sequence length

mismatch_conditions <- list("2/")

generate_urls <- function(aso_sequence, mismatch_conditions) {
  base_url <- "https://gggenome.dbcls.jp/hg38"
  conditions <- c("_RefSeqCurated_prespliced_d3g2202/",
                "_RefSeqCurated_spliced_d3g2202/", ".p13_d3g2202/")
  urls = c()
  for (mc in mismatch_conditions) {
    for (cond in conditions) {
      url <- paste0(base_url, cond, mc, aso_sequence, ".txt")
      urls <- c(urls, url)
    }
  }
  return(urls)
}

filter_data <- function(data, criterion) {
  filtered_data <- c()
  for (line in str_split(data, "\n")[[1]]) {
    if (any(str_starts(line, c("NM_", "chr", "NR_")))) {
      equal_signs <- str_count(line, "=")
      if (equal_signs >= criterion) {
        filtered_data <- c(filtered_data, line)
      }
    }
  }
  return(filtered_data)
}