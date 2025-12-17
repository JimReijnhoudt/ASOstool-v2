library(httr)
library(jsonlite)
library(dplyr)

omim_key <- readLines("../OMIM_key")

OMIM_search <- function(gene) {
  base_url <- "https://api.omim.org/api/entry/search?search="

  params <- paste0("approved_gene_symbol:", gene, "&start=0&limit=1&include=text&include=geneMap&format=json&apiKey=", omim_key)
  full_url <- paste0(base_url, params)
  
  dat <- fromJSON(full_url)
  
  if (length(dat$omim$searchResponse$entryList) == 0) {
    return(tibble(
      gene = gene,
      mim_number = NA,
      n_diseases = NA
    ))
  }
  
  entry <- dat$omim$searchResponse$entryList[1]$entry
  # OMIM entry nummer (gen/ziekte)
  mim_number <- entry$mimNumber
  
  # aantal ziekten (phenotypes) die aan dit gen gekoppeld zijn
  n_diseases <- length(entry$geneMap$phenotypeMapList)
  
  res <- tibble(
    gene = gene,
    mim_number = mim_number,
    n_diseases = n_diseases
  )
  return(res)
}
