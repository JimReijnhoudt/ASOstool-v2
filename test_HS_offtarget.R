library(shinythemes)
library(shiny)
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
library(txdbmaker)

library(tibble)



txdb_hsa <- loadDb("txdb_hsa_biomart.db")

print("milestone1")

# Extract the genes
gdb_hsa <- genes(txdb_hsa)

# Define the chromosomes to keep
chr_to_keep <- c(as.character(1:22),'X','Y','MT')

# filtert Hsapiens zodat het alleen de aangegeven chromosomen pakt
Hsapiens@user_seqnames <- setNames(chr_to_keep,chr_to_keep)
Hsapiens@seqinfo <- Hsapiens@seqinfo[chr_to_keep]

# Subset genes to keep only those on specified chromosomes
gdb_hsa <- gdb_hsa[seqnames(gdb_hsa) %in% chr_to_keep]

print("milestone2")

# Get the sequences *
HS <- getSeq(Hsapiens, gdb_hsa)

target_annotation <- tibble::tibble(
  start = c(1, 3, 4, 5, 6, 7, 8, 11, 12, 15),
  length = c(15, 15, 15, 15, 15, 15, 15, 15, 15, 15),
  end = c(15, 17, 18, 19, 20, 21, 22, 25, 26, 29),
  name = c(
    "GCTTTTTGCCATCCT",
    "TTTTTGCCATCCTGG",
    "TTTTGCCATCCTGGG",
    "TTTGCCATCCTGGGC",
    "TTGCCATCCTGGGCG",
    "TGCCATCCTGGGCGC",
    "GCCATCCTGGGCGCT",
    "ATCCTGGGCGCTAGA",
    "TCCTGGGCGCTAGAG",
    "TGGGCGCTAGAGGGA"
  ),
  NoRepeats = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
)

uni_tar = dplyr::select(target_annotation, name, length)%>%
  unique() %>%
  split(.,.$length)

offtarget_results_list <- list()

uni_tar <- lapply(seq_along(uni_tar), function(k) {
  X <- uni_tar[[k]]
  
  dict0 <- PDict(X$name, max.mismatch = 0)
  dict1 <- PDict(X$name, max.mismatch = 1)
  
  pm <- vwhichPDict(dict0, HS, max.mismatch = 0, min.mismatch = 0)
  
  
  X$gene_hits_pm <- tabulate(unlist(pm), nbins = nrow(X))
  
  mm1 <- vwhichPDict(dict1, HS, max.mismatch = 1, min.mismatch = 1)
  X$gene_hits_1mm <- tabulate(unlist(mm1), nbins = nrow(X))
  
  # ---- Nieuw stuk: sla per-hit info op ----
  pm_df <- lapply(seq_along(pm), function(i) {
    if (length(pm[[i]]) > 0) {
      tibble(
        query_seq = X$name[i],
        ensembl_id = names(HS)[pm[[i]]],
        mismatch = 0
      )
    } else NULL
  }) %>% bind_rows()
  
  mm1_df <- lapply(seq_along(mm1), function(i) {
    if (length(mm1[[i]]) > 0) {
      tibble(
        query_seq = X$name[i],
        ensembl_id = names(HS)[mm1[[i]]],
        mismatch = 1
      )
    } else NULL
  }) %>% bind_rows()
  
  # Combineer 0- en 1-mismatch hits
  offtarget_results_list[[k]] <- bind_rows(pm_df, mm1_df)
  
  return(X)
})

# Combineer alles weer zoals origineel
uni_tar <- bind_rows(uni_tar)

# Combineer alle off-target hits in één dataframe
offtarget_results <- bind_rows(offtarget_results_list) %>%
  distinct()

print(offtarget_results)
print(uni_tar)

# -------------------------------------------------------

dict <- PDict("GCTTTTTGCCATCCT", max.mismatch = 0)
pm <- vwhichPDict(dict, HS)
names(HS)[unlist(pm)]

dict0 <- PDict(target_annotation$name, max.mismatch = 0)
pm <- vwhichPDict(dict0, HS)
lengths(pm)

grep("[^ACGT]", target_annotation$name, value = TRUE)
str(target_annotation$name)

dict_fwd <- PDict(target_annotation$name, max.mismatch = 0)
dict_rev <- PDict(as.character(reverseComplement(DNAStringSet(target_annotation$name))), max.mismatch = 0)

pm_fwd <- vwhichPDict(dict_fwd, HS)
pm_rev <- vwhichPDict(dict_rev, HS)

# Combineer beide richtingen
pm_all <- Map(c, pm_fwd, pm_rev)
names(HS)[unique(unlist(pm_all))]