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
print("HS: ")
print(HS)

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

# uni_tar = dplyr::select(target_annotation, name, length)%>%
#   unique() %>%
#   split(.,.$length)
# 
# offtarget_results_list <- list()
# 
# uni_tar <- lapply(seq_along(uni_tar), function(k) {
#   
#   X <- uni_tar[[k]]
#   
#   fwd <- DNAStringSet(X$name)
#   rev <- reverseComplement(fwd)
#   
#   dict0_fwd <- PDict(fwd, max.mismatch = 0)
#   dict0_rev <- PDict(rev, max.mismatch = 0)
#   
#   dict1_fwd <- PDict(fwd, max.mismatch = 1)
#   dict1_rev <- PDict(rev, max.mismatch = 1)
#   
#   pm_fwd  <- vwhichPDict(dict0_fwd, HS, max.mismatch = 0, min.mismatch = 0)
#   mm1_fwd <- vwhichPDict(dict1_fwd, HS, max.mismatch = 1, min.mismatch = 1)
#   
#   pm_rev  <- vwhichPDict(dict0_rev, HS, max.mismatch = 0, min.mismatch = 0)
#   mm1_rev <- vwhichPDict(dict1_rev, HS, max.mismatch = 1, min.mismatch = 1)
# 
#   X$gene_hits_pm  <- tabulate(unlist(pm_fwd),  nbins=nrow(X)) +
#     tabulate(unlist(pm_rev),  nbins=nrow(X))
#   
#   X$gene_hits_1mm <- tabulate(unlist(mm1_fwd), nbins=nrow(X)) +
#     tabulate(unlist(mm1_rev), nbins=nrow(X))
#   
# 
#   make_df <- function(hitlist, strand, mismatch) {
#     lapply(seq_along(hitlist), function(i) {
#       if (length(hitlist[[i]]) > 0) {
#         tibble(
#           query_seq = X$name[i],
#           ensembl_id = names(HS)[hitlist[[i]]],
#           strand = strand,
#           mismatch = mismatch
#         )
#       } else NULL
#     }) %>% bind_rows()
#   }
#   
# 
#   df_pm_fwd  <- make_df(pm_fwd,  "+", 0)
#   df_mm1_fwd <- make_df(mm1_fwd, "+", 1)
#   
# 
#   df_pm_rev  <- make_df(pm_rev,  "-", 0)
#   df_mm1_rev <- make_df(mm1_rev, "-", 1)
#   
#   offtarget_results_list[[k]] <- bind_rows(
#     df_pm_fwd, df_mm1_fwd,
#     df_pm_rev, df_mm1_rev
#   )
#   
#   return(X)
# })
# 
# # Combineer alles weer zoals origineel
# uni_tar <- bind_rows(uni_tar)
# 
# # Combineer alle off-target hits in één dataframe
# offtarget_results <- bind_rows(offtarget_results_list) %>%
#   distinct()
# 
# print("offtarget_results: ")
# print(offtarget_results)
# print("uni_tar: ")
# print(uni_tar)

# -------------------------------------------------------
uni_tar = dplyr::select(target_annotation, name, length) %>%
  unique() %>%
  split(., .$length)

print("uni_tar: ")
print(uni_tar)

offtarget_results_list <- list()

uni_tar <- lapply(seq_along(uni_tar), function(k) {
  print("k")
  print(k)

  X <- uni_tar[[k]]
  print("X:")
  print(head(X))

  X$name <- as.character(X$name)
  attributes(X$name) <- NULL

  # Reverse complement names
  X$name_rev <- as.character(
    reverseComplement(DNAStringSet(X$name))
  )

  # -----------------------------
  # Build dictionaries
  # -----------------------------
  dict_fwd0 <- PDict(X$name, max.mismatch = 0)
  dict_rev0 <- PDict(X$name_rev, max.mismatch = 0)

  dict_fwd1 <- PDict(X$name, max.mismatch = 1)
  dict_rev1 <- PDict(X$name_rev, max.mismatch = 1)

  # -----------------------------
  # Search HS (forward + reverse)
  # -----------------------------
  pm_fwd  <- vwhichPDict(dict_fwd0, HS)
  pm_rev  <- vwhichPDict(dict_rev0, HS)
  mm1_fwd <- vwhichPDict(dict_fwd1, HS)
  mm1_rev <- vwhichPDict(dict_rev1, HS)

  # Combine forward+reverse hits per query
  pm_all  <- Map(c, pm_fwd, pm_rev)
  mm1_all <- Map(c, mm1_fwd, mm1_rev)

  # -----------------------------
  # Store hit counts
  # -----------------------------
  X$gene_hits_pm  <- tabulate(unlist(pm_all),  nbins=nrow(X))
  X$gene_hits_1mm <- tabulate(unlist(mm1_all), nbins=nrow(X))

  # -----------------------------
  # Convert hits into a dataframe
  # -----------------------------
  make_df <- function(hitlist, mm_value) {
    bind_rows(lapply(seq_along(hitlist), function(i) {
      if (length(hitlist[[i]]) > 0) {
        tibble(
          query_seq = X$name[i],
          ensembl_id = names(HS)[hitlist[[i]]],
          mismatch = mm_value
        )
      } else NULL
    }))
  }

  pm_df  <- make_df(pm_all, 0)
  mm1_df <- make_df(mm1_all, 1)

  offtarget_results_list[[k]] <- bind_rows(pm_df, mm1_df)

  return(X)
})

uni_tar <- bind_rows(uni_tar)

offtarget_results <- bind_rows(offtarget_results_list) %>%
  distinct()

print("offtarget_results:")
print(offtarget_results)
# -------------------------------------------------------

# dict <- PDict("GCTTTTTGCCATCCT", max.mismatch = 0)
# pm <- vwhichPDict(dict, HS)
# names(HS)[unlist(pm)]
# 
# dict0 <- PDict(target_annotation$name, max.mismatch = 0)
# pm <- vwhichPDict(dict0, HS)
# lengths(pm)
# 
# grep("[^ACGT]", target_annotation$name, value = TRUE)
# str(target_annotation$name)
# 
# dict_fwd_ <- PDict(target_annotation$name, max.mismatch = 0)
# dict_rev_ <- PDict(as.character(reverseComplement(DNAStringSet(target_annotation$name))), max.mismatch = 0)
# 
# pm_fwd_ <- vwhichPDict(dict_fwd_, HS)
# pm_rev_ <- vwhichPDict(dict_rev_, HS)
# 
# # Combineer beide richtingen
# pm_all_ <- Map(c, pm_fwd_, pm_rev_)
# names(HS)[unique(unlist(pm_all_))]