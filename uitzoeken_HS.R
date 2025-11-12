setwd("C:/Users/fayef/Documents/BI/BI3/periode_1/XEXT/ASOstool-v2/")

library(GenomicFeatures)

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

HS
