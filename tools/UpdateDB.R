# setwd("/mnt/data/Jeremy/txdb/")
# path = "/mnt/data/Jeremy/txdb/"

setwd("C:/Users/fayef/Documents/BI/BI3/periode_1/XEXT/ASOstool-v2")
path <- "C:/Users/fayef/Documents/BI/BI3/periode_1/XEXT/ASOstool-v2"

################################################Get human transcript annotations
library(GenomicFeatures)
library(AnnotationDbi)

# obtain dataset as a db file from ensembl
# txdb_hsa = makeTxDbFromBiomart(
#   biomart="ENSEMBL_MART_ENSEMBL",
#   dataset="hsapiens_gene_ensembl",
#   taxonomyId=9606)

txdb <- txdbmaker::makeTxDbFromGFF("Homo_sapiens.GRCh38.112.gtf.gz") 

# Save Db file (txdb_hsa) locally to path
# saveDb(txdb_hsa,paste0(path, "txdb_hsa_biomart.db"))

saveDb(txdb, "txdb_hsa_biomart.db") 
print("done")
