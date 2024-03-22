setwd("/mnt/data/Jeremy/txdb/")
path = "/mnt/data/Jeremy/txdb/"

################################################Get human transcript annotations
library(GenomicFeatures)
library(AnnotationDbi)

# obtain dataset as a db file from ensembl
#txdb_hsa = makeTxDbFromBiomart(
  #biomart="ENSEMBL_MART_ENSEMBL",
 # dataset="hsapiens_gene_ensembl",
  #taxonomyId=9606)

# Save Db file (txdb_hsa) locally to path
#saveDb(txdb_hsa,paste0(path, "txdb_hsa_biomart.db"))
