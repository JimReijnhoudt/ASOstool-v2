# Function to install packages if not already installed
install_if_missing <- function(package_name) {
  if (!package_name %in% installed.packages()[, "Package"]) {
    message(paste("Installing", package_name, "..."))
    tryCatch({
      BiocManager::install(package_name, ask = FALSE)
    }, error = function(e) {
      # If BiocManager fails (e.g., package not available on Bioconductor),
      # install from CRAN using install.packages
      install.packages(package_name, dependencies = TRUE)
    })
    message(paste(package_name, "installed successfully."))
  } else {
    message(paste(package_name, "is already installed."))
  }
}

# List of packages to check and install
packages <- c("shinythemes", "shiny", "GenomicFeatures", "AnnotationDbi",
              "BSgenome.Hsapiens.NCBI.GRCh38", "biomaRt", "Biostrings",
              "tidyverse", "cluster", "rlang", "dplyr", "purrr", "DT", "shinyBS", "txdbmaker", "openxlsx")

# Install BiocManager if not already installed
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Check and install each package
lapply(packages, install_if_missing)

# Load the packages
library(DT)
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
library(purrr)
library(shinyBS)
library(txdbmaker)
library(openxlsx)

print("Done")
 