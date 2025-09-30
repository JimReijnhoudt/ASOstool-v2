path <- getwd()

# Build a subfolder path dynamically.
location <- file.path(path, "Data")
db_file <- file.path(location, "txdb_hsa.db")

# Create the folder if it doesn't exist.
if (!dir.exists(location)) {
  dir.create(location, recursive = TRUE)
}

# Build the full destination path for the file.
destination <- file.path(location, "Homo_sapiens.GRCh38.115.gtf.gz")

url <- "https://ftp.ensembl.org/pub/current_gtf/homo_sapiens/Homo_sapiens.GRCh38.115.gtf.gz"

# Download the file to the specified location.
download.file(url, destfile = destination, mode = "wb")

# Check if file is downloaded and if it is downloaded make a db file with it.
if (file.exists(destination)) {
  message("Download successful: ", destination)
  
  txdb_hsa <- makeTxDbFromGFF(destination, format="gtf")
  saveDb(txdb_hsa, db_file)
  message("Done, saved to: ", db_file)
  
} else {
  warning("Download failed.")
}