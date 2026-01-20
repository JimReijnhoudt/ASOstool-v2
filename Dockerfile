FROM rocker/tidyverse:4.5.1

LABEL maintainer="you@example.com"

RUN /rocker_scripts/install_shiny_server.sh 1.5.23.1030


# 1. System dependencies
##############################################
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    wget \
    libxml2-dev \
    libssl-dev \
    libcurl4-openssl-dev \
    libgsl-dev \
    python3 \
    zlib1g-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install GIT
RUN apt-get -y install git

# Install SSH client tools for GitHub

RUN apt-get update && apt-get install -y \
    openssh-client \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*


# Configure git for RStudio user
##############################################
RUN git config --system user.name "rstudio" && \
    git config --system user.email "rstudio@example.com"

# ViennaRNA 2.7.0 installation
##############################################
RUN wget https://www.tbi.univie.ac.at/RNA/download/sourcecode/2_7_x/ViennaRNA-2.7.0.tar.gz && \
    tar xzf ViennaRNA-2.7.0.tar.gz && \
    cd ViennaRNA-2.7.0 && \
    ./configure && \
    make -j$(nproc) && \
    make install && \
    ldconfig && \
    cd .. && \
    rm -rf ViennaRNA-2.7.0 ViennaRNA-2.7.0.tar.gz


# Install R packages (CRAN + Bioconductor)
##############################################
RUN R -e "install.packages('BiocManager', repos='https://cloud.r-project.org')" \
 && R -e "BiocManager::install(c( \
        'GenomicFeatures', \
        'AnnotationDbi', \
        'BSgenome.Hsapiens.NCBI.GRCh38', \
        'biomaRt', \
        'Biostrings', \
        'txdbmaker' \
     ), ask = FALSE)" \
 && R -e "install.packages(c( \
        'shinythemes', 'shiny', 'tidyverse', 'cluster', \
        'rlang', 'dplyr', 'DT', 'shinyBS', 'shinydashboard', 'openxlsx', \
        'future', 'future.apply' \
     ), repos='https://cloud.r-project.org', dependencies=TRUE)"


# Make necesary directories
##############################################
RUN mkdir -p /home/rstudio/ASOstool-v2 && \
    chown -R rstudio:rstudio /home/rstudio

RUN mkdir -p /srv/shiny-server/ASOstool-v2 && \
    chown -R shiny:shiny /srv/shiny-server/ && \
    chmod -R 775 /srv/shiny-server/ASOstool-v2

RUN mkdir -p /home/rstudio/.ssh \
 && chown -R rstudio:rstudio /home/rstudio/.ssh \
 && chmod 700 /home/rstudio/.ssh


# Add run as rstudio for shiny-server
##############################################
RUN sed -i 's/^run_as.*/run_as rstudio;/' /etc/shiny-server/shiny-server.conf || \
    sed -i '1irun_as rstudio;' /etc/shiny-server/shiny-server.conf
    

# Download and create txdb database
##############################################
RUN mkdir -p /opt/ASOstool-v2


RUN R -e "\
  library(txdbmaker); \
  gtf <- 'Homo_sapiens.GRCh38.115.gtf.gz'; \
  gtf_url <- paste0('https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/', gtf); \
  download.file(gtf_url, gtf, mode = 'wb'); \
  txdb <- makeTxDbFromGFF(gtf); \
  saveDb(txdb, '/opt/ASOstool-v2/txdb_hsa_biomart.db'); \
  unlink(gtf); \
"


# 5. Expose ports
##############################################
EXPOSE 3838 8787

# 6. Default user credentials
##############################################
ENV USER=rstudio
ENV PASSWORD=rstudio


# 7. Start RStudio Server
##############################################
CMD ["/init"]
