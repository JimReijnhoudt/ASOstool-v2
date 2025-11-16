FROM rocker/tidyverse:4.5.1

LABEL maintainer="you@example.com"

RUN /rocker_scripts/install_shiny_server.sh 1.5.23.1030


# 1. System dependencies

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
RUN git config --system user.name "rstudio" && \
    git config --system user.email "rstudio@example.com"

##############################################
# 2. ViennaRNA 2.7.0 installation
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

##############################################
# 3. Install R packages (CRAN + Bioconductor)
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
        'rlang', 'dplyr', 'DT', 'shinyBS', 'shinydashboard' \
     ), repos='https://cloud.r-project.org', dependencies=TRUE)"

##############################################
# 4. Make necesary directories
##############################################
RUN mkdir -p /home/rstudio && \
    chown -R rstudio:rstudio /home/rstudio

# Create shiny-server directory structure
RUN mkdir -p /srv/shiny-server/app && \
    chown -R shiny:shiny /srv/shiny-server/app

RUN mkdir -p /home/rstudio/.ssh \
 && chown -R rstudio:rstudio /home/rstudio/.ssh \
 && chmod 700 /home/rstudio/.ssh

##############################################
# 5. Expose ports
##############################################
EXPOSE 3838 8787

##############################################
# 6. Default user credentials
##############################################
ENV USER=rstudio
ENV PASSWORD=rstudio

##############################################
# 7. Start RStudio Server
##############################################
CMD ["/init"]
