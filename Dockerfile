FROM rocker/r-ver:4.4.1

ARG P3MVER='jammy/2024-09-11'
# Add P3MVER as default repository for CRAN in .Rprofile
RUN touch /root/.Rprofile \
    && echo "options(repos=c(CRAN='"https://p3m.dev/cran/__linux__/${P3MVER}"'))" >> /root/.Rprofile
    
ARG BIOCONDUCTORVER='3.19'
# Set BIOCONDUCTORVER as biocondutor version and add P3MVER as default repository for BIOCONDUCTOR in .Rprofile
RUN mkdir -p /tmp/downloaded_packages/ && echo ${BIOCONDUCTORVER} > /tmp/downloaded_packages/BIOCONDUCTORVER.txt \
    && tmp="https://packagemanager.posit.co/bioconductor/"$( echo ${P3MVER} | cut -f2 -d/ ) \
    && echo "options(BioC_mirror='${tmp}')" >> /root/.Rprofile \
    && echo "options(BIOCONDUCTOR_CONFIG_FILE='${tmp}/config.yaml')" >> /root/.Rprofile
# ------------------------------------------------------------
# System libraries
# ------------------------------------------------------------
RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev libssl-dev libxml2-dev \
    libpng-dev libjpeg-dev libtiff-dev \
    libfreetype6-dev libharfbuzz-dev libfribidi-dev \
    libcairo2-dev libxt-dev \
    libglpk-dev libgmp-dev \
    libgeos-dev libgdal-dev libudunits2-dev \
    libhdf5-dev libfontconfig1-dev zlib1g-dev \
    git \
  && rm -rf /var/lib/apt/lists/*

# ------------------------------------------------------------
# CRAN packages (binaries from P3M -> fast, no compilation)
# ------------------------------------------------------------
RUN R -e "install.packages(c( \
    'shiny','ggplot2','dplyr','tibble','tidyr','plotly','DT','qs', \
    'enrichR','patchwork','scales','stringr','ggrepel','visNetwork', \
    'remotes','BiocManager','R.utils'))"

# Seurat — separate layer for cache
RUN R -e "remotes::install_version('SeuratObject', version = '5.0.2', \
          repos = 'https://cloud.r-project.org', upgrade = 'never')"

# Seurat 5.1.0
RUN R -e "remotes::install_version('Seurat', version = '5.1.0', \
          repos = 'https://cloud.r-project.org', upgrade = 'never')"

# Bioconductor deps for hdWGCNA
RUN R -e "BiocManager::install(c('WGCNA','GeneOverlap'), ask = FALSE, update = FALSE)"

# ------------------------------------------------------------
# GitHub installs — use PAT to avoid rate limit
# ------------------------------------------------------------
ARG GITHUB_PAT
ENV GITHUB_PAT=${GITHUB_PAT}

# Pin hdWGCNA to a specific commit for reproducibility
# (replace SHA with whatever commit you've tested working)
RUN R -e "options(install.packages.check.source = 'no'); \
          remotes::install_github('smorabit/hdWGCNA', \
            ref = 'dev', \
            upgrade = 'never', \
            quiet = FALSE, \
            verbose = TRUE, \
            dependencies = TRUE); \
          if (!requireNamespace('hdWGCNA', quietly = TRUE)) \
            stop('hdWGCNA installation FAILED')"

ARG PKG_REF=main
RUN R -e "remotes::install_github('francescopatane96/scTAMsExplorer', \
          ref = '${PKG_REF}', upgrade = 'never')"

# ------------------------------------------------------------
# Entry point
# ------------------------------------------------------------
COPY run_app.R /usr/local/bin/run_app.R

ENV SEURAT_RDS=/data/atlas.rds
ENV SHINY_PORT=3838
ENV SHINY_HOST=0.0.0.0

EXPOSE 3838
CMD ["Rscript", "/usr/local/bin/run_app.R"]
