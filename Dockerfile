FROM rocker/r-ver:4.4.1

# ------------------------------------------------------------
# System libraries needed by R deps:
#   - Seurat / SeuratObject  -> HDF5, BLAS, igraph (glpk, gmp)
#   - visNetwork / plotly    -> standard web deps
#   - ggplot2 / patchwork    -> Cairo, fonts
#   - hdWGCNA (Suggests)     -> GEOS, GDAL, udunits
# ------------------------------------------------------------
RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libpng-dev \
    libjpeg-dev \
    libtiff5-dev \
    libfreetype6-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libcairo2-dev \
    libxt-dev \
    libglpk-dev \
    libgmp-dev \
    libgeos-dev \
    libgdal-dev \
    libudunits2-dev \
    libhdf5-dev \
    libfontconfig1-dev \
    zlib1g-dev \
    git \
  && rm -rf /var/lib/apt/lists/*

# ------------------------------------------------------------
# Install CRAN packages needed by scTAMsExplorer
# ------------------------------------------------------------
RUN R -e "install.packages(c( \
    'shiny', \
    'ggplot2', \
    'dplyr', \
    'tibble', \
    'tidyr', \
    'plotly', \
    'DT', \
    'qs', \
    'enrichR', \
    'patchwork', \
    'scales', \
    'stringr', \
    'ggrepel', \
    'visNetwork', \
    'remotes', \
    'BiocManager', \
    'R.utils' \
  ), repos = 'https://cloud.r-project.org')"

# Seurat is heavy — separate layer for better caching
RUN R -e "install.packages('Seurat', repos = 'https://cloud.r-project.org')"

# Bioconductor dependencies used (indirectly) by hdWGCNA
RUN R -e "BiocManager::install(c('WGCNA','GeneOverlap'), ask = FALSE, update = FALSE)"

# hdWGCNA — required by the Modules tab of the app
RUN R -e "remotes::install_github('smorabit/hdWGCNA', ref = 'dev', upgrade = 'never')"

# ------------------------------------------------------------
# Install scTAMsExplorer from GitHub
# ------------------------------------------------------------
ARG PKG_REF=main
RUN R -e "remotes::install_github('francescopatane96/scTAMsExplorer', \
         ref = '${PKG_REF}', upgrade = 'never')"

# ------------------------------------------------------------
# Entry point
# ------------------------------------------------------------
COPY run_app.R /usr/local/bin/run_app.R
EXPOSE 3838

ENV SEURAT_RDS=/data/atlas.rds
ENV SHINY_PORT=3838
ENV SHINY_HOST=0.0.0.0

CMD ["Rscript", "/usr/local/bin/run_app.R"]
