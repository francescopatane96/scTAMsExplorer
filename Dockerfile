FROM rocker/r-ver:4.4.1

# System libraries needed by your R deps (adatta alle tue)
RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev libssl-dev libxml2-dev \
    libpng-dev libjpeg-dev libfreetype6-dev \
    libcairo2-dev libharfbuzz-dev libfribidi-dev \
    git \
  && rm -rf /var/lib/apt/lists/*

# Install your package's R dependencies first (better cache)
RUN R -e "install.packages(c('shiny','remotes','BiocManager'), \
         repos='https://cloud.r-project.org')"

# Install your package from GitHub
ARG PKG_REF=main
RUN R -e "remotes::install_github('francescopatane96/scTAMsExplorer', \
         ref='${PKG_REF}', upgrade='never')"

# Entry point
COPY run_app.R /usr/local/bin/run_app.R
EXPOSE 3838
CMD ["Rscript", "/usr/local/bin/run_app.R"]
