# 🧬 scTAMsExplorer

> An interactive Shiny dashboard for exploring a single cell transcriptomic atlas of Mono&Macro extracted from humanized MISTRG healthy and tumoral tissues.  — developed by the **Goriely Lab**, Institute of Medical Immunology (IMI), University of Brussels (ULB)

---

## Features

| Tab | What you can do |
|-----|----------------|
| 🗺️ **UMAP + Expression** | UMAP by any metadata · FeaturePlot · VlnPlot · DotPlot · Heatmap · blend & split modes · marker genes |
| 🔬 **DEGs + Volcano** | FindMarkers between any two groups · interactive volcano (ggrepel) · separate EnrichR enrichment step · geneset tables |
| 🌐 **TF Network** | Interactive visNetwork graph · 5 layouts · Gain/top-N filters · reverse search · TF–TF only mode · CSV downloads |
| 🔗 **Co-expression Modules** | hdWGCNA hub genes · module sizes · gene → module lookup · per-module EnrichR enrichment |
| 🎯 **TF Regulon Heatmap** | Dot heatmap of TF expression vs positive/negative regulon activity · module filter · per-TF EnrichR enrichment |

All plots have adjustable **display size**, **download resolution (DPI + inches)**, and **font size**. All tables are sortable, filterable, and paginated via DataTables.

---
## Installation - Docker (recommended way)

### 1 · Install Docker Desktop
### 2 · Change the docker memory limit up to 16 GB
Settings/Resources/Resource Allocation/Memory Limit

### 3 · Download the image and run the container via terminal (adapt "/path/to/your/seurat_object.rds")

```
docker run --rm -p 3838:3838 \
     -v /path/to/your/seurat_object.rds:/data/atlas.rds:ro \
     ghcr.io/francescopatane96/sctamsexplorer:latest
```
wait until the seurat object has been completely loaded

### 3 · open a browser and go to http://localhost:3838 

## Installation

### 1 · Install from GitHub

```r
# Install remotes if needed
install.packages("remotes")

# Install SeuratAtlasExplorer
remotes::install_github("francescopatane96/scTAMsExplorer")
```

### 2 · Install required dependencies

Most dependencies are installed automatically. If any are missing, install them manually:

```r
# CRAN packages
install.packages(c(
  "shiny", "ggplot2", "dplyr", "tidyr", "tibble",
  "plotly", "DT", "enrichR", "patchwork", "scales",
  "stringr", "ggrepel", "visNetwork"
))

# Bioconductor / GitHub packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Seurat")   # if not already installed

# hdWGCNA (required for co-expression module tab)
remotes::install_github("smorabit/hdWGCNA", ref = "dev")
```

### 3 · Import required dependencies

```
library(shiny)
library(Seurat)
library(ggplot2)
library(dplyr)
library(tibble)
library(tidyr)
library(plotly)
library(DT)
library(enrichR)
library(patchwork)
library(scales)
library(stringr)
library(ggrepel)
library(visNetwork)
```

---

## Quick start

```r
library(scTAMsExplorer)

# Load your Seurat object
seurat_obj <- readRDS("path/to/my_atlas.rds")

# Launch the app — opens automatically in your default browser
launch_explorer(seurat_obj)
```

That's it. The app opens in your browser and all analyses run on your object.

---

## Function reference

### `launch_explorer()`

```r
launch_explorer(
  seurat_obj,            # Your Seurat object (required)
  port           = NULL, # Fixed port, e.g. 4242 (optional)
  launch.browser = TRUE, # Open browser automatically?
  host           = "127.0.0.1", # "0.0.0.0" to expose on LAN
  ...                    # Passed to shiny::runApp()
)
```

#### Arguments

| Argument | Type | Default | Description |
|----------|------|---------|-------------|
| `seurat_obj` | Seurat | — | **Required.** Your processed Seurat atlas object. |
| `port` | integer | `NULL` | Fixed port number. Leave `NULL` for automatic selection. |
| `launch.browser` | logical | `TRUE` | Open the app in the default browser? |
| `host` | character | `"127.0.0.1"` | Set to `"0.0.0.0"` to make accessible on your local network (e.g., from another machine). |
| `...` | — | — | Additional arguments forwarded to `shiny::runApp()`. |

#### Examples

```r
# Standard usage
launch_explorer(seurat_obj)

# Fixed port, no browser (e.g. running on a remote server)
launch_explorer(seurat_obj, port = 8080, launch.browser = FALSE)

# Expose on local network (accessible from other devices)
launch_explorer(seurat_obj, host = "0.0.0.0", port = 3838)
```

---

## Seurat object requirements

The app is designed to work with **any Seurat object**, but some tabs assume specific slots:

| Tab | Requirement |
|-----|-------------|
| All tabs | Any `factor` or `character` column in `@meta.data` is available as a grouping variable |
| UMAP | A reduction named `umap.harmony` (rename yours if needed: `seurat_obj[["umap.harmony"]] <- seurat_obj[["umap"]]`) |
| TF Network | A `GetTFNetwork()` function available (requires a compatible GRN object) |
| Co-expression | hdWGCNA results stored in the Seurat object (`GetModules()`, `GetHubGenes()`) |
| TF Regulon Heatmap | A `Population_level3` metadata column for `AverageExpression()` grouping (configurable in the source) |

### Rename your UMAP reduction (if needed)

```r
# If your UMAP is stored as "umap" rather than "umap.harmony":
seurat_obj[["umap.harmony"]] <- seurat_obj[["umap"]]
```

---

## Running on a server (e.g. RStudio Server / HPC)

```r
# On a server, disable browser auto-launch and bind to all interfaces
launch_explorer(
  seurat_obj,
  port           = 3838,
  launch.browser = FALSE,
  host           = "0.0.0.0"
)
# Then open http://<server-ip>:3838 in your browser
```

For **RStudio Server**, simply run `launch_explorer(seurat_obj)` and RStudio will open the app in the Viewer pane or a new browser tab automatically.

---

## Citation

If you use this app in your research, please cite:

> Detavernier, Donckier de Donceel, Patanè, Pedron et al. Regulatory Networks Shaping Human Tumor-Associated Macrophages In Vivo Identify MAFB as a Key Transcriptional Driver, 18 January 2026, PREPRINT (Version 1) available at Research Square [https://doi.org/10.21203/rs.3.rs-8562817/v1]

> GitHub: https://github.com/francescopatane96/SeuratAtlasExplorer

---

## License

MIT © Goriely Lab, University of Brussels
