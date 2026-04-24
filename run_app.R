
#!/usr/bin/env Rscript

rds_path <- Sys.getenv("SEURAT_RDS", unset = "/data/atlas.rds")
port     <- as.integer(Sys.getenv("SHINY_PORT", unset = "3838"))
host     <- Sys.getenv("SHINY_HOST", unset = "0.0.0.0")

cat("scTAMsExplorer — containerised launcher\n")
cat("  RDS path:", rds_path, "\n")
cat("  Host    :", host,     "\n")
cat("  Port    :", port,     "\n")

if (!file.exists(rds_path)) {
  stop("Seurat .rds not found at: ", rds_path,
       "\nMount it with -v /host/path.rds:", rds_path, ":ro",
       call. = FALSE)
}

suppressPackageStartupMessages({
  library(scTAMsExplorer)
  library(Seurat)
})

cat("Loading Seurat object...\n")
seurat_obj <- readRDS(rds_path)

# Alias umap 
if (!"umap.harmony" %in% names(seurat_obj@reductions) &&
     "umap"         %in% names(seurat_obj@reductions)) {
  seurat_obj[["umap.harmony"]] <- seurat_obj[["umap"]]
}

scTAMsExplorer::launch_explorer(
  seurat_obj      = seurat_obj,     
  port            = port,
  host            = host,
  launch.browser  = FALSE
)
