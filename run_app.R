
#!/usr/bin/env Rscript

qs_path <- Sys.getenv("SEURAT_QS", unset = "/data/atlas.qs")
port     <- as.integer(Sys.getenv("SHINY_PORT", unset = "3838"))
host     <- Sys.getenv("SHINY_HOST", unset = "0.0.0.0")

cat("scTAMsExplorer — containerised launcher\n")
cat("  QS path:", qs_path, "\n")
cat("  Host    :", host,     "\n")
cat("  Port    :", port,     "\n")

if (!file.exists(qs_path)) {
  stop("Seurat .qs not found at: ", qs_path,
       "\nMount it with -v /host/path.qs:", qs_path, ":ro",
       call. = FALSE)
}

suppressPackageStartupMessages({
  library(scTAMsExplorer)
  library(Seurat)
  library(qs)
})

cat("Loading Seurat object...\n")
seurat_obj <- qread(qs_path)

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
