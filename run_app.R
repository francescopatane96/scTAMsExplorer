
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
#
required_pkgs <- c("qs", "shiny", "DT", "ggplot2", "dplyr", "tidyr", "tibble",
                     "plotly", "patchwork", "ggrepel", "visNetwork",
                     "stringr", "scales", "enrichR", "Seurat", "scTAMsExplorer")
  for (pkg in required_pkgs) {
    print("ciao")
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package '", pkg, "' is required but not installed.", call. = FALSE)
    }
    
    if (!paste0("package:", pkg) %in% search()) {
      suppressMessages(suppressWarnings(attachNamespace(pkg)))
    }
  }

cat("Loading Seurat object...\n")
seurat_obj <- qread(qs_path, nthreads=1)

# Alias umap 
if (!"umap.harmony" %in% names(seurat_obj@reductions) &&
     "umap"         %in% names(seurat_obj@reductions)) {
  seurat_obj[["umap.harmony"]] <- seurat_obj[["umap"]]
}

print("ready to launch")

scTAMsExplorer::launch_explorer(
  seurat_obj      = seurat_obj,     
  port            = port,
  host            = host,
  launch.browser  = FALSE
)
