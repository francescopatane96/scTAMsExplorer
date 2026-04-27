# ============================================================
# launch_explorer.R -- Main exported function
# ============================================================

#' Launch the Seurat Atlas Explorer Shiny app
#'
#' Starts an interactive Shiny dashboard for exploring a Seurat
#' single-cell RNA-seq atlas. The dashboard provides five modules:
#' \itemize{
#'   \item \strong{UMAP + Expression} -- UMAP, FeaturePlot, VlnPlot,
#'         DotPlot, DoHeatmap with optional blend and split modes.
#'   \item \strong{DEGs + Volcano + Enrichment} -- differential expression
#'         with interactive volcano plot and EnrichR pathway enrichment.
#'   \item \strong{TF Network} -- interactive transcription factor network.
#'   \item \strong{Co-expression Modules} -- hdWGCNA hub genes and enrichment.
#'   \item \strong{TF Regulon Heatmap} -- dot heatmap of TF expression vs
#'         positive/negative regulon activity, with per-TF enrichment.
#' }
#'
#' @param seurat_obj A Seurat object. Must contain a \code{umap.harmony}
#'   reduction. The regulon heatmap tab also requires a
#'   \code{Population_level3} metadata column.
#' @param port Integer. Port to listen on (default: random available port).
#' @param launch.browser Logical. Open in the default browser? Default TRUE.
#' @param host Character. Network interface to bind to. Default
#'   \code{"127.0.0.1"} (localhost). Use \code{"0.0.0.0"} to expose on LAN.
#' @param ... Additional arguments passed to \code{\link[shiny]{runApp}}.
#'
#' @return Called for its side-effect. Returns invisibly when the app stops.
#'
#' @examples
#' \dontrun{
#' library(SeuratAtlasExplorer)
#' obj <- readRDS("my_atlas.rds")
#'
#' # Standard launch
#' launch_explorer(obj)
#'
#' # Fixed port, no browser (remote server)
#' launch_explorer(obj, port = 4242, launch.browser = FALSE)
#'
#' # Expose on local network
#' launch_explorer(obj, host = "0.0.0.0", port = 3838)
#' }
#'
#' @export
launch_explorer <- function(seurat_obj,
                             port           = NULL,
                             launch.browser = TRUE,
                             host           = "127.0.0.1",
                             ...) {
  required_pkgs <- c("qs", "shiny", "DT", "ggplot2", "dplyr", "tidyr", "tibble",
                     "plotly", "patchwork", "ggrepel", "visNetwork",
                     "stringr", "scales", "enrichR", "Seurat")
  for (pkg in required_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package '", pkg, "' is required but not installed.", call. = FALSE)
    }
    library(pkg)
    print(paste("library", pkg, "loaded"))
    if (!paste0("package:", pkg) %in% search()) {
      suppressMessages(suppressWarnings(attachNamespace(pkg)))
    }
  }

  if (!inherits(seurat_obj, "Seurat")) {
    stop("`seurat_obj` must be a Seurat object.", call. = FALSE)
  }
  if (ncol(seurat_obj) == 0L) {
    stop("`seurat_obj` contains no cells.", call. = FALSE)
  }

  missing_pkg <- character(0)
  for (pkg in c("Seurat","enrichR","visNetwork","DT","plotly","ggrepel","patchwork")) {
    if (!requireNamespace(pkg, quietly = TRUE)) missing_pkg <- c(missing_pkg, pkg)
  }
  if (length(missing_pkg) > 0L) {
    stop(
      "The following packages are required but not installed:\n  ",
      paste(missing_pkg, collapse = ", "),
      "\nInstall with: install.packages(c(",
      paste0('"', missing_pkg, '"', collapse = ", "), "))",
      call. = FALSE
    )
  }

  message("Launching Seurat Atlas Explorer ...")
  message("  Cells:    ", ncol(seurat_obj))
  message("  Features: ", nrow(seurat_obj))

  mc <- get_metadata_choices(seurat_obj)
  message("  Metadata columns: ", paste(mc, collapse = ", "))

  ui     <- atlas_ui(mc)
  server <- atlas_server(seurat_obj, mc)

  opts <- list(launch.browser = launch.browser, host = host)
  if (!is.null(port)) opts$port <- as.integer(port)

  shiny::shinyApp(ui = ui, server = server, options = opts)
}
