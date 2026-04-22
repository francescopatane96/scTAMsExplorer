# ============================================================
# launch_explorer.R — Main user-facing function
# ============================================================

#' Launch the Seurat Atlas Explorer Shiny app
#'
#' Starts an interactive Shiny dashboard for exploring a Seurat single-cell
#' RNA-seq atlas. The dashboard provides five analytical modules:
#' \itemize{
#'   \item \strong{UMAP + Expression} — UMAP, FeaturePlot, VlnPlot, DotPlot,
#'         DoHeatmap with optional blend and split modes.
#'   \item \strong{DEGs + Volcano + Enrichment} — differential expression
#'         between any two groups with interactive volcano plot and EnrichR
#'         pathway enrichment.
#'   \item \strong{TF Network} — interactive transcription factor interaction
#'         network (visNetwork) with multiple layout algorithms.
#'   \item \strong{Co-expression Modules} — hdWGCNA module hub genes, module
#'         sizes, gene-to-module lookup, and per-module enrichment.
#'   \item \strong{TF Regulon Heatmap} — dot heatmap of TF expression vs
#'         positive/negative regulon activity across clusters, with per-TF
#'         EnrichR enrichment.
#' }
#'
#' @param seurat_obj A Seurat object. Must contain a \code{umap.harmony}
#'   reduction and a \code{Population_level3} metadata column for the
#'   regulon heatmap module. All other modules are flexible.
#' @param port       Integer. Port to run the app on. Defaults to a random
#'   available port chosen by Shiny.
#' @param launch.browser Logical. Open the app in the default browser
#'   automatically? Default \code{TRUE}.
#' @param host       Character. Network interface to listen on.
#'   Default \code{"127.0.0.1"} (localhost only). Set to \code{"0.0.0.0"}
#'   to expose on the local network.
#' @param ...        Additional arguments passed to \code{\link[shiny]{runApp}}.
#'
#' @return Called for its side-effect (launches a Shiny app). Returns
#'   invisibly when the app is stopped.
#'
#' @examples
#' \dontrun{
#' # Load your Seurat object then launch the explorer
#' library(SeuratAtlasExplorer)
#' seurat_obj <- readRDS("my_atlas.rds")
#' launch_explorer(seurat_obj)
#'
#' # Use a fixed port and keep the browser closed (e.g. in a server session)
#' launch_explorer(seurat_obj, port = 4242, launch.browser = FALSE)
#' }
#'
#' @export
launch_explorer <- function(seurat_obj,
                             port           = NULL,
                             launch.browser = TRUE,
                             host           = "127.0.0.1",
                             ...) {

  # ---- Input validation -------------------------------------------
  if (!inherits(seurat_obj, "Seurat")) {
    stop("`seurat_obj` must be a Seurat object.", call. = FALSE)
  }

  if (ncol(seurat_obj) == 0L) {
    stop("`seurat_obj` contains no cells.", call. = FALSE)
  }

  # ---- Optional dependency checks ---------------------------------
  missing_pkg <- character(0)
  for (pkg in c("Seurat", "enrichR", "visNetwork", "DT", "plotly",
                "ggrepel", "patchwork")) {
    if (!requireNamespace(pkg, quietly = TRUE)) missing_pkg <- c(missing_pkg, pkg)
  }
  if (length(missing_pkg)) {
    stop("The following packages are required but not installed:\n  ",
         paste(missing_pkg, collapse = ", "),
         "\nInstall them with: install.packages(c(",
         paste0('"', missing_pkg, '"', collapse = ", "), "))",
         call. = FALSE)
  }

  message("Launching Seurat Atlas Explorer...")
  message("  Cells:    ", ncol(seurat_obj))
  message("  Features: ", nrow(seurat_obj))
  message("  Metadata: ", paste(get_metadata_choices(seurat_obj), collapse = ", "))

  # ---- Build UI and server ----------------------------------------
  metadata_choices <- get_metadata_choices(seurat_obj)
  ui     <- atlas_ui(metadata_choices)
  server <- atlas_server(seurat_obj, metadata_choices)

  # ---- Run --------------------------------------------------------
  options_list <- list(launch.browser = launch.browser, host = host)
  if (!is.null(port)) options_list$port <- as.integer(port)

  shiny::shinyApp(
    ui      = ui,
    server  = server,
    options = options_list
  )
}
