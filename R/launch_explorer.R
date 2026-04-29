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
#' @param seurat_obj A Seurat object. Must contain a `umap.harmony`
#'   reduction. The regulon heatmap tab also requires a
#'   `Population_level3` metadata column.
#' @param port Integer. Port to listen on. Default `NULL` lets
#'   Shiny pick a random available port.
#' @param launch.browser Logical. Open the app in the default
#'   browser? Default `TRUE`.
#' @param host Character. Network interface to bind to.
#'   Default `"127.0.0.1"` (localhost only). Use `"0.0.0.0"` to
#'   expose the app on the local network.
#' @param ... Additional arguments forwarded to
#'   [shiny::runApp()] via the `options` argument of
#'   [shiny::shinyApp()].
#'
#' @return Called for its side-effect (starts the Shiny app).
#'   Returns the shinyApp object invisibly.
#'
#' @examples
#' \dontrun{
#' library(scTAMsExplorer)
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
#' @importFrom shiny shinyApp
#'
#' @export
launch_explorer <- function(seurat_obj,
                            port           = NULL,
                            launch.browser = TRUE,
                            host           = "127.0.0.1",
                            ...) {

  # ---- Input validation ----------------------------------
  if (!inherits(seurat_obj, "Seurat")) {
    stop("`seurat_obj` must be a Seurat object.", call. = FALSE)
  }
  if (ncol(seurat_obj) == 0L) {
    stop("`seurat_obj` contains no cells.", call. = FALSE)
  }

  # ---- Startup banner ------------------------------------
  message("Launching Seurat Atlas Explorer ...")
  message("  Cells:    ", ncol(seurat_obj))
  message("  Features: ", nrow(seurat_obj))

  mc <- get_metadata_choices(seurat_obj)
  message("  Metadata columns: ", paste(mc, collapse = ", "))

  # ---- Build app ----------------------------------------
  ui     <- atlas_ui(mc)
  server <- atlas_server(seurat_obj, mc)

  opts <- list(launch.browser = launch.browser, host = host)
  if (!is.null(port)) opts$port <- as.integer(port)

  shiny::shinyApp(ui = ui, server = server, options = opts)
}
