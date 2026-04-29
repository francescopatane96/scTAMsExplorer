# ============================================================
# scTAMsExplorer-package.R
# Centralised package-level documentation and imports.
#
# This file is the conventional place to declare imports that
# would be tedious to repeat on every function (e.g. all of
# ggplot2 or shiny). Function-specific imports stay near each
# function via @importFrom tags.
# ============================================================

#' scTAMsExplorer: Interactive Single-Cell RNA-seq Atlas Explorer
#'
#' A modular Shiny dashboard for exploring single-cell RNA-seq
#' atlases built with Seurat. Provides interactive tools for UMAP
#' visualisation, differential gene expression, volcano plots,
#' pathway enrichment (EnrichR), transcription factor network
#' analysis, co-expression module exploration (hdWGCNA), and
#' TF regulon heatmaps.
#'
#' The main entry point is [launch_explorer()].
#'
#' @keywords internal
#' @name scTAMsExplorer-package
#' @aliases scTAMsExplorer
#'
#' ---------------------------------------------------------
#' Package-level imports ("in blocco")
#' ---------------------------------------------------------
#'
#' shiny is imported in full because the server function uses
#' dozens of its symbols. We exclude `dataTableOutput` and
#' `renderDataTable`: those names also exist in DT and we want
#' the DT versions (richer features, used elsewhere in the app).
#'
#' @rawNamespace import(shiny, except = c(dataTableOutput, renderDataTable))
#'
#' DT is imported in full so that `dataTableOutput` and
#' `renderDataTable` resolve to DT's versions, and `datatable()`
#' / `formatSignif()` are available unprefixed.
#'
#' @import DT
#'
#' ggplot2 is imported in full because virtually every plot in
#' the app uses several ggplot2 symbols (ggplot, aes, geom_*,
#' theme_*, scale_*, labs, etc.). Importing piecemeal would
#' bloat every roxygen header.
#'
#' @import ggplot2
#'
#' dplyr/tidyr verbs are used pervasively in reactive pipelines.
#'
#' @import dplyr
#' @import tidyr
#'
#' The pipe operator from magrittr is used throughout
#' (alongside R 4.1's native `|>`).
#'
#' @importFrom magrittr %>%
#'
#' Common base/utils helpers used across reactives.
#'
#' @importFrom utils head write.csv
#' @importFrom stats reorder
#'
"_PACKAGE"
