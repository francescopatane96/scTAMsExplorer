# ============================================================
# helpers.R — Shared utilities for SeuratAtlasExplorer
# ============================================================

#' Extract suitable metadata columns from a Seurat object
#'
#' Returns the names of columns in \code{seurat_obj@meta.data} that are
#' either factors or character vectors, making them suitable for grouping /
#' colouring cells in UMAP and expression plots.
#'
#' @param obj A Seurat object.
#' @return A character vector of column names.
#' @keywords internal
get_metadata_choices <- function(obj) {
  meta     <- obj@meta.data
  suitable <- vapply(meta, function(x) is.factor(x) || is.character(x), logical(1))
  names(meta)[suitable]
}

# ---- Sidebar plot-size control block -----------------------------------

#' Build a reusable plot-size control panel for the sidebar
#'
#' Produces a set of \code{numericInput} widgets for controlling the
#' on-screen pixel dimensions of a plot, its download dimensions (in
#' inches), DPI, and font size.
#'
#' @param id_prefix Character. Prefix appended to each input id.
#' @param default_w,default_h Default display width / height in pixels.
#' @param default_pt Default font size in points.
#' @return A \code{tagList} of Shiny UI elements.
#' @keywords internal
plot_size_controls <- function(id_prefix,
                               default_w  = 800,
                               default_h  = 560,
                               default_pt = 12) {
  shiny::tagList(
    shiny::tags$div(class = "ctrl-section",
      shiny::tags$p(class = "ctrl-label", "Display size (px)"),
      shiny::fluidRow(
        shiny::column(6, shiny::numericInput(paste0(id_prefix, "_width"),
                                             "Width",  value = default_w, min = 200, step = 50)),
        shiny::column(6, shiny::numericInput(paste0(id_prefix, "_height"),
                                             "Height", value = default_h, min = 200, step = 50))
      ),
      shiny::tags$p(class = "ctrl-label", "Download (inches \u00b7 DPI)"),
      shiny::fluidRow(
        shiny::column(4, shiny::numericInput(paste0(id_prefix, "_dl_w"), "W",
                                             value = 10, min = 3, step = 1)),
        shiny::column(4, shiny::numericInput(paste0(id_prefix, "_dl_h"), "H",
                                             value = 8,  min = 3, step = 1)),
        shiny::column(4, shiny::numericInput(paste0(id_prefix, "_dpi"),  "DPI",
                                             value = 300, min = 72, step = 50))
      ),
      shiny::numericInput(paste0(id_prefix, "_pt"), "Font size (pt)",
                          value = default_pt, min = 6, max = 30)
    )
  )
}

# ---- EnrichR enrichment helper -----------------------------------------

#' Run EnrichR and return a tidy data frame
#'
#' Calls \code{enrichR::enrichr}, waits 1 second to avoid API rate limits,
#' and returns a tidy result sorted by adjusted p-value with a
#' \code{neg_log10_p} column added and Term strings truncated to 65 chars.
#'
#' @param genes Character vector of gene symbols.
#' @param db    EnrichR database name.
#' @return A data frame, or \code{NULL} if no results.
#' @keywords internal
run_enrichr_tidy <- function(genes, db) {
  Sys.sleep(1)
  res <- enrichR::enrichr(genes, databases = db)
  df  <- res[[db]]
  if (is.null(df) || nrow(df) == 0L) return(NULL)
  df %>%
    dplyr::arrange(.data$Adjusted.P.value) %>%
    dplyr::mutate(
      neg_log10_p = -log10(.data$Adjusted.P.value + 1e-300),
      Term        = stringr::str_trunc(.data$Term, 65)
    )
}

# ---- Enrichment bar-plot helper ----------------------------------------

#' Build an EnrichR bar plot
#'
#' @param df         Tidy EnrichR data frame (from \code{run_enrichr_tidy}).
#' @param color_high High colour for the fill gradient.
#' @param title_txt  Plot title string.
#' @param n_terms    Number of terms to show.
#' @param pt         Base font size (pt).
#' @return A \code{ggplot} object.
#' @keywords internal
make_enrich_plot <- function(df, color_high, title_txt, n_terms, pt) {
  df_top <- utils::head(df, n_terms)
  ggplot2::ggplot(df_top,
                  ggplot2::aes(x = reorder(.data$Term, .data$Combined.Score),
                               y = .data$Combined.Score,
                               fill = -log10(.data$Adjusted.P.value + 1e-300))) +
    ggplot2::geom_col(color = "white", linewidth = 0.25) +
    ggplot2::scale_fill_gradient(low = "lightyellow", high = color_high,
                                 name = "-log10\n(adj.p)") +
    ggplot2::coord_flip() +
    ggplot2::theme_minimal(base_size = pt) +
    ggplot2::labs(title = title_txt, x = NULL, y = "Combined Score") +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
}
