# ============================================================
# helpers.R -- Shared utilities for SeuratAtlasExplorer
# ============================================================

GetTFNetwork <- function(obj, wgcna_name=NULL){
  if(is.null(wgcna_name)){wgcna_name <- obj@misc$active_wgcna}
  #CheckWGCNAName(seurat_obj, wgcna_name)
  obj@misc[[wgcna_name]]$tf_net 
}

#' Extract suitable metadata columns from a Seurat object
#' @param obj A Seurat object.
#' @return Character vector of column names.
#' @keywords internal
get_metadata_choices <- function(obj) {
  meta     <- obj@meta.data
  suitable <- vapply(meta, function(x) is.factor(x) || is.character(x), logical(1))
  names(meta)[suitable]
}

#' Sidebar plot-size control block
#' @param id_prefix    Character prefix for input IDs.
#' @param default_w,default_h Default display width/height in pixels.
#' @param default_pt   Default font size in points.
#' @return A tagList of Shiny UI elements.
#' @keywords internal
plot_size_controls <- function(id_prefix,
                               default_w  = 800,
                               default_h  = 560,
                               default_pt = 12) {
  shiny::tagList(
    shiny::tags$div(class = "ctrl-section",
      shiny::tags$p(class = "ctrl-label", "Display size (px)"),
      shiny::fluidRow(
        shiny::column(6, shiny::numericInput(
          paste0(id_prefix, "_width"),  "Width",
          value = default_w, min = 200, step = 50)),
        shiny::column(6, shiny::numericInput(
          paste0(id_prefix, "_height"), "Height",
          value = default_h, min = 200, step = 50))
      ),
      shiny::tags$p(class = "ctrl-label", "Download (inches / DPI)"),
      shiny::fluidRow(
        shiny::column(4, shiny::numericInput(
          paste0(id_prefix, "_dl_w"), "W",   value = 10,  min = 3, step = 1)),
        shiny::column(4, shiny::numericInput(
          paste0(id_prefix, "_dl_h"), "H",   value = 8,   min = 3, step = 1)),
        shiny::column(4, shiny::numericInput(
          paste0(id_prefix, "_dpi"),  "DPI", value = 300, min = 72, step = 50))
      ),
      shiny::numericInput(paste0(id_prefix, "_pt"), "Font size (pt)",
                          value = default_pt, min = 6, max = 30)
    )
  )
}

#' Run EnrichR and return a tidy data frame
#' @param genes Character vector of gene symbols.
#' @param db    EnrichR database name.
#' @return A data frame or NULL.
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

#' Build an EnrichR bar plot
#' @param df         Tidy EnrichR data frame.
#' @param color_high High colour for fill gradient.
#' @param title_txt  Plot title.
#' @param n_terms    Number of terms.
#' @param pt         Base font size.
#' @return A ggplot object.
#' @keywords internal
make_enrich_plot <- function(df, color_high, title_txt, n_terms, pt) {
  df_top <- utils::head(df, n_terms)
  ggplot2::ggplot(
    df_top,
    ggplot2::aes(
      x    = reorder(.data$Term, .data$Combined.Score),
      y    = .data$Combined.Score,
      fill = -log10(.data$Adjusted.P.value + 1e-300)
    )
  ) +
    ggplot2::geom_col(color = "white", linewidth = 0.25) +
    ggplot2::scale_fill_gradient(low = "lightyellow", high = color_high,
                                 name = "-log10\n(adj.p)") +
    ggplot2::coord_flip() +
    ggplot2::theme_minimal(base_size = pt) +
    ggplot2::labs(title = title_txt, x = NULL, y = "Combined Score") +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
}
