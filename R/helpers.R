# ============================================================
# helpers.R -- Shared utilities for SeuratAtlasExplorer
# ============================================================

#' @keywords internal
GetModules <- function(obj, wgcna_name=NULL){
  
  # get data from active assay if wgcna_name is not given
  if(is.null(wgcna_name)){wgcna_name <- obj@misc$active_wgcna}
  
  obj@misc[[wgcna_name]]$wgcna_modules
}



#' @keywords internal
GetHubGenes <- function(
    obj,
    n_hubs = 10,
    mods = NULL,
    wgcna_name=NULL
){
  
  if(is.null(wgcna_name)){wgcna_name <- obj@misc$active_wgcna}
  
  # get the modules table
  modules <- GetModules(obj, wgcna_name) %>% subset(module != 'grey')
  
  if(is.null(mods)){
    mods <- levels(modules$module); mods <- mods[mods != 'grey']
  } else{
    if(!all(mods %in% modules$module)){
      stop("Invalid selection for mods.")
    }
  }
  
  #get hub genes:
  hub_df <- do.call(rbind, lapply(mods, function(cur_mod){
    cur <- subset(modules, module == cur_mod)
    cur <- cur[,c('gene_name', 'module', paste0('kME_', cur_mod))]
    names(cur)[3] <- 'kME'
    cur <- dplyr::arrange(cur, desc(kME))
    cur %>% dplyr::slice_max(n=n_hubs, order_by=kME)
  }))
  rownames(hub_df) <- 1:nrow(hub_df)
  hub_df
  
}

#' Extract the TF regulatory network from a Seurat object
#'
#' Retrieves the transcription factor network stored in the misc
#' slot of a Seurat object by hdWGCNA. If `wgcna_name` is NULL,
#' the function uses the currently active hdWGCNA experiment
#' (`obj@misc$active_wgcna`).
#'
#' @param obj A Seurat object with an hdWGCNA experiment that
#'   includes a `tf_net` element.
#' @param wgcna_name Character. Name of the hdWGCNA experiment
#'   to query. Defaults to the active one.
#'
#' @return A data frame describing TF -> target relationships
#'   (typically with columns `tf`, `gene`, `Gain`, `Cor`).
#'
#' @keywords internal
GetTFNetwork <- function(obj, wgcna_name = NULL) {
  if (is.null(wgcna_name)) {
    wgcna_name <- obj@misc$active_wgcna
  }
  obj@misc[[wgcna_name]]$tf_net
}

#' Extract suitable metadata columns from a Seurat object
#'
#' Returns the names of metadata columns that can be used as
#' grouping variables in the dashboard. A column is "suitable"
#' if it is a factor or a character vector.
#'
#' @param obj A Seurat object.
#'
#' @return Character vector of column names.
#'
#' @keywords internal
get_metadata_choices <- function(obj) {
  meta     <- obj@meta.data
  suitable <- vapply(meta, function(x) is.factor(x) || is.character(x), logical(1))
  names(meta)[suitable]
}

#' Build a sidebar block of plot-size controls
#'
#' Generates a tagList with paired display-size and
#' download-size inputs (width, height, DPI, font size). Used
#' across all five tabs to keep plot-control UI consistent.
#'
#' @param id_prefix Character prefix for the input IDs. The
#'   resulting inputs are named `<prefix>_width`, `<prefix>_height`,
#'   `<prefix>_dl_w`, `<prefix>_dl_h`, `<prefix>_dpi`, `<prefix>_pt`.
#' @param default_w Default display width in pixels.
#' @param default_h Default display height in pixels.
#' @param default_pt Default font size in points.
#'
#' @return A `shiny::tagList()` of UI elements.
#'
#' @importFrom shiny tagList tags fluidRow column numericInput
#'
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
#'
#' Wraps `enrichR::enrichr()` to query a single database and
#' return the result as a tidy data frame, with a `neg_log10_p`
#' column added and term names truncated for plot readability.
#' Adds a 1-second sleep before the API call to be polite.
#'
#' @param genes Character vector of gene symbols.
#' @param db Character. Name of an EnrichR database (e.g.
#'   `"GO_Biological_Process_2023"`).
#'
#' @return A data frame with one row per enriched term, or
#'   `NULL` if the query returned nothing.
#'
#' @importFrom enrichR enrichr
#' @importFrom dplyr arrange mutate
#' @importFrom rlang .data
#' @importFrom stringr str_trunc
#'
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
#'
#' Renders a horizontal bar chart of the top enriched terms,
#' sorted by Combined Score and coloured by significance.
#'
#' @param df Tidy EnrichR data frame, typically the output of
#'   [run_enrichr_tidy()].
#' @param color_high Character. The high colour for the fill
#'   gradient (e.g. `"firebrick"`).
#' @param title_txt Character. Plot title.
#' @param n_terms Integer. Number of top terms to plot.
#' @param pt Numeric. Base font size for the theme.
#'
#' @return A ggplot object.
#'
#' @importFrom utils head
#' @importFrom rlang .data
#' @importFrom ggplot2 ggplot aes geom_col scale_fill_gradient
#' @importFrom ggplot2 coord_flip theme_minimal labs theme element_text
#'
#' @keywords internal
make_enrich_plot <- function(df, color_high, title_txt, n_terms, pt) {
  df_top <- utils::head(df, n_terms)
  ggplot2::ggplot(
    df_top,
    ggplot2::aes(
      x    = stats::reorder(.data$Term, .data$Combined.Score),
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
