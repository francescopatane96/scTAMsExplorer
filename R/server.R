# ============================================================
# server.R -- Server builder for SeuratAtlasExplorer
# ============================================================

#' Build the Shiny server for the Atlas Explorer
#'
#' Returns a standard Shiny server function closed over
#' `seurat_obj` and `metadata_choices`, so no global variables
#' are required. The returned closure wires up reactives,
#' outputs and download handlers for all five tabs of the
#' dashboard:
#'
#' \enumerate{
#'   \item UMAP + Expression (DimPlot, FeaturePlot, VlnPlot,
#'         DotPlot, DoHeatmap, marker discovery).
#'   \item DEGs + Volcano + Enrichment.
#'   \item TF interaction network (visNetwork).
#'   \item hdWGCNA co-expression modules.
#'   \item TF regulon dot heatmap with per-TF enrichment.
#' }
#'
#' This function is not exported. End users start the app via
#' [launch_explorer()], which wires this server together with
#' [atlas_ui()] inside a [shiny::shinyApp()] call.
#'
#' @param seurat_obj A Seurat object. Must contain a
#'   `umap.harmony` reduction. The regulon heatmap tab also
#'   requires a `Population_level3` metadata column and an
#'   hdWGCNA experiment with a TF network.
#' @param metadata_choices Character vector of metadata column
#'   names, typically from [get_metadata_choices()]. Used to
#'   resolve user-selected grouping variables.
#'
#' @return A Shiny server function `function(input, output, session)`.
#'
#' @section Imports:
#' Almost the entire ggplot2 / dplyr / shiny / DT / plotly /
#' visNetwork surface is used here. To keep this header
#' readable, broad imports are declared once at the package
#' level in `R/scTAMsExplorer-package.R` (shiny, DT, ggplot2,
#' dplyr, tidyr, magrittr). Below we only list function-specific
#' imports that are NOT covered by those package-level imports.
#'
#' @importFrom Seurat DimPlot FeaturePlot VlnPlot DotPlot DoHeatmap
#' @importFrom Seurat FindMarkers AverageExpression
#' @importFrom Seurat Idents DefaultAssay
#' @importFrom SeuratObject Idents<-
#' @importFrom tibble rownames_to_column
#' @importFrom plotly ggplotly renderPlotly
#' @importFrom visNetwork visNetwork visEdges visNodes visOptions
#' @importFrom visNetwork visInteraction visPhysics visLayout
#' @importFrom visNetwork renderVisNetwork
#' @importFrom ggrepel geom_text_repel
#' @importFrom patchwork plot_annotation
#' @importFrom scales rescale
#' @importFrom stringr str_trunc
#' @importFrom enrichR enrichr
#' @importFrom grid unit
#' @importFrom utils write.csv head
#' @importFrom stats reorder
#'
#' @keywords internal
atlas_server <- function(seurat_obj, metadata_choices) {
  function(input, output, session) {

  # ================================================
  # TAB 1
  # ================================================
    library(patchwork)

  output$cluster_selector <- renderUI({
    req(input$umap_groupby)
    choices <- sort(unique(as.character(seurat_obj@meta.data[[input$umap_groupby]])))
    selectInput("cluster", paste0("Cluster (", input$umap_groupby, "):"),
                choices = choices, selected = choices[1])
  })

  umap_plot_obj <- reactive({
    req(input$umap_groupby)
    DimPlot(seurat_obj, reduction = "umap.harmony",
            group.by = input$umap_groupby, label = TRUE,
            label.size = input$umap_pt / 3) +
      theme_minimal(base_size = input$umap_pt)
  })

  output$umap_container <- renderUI({
    req(input$umap_width, input$umap_height)
    plotOutput("umap", width = paste0(input$umap_width, "px"),
                       height = paste0(input$umap_height, "px"))
  })
  output$umap <- renderPlot({ umap_plot_obj() })

  output$download_umap <- downloadHandler(
    filename = function() paste0("umap_", input$umap_groupby, ".png"),
    content  = function(file)
      ggsave(file, umap_plot_obj(),
             width = input$umap_dl_w, height = input$umap_dl_h, dpi = input$umap_dpi)
  )

  markers_data <- eventReactive(input$markers, {
    req(input$cluster, input$umap_groupby)
    Idents(seurat_obj) <- input$umap_groupby
    FindMarkers(seurat_obj, ident.1 = input$cluster, only.pos = TRUE, max.cells.per.ident = 5000) |>
      rownames_to_column("gene") |>
      arrange(desc(avg_log2FC)) |>
      head(input$n_markers)
  })

  output$marker_table <- DT::renderDataTable({
    req(markers_data())
    datatable(markers_data(), filter = "top",
              options = list(pageLength = 10, scrollX = TRUE))
  })
  output$download_markers <- downloadHandler(
    filename = function() paste0("markers_", input$cluster, ".csv"),
    content  = function(file) write.csv(markers_data(), file, row.names = FALSE)
  )

  expr_plot_obj <- reactive({
    req(input$umap_groupby)
    Idents(seurat_obj) <- input$umap_groupby
    pt      <- input$expr_pt
    do_split <- isTRUE(input$use_split) && nchar(trimws(input$split_by)) > 0

    if (input$plot_type == "feature") {
      feats <- if (isTRUE(input$blend_mode))
        c(input$feature, input$feature2) else input$feature
      args <- list(seurat_obj, features = feats, reduction = "umap.harmony",
                   pt.size = 1, blend = isTRUE(input$blend_mode), combine=FALSE)
      if (do_split && !isTRUE(input$blend_mode)) args$split.by <- input$split_by
      do.call(FeaturePlot, args)

    } else if (input$plot_type == "vln") {
      args <- list(seurat_obj, features = input$feature,
                   group.by = input$umap_groupby, pt.size = 0)
      if (do_split) args$split.by <- input$split_by
      p <- do.call(VlnPlot, args)
      p + theme(axis.text.x = element_text(size = pt, angle = 45, hjust = 1))

    } else if (input$plot_type == "dot") {
      DotPlot(seurat_obj,
              features = unique(c(input$feature, input$feature2)),
              group.by = input$umap_groupby) +
        theme_minimal(base_size = pt) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))

    } else if (input$plot_type == "heatmap") {
      DoHeatmap(seurat_obj,
                features = unique(c(input$feature, input$feature2)),
                group.by = input$umap_groupby, angle = 45, size = pt / 3)
    }
  })

  output$expr_container <- renderUI({
    req(input$expr_width, input$expr_height)
    plotOutput("expr_plot", width  = paste0(input$expr_width,  "px"),
                            height = paste0(input$expr_height, "px"))
  })
  output$expr_plot <- renderPlot(expr_plot_obj())

  output$download_expr_plot <- downloadHandler(
    filename = function() paste0("expr_", input$feature, ".png"),
    content  = function(file)
      ggsave(file, expr_plot_obj(),
             width = input$expr_dl_w, height = input$expr_dl_h, dpi = input$expr_dpi)
  )

  # ================================================
  # TAB 2
  # ================================================

  output$deg_subset_value_selector <- renderUI({
    req(input$deg_subset_meta)
    choices <- sort(unique(as.character(seurat_obj@meta.data[[input$deg_subset_meta]])))
    selectInput("deg_subset_value", paste0("Subset to (", input$deg_subset_meta, "):"),
                choices = choices, selected = choices[1])
  })

  deg_available_groups <- reactive({
    req(input$deg_subset_meta, input$deg_subset_value, input$deg_compare_meta)
    Idents(seurat_obj) <- input$deg_subset_meta
    sub <- subset(seurat_obj, idents = input$deg_subset_value)
    sort(unique(as.character(sub@meta.data[[input$deg_compare_meta]])))
  })

  output$deg_group1_selector <- renderUI({
    req(deg_available_groups())
    groups <- deg_available_groups()
    selectInput("deg_group1", "Group 1 (ident.1):", choices = groups, selected = groups[1])
  })
  output$deg_group2_selector <- renderUI({
    req(deg_available_groups(), input$deg_group1)
    groups  <- deg_available_groups()
    default <- if (length(groups) >= 2) groups[groups != input$deg_group1][1] else groups[1]
    selectInput("deg_group2", "Group 2 (ident.2):", choices = groups, selected = default)
  })

  degs_data <- eventReactive(input$deg, {
    req(input$deg_subset_meta, input$deg_subset_value,
        input$deg_compare_meta, input$deg_group1, input$deg_group2)
    validate(need(input$deg_group1 != input$deg_group2, "Groups must be different."))
    Idents(seurat_obj) <- input$deg_subset_meta
    sub <- subset(seurat_obj, idents = input$deg_subset_value)
    Idents(sub) <- input$deg_compare_meta
    FindMarkers(sub, ident.1 = input$deg_group1, ident.2 = input$deg_group2,
                fc.slot = "counts") |>
      rownames_to_column("gene") |>
      mutate(direction = ifelse(avg_log2FC > 0, "UP", "DOWN"),
             log_p     = -log10(p_val_adj + 1e-300),
             abs_fc    = abs(avg_log2FC))
  })

  output$deg_table <- DT::renderDataTable({
    req(degs_data())
    df <- degs_data() |>
      filter(p_val_adj < input$pval_cut, abs_fc >= input$lfc_cut) |>
      arrange(p_val_adj) |> head(input$n_deg_table)
    datatable(df, filter = "top",
              options = list(pageLength = 15, scrollX = TRUE)) |>
      formatSignif(columns = c("p_val", "p_val_adj", "avg_log2FC"), digits = 3)
  })
  output$download_degs <- downloadHandler(
    filename = function() paste0("degs_", input$deg_subset_value, "_",
                                 input$deg_group1, "_vs_", input$deg_group2, ".csv"),
    content  = function(file) write.csv(degs_data(), file, row.names = FALSE)
  )

  volcano_plot_obj <- reactive({
    req(degs_data())
    df        <- degs_data() |> filter(p_val_adj < input$pval_cut, abs_fc >= input$lfc_cut)
    top_label <- df %>% arrange(desc(abs_fc)) %>% head(input$n_label)
    pt        <- input$volcano_pt
    ggplot(df, aes(x = avg_log2FC, y = log_p, color = direction,
                   text = paste0("Gene: ", gene,
                                 "<br>log2FC: ", round(avg_log2FC, 3),
                                 "<br>p.adj: ",  signif(p_val_adj, 3)))) +
      geom_point(alpha = 0.55, size = 1.5) +
      scale_color_manual(values = c("UP" = "#c0392b", "DOWN" = "#2980b9")) +
      geom_text_repel(data = top_label, aes(label = gene), size = pt / 3,
                      max.overlaps = 30, segment.color = "grey60", show.legend = FALSE) +
      geom_hline(yintercept = -log10(input$pval_cut), linetype = "dashed", color = "grey60", linewidth = 0.4) +
      geom_vline(xintercept = c(-input$lfc_cut, input$lfc_cut),
                 linetype = "dashed", color = "grey60", linewidth = 0.4) +
      labs(title = paste0(input$deg_group1, " vs ", input$deg_group2,
                          "  |  ", input$deg_subset_value),
           x = "log2 Fold Change", y = "-log10(adj. p-value)", color = NULL) +
      theme_minimal(base_size = pt) +
      theme(legend.position = "top",
            plot.title = element_text(face = "bold", size = pt + 1))
  })

  output$volcano_container <- renderUI({
    req(input$volcano_width, input$volcano_height)
    plotlyOutput("volcano_plot",
                 width  = paste0(input$volcano_width,  "px"),
                 height = paste0(input$volcano_height, "px"))
  })
  output$volcano_plot <- renderPlotly({ ggplotly(volcano_plot_obj(), tooltip = "text") })

  output$download_volcano <- downloadHandler(
    filename = function() paste0("volcano_", input$deg_subset_value, ".png"),
    content  = function(file)
      ggsave(file, volcano_plot_obj(),
             width = input$volcano_dl_w, height = input$volcano_dl_h, dpi = input$volcano_dpi)
  )

  enrich_results <- eventReactive(input$run_enrich, {
    req(degs_data())
    df         <- degs_data() |> filter(p_val_adj < input$pval_cut, abs_fc >= input$lfc_cut)
    up_genes   <- df %>% filter(direction == "UP")   %>% pull(gene) %>% unique()
    down_genes <- df %>% filter(direction == "DOWN") %>% pull(gene) %>% unique()
    db         <- input$enrich_db
    up_res <- if (length(up_genes) >= 5) {
      Sys.sleep(1); enrichR::enrichr(up_genes, db)[[db]]
    } else NULL
    down_res <- if (length(down_genes) >= 5) {
      Sys.sleep(1); enrichR::enrichr(down_genes, db)[[db]]
    } else NULL
    list(up = up_res, down = down_res)
  })

  make_enrich_plot <- function(df, color_high, title_txt, n_terms, pt) {
    df_top <- df %>% head(n_terms)
    ggplot(df_top, aes(x = reorder(Term, Combined.Score), y = Combined.Score,
                       fill = -log10(Adjusted.P.value + 1e-300))) +
      geom_col(color = "white", linewidth = 0.25) +
      scale_fill_gradient(low = "lightyellow", high = color_high, name = "-log10(adj.p)") +
      coord_flip() +
      theme_minimal(base_size = pt) +
      labs(title = title_txt, x = NULL, y = "Combined Score") +
      theme(plot.title = element_text(face = "bold"))
  }

  enrich_up_plot_obj   <- reactive({ req(enrich_results()$up);
    make_enrich_plot(enrich_results()$up,   "#c0392b", paste0("UP in ",   input$deg_group1), input$n_enrich, input$enrich_pt) })
  enrich_down_plot_obj <- reactive({ req(enrich_results()$down);
    make_enrich_plot(enrich_results()$down, "#2980b9", paste0("DOWN in ", input$deg_group1), input$n_enrich, input$enrich_pt) })

  output$enrich_up_plot   <- renderPlot({ req(enrich_up_plot_obj());   enrich_up_plot_obj()   },
    width = reactive(input$enrich_width), height = reactive(input$enrich_height))
  output$enrich_down_plot <- renderPlot({ req(enrich_down_plot_obj()); enrich_down_plot_obj() },
    width = reactive(input$enrich_width), height = reactive(input$enrich_height))

  enrich_tidy <- function(df, n) {
    if (is.null(df)) return(NULL)
    df %>% head(n) %>%
      select(Term, Overlap, P.value, Adjusted.P.value, Combined.Score, Genes) %>%
      mutate(across(where(is.numeric), ~ signif(.x, 3)))
  }
  output$enrich_up_table <- DT::renderDataTable({
    req(enrich_results()$up)
    datatable(enrich_tidy(enrich_results()$up, input$n_enrich),
              options = list(pageLength = 8, scrollX = TRUE), filter = "top")
  })
  output$enrich_down_table <- DT::renderDataTable({
    req(enrich_results()$down)
    datatable(enrich_tidy(enrich_results()$down, input$n_enrich),
              options = list(pageLength = 8, scrollX = TRUE), filter = "top")
  })
  output$download_enrich_up <- downloadHandler(
    filename = function() paste0("enrichment_UP_",   input$deg_subset_value, ".png"),
    content  = function(file) ggsave(file, enrich_up_plot_obj(),
      width = input$enrich_dl_w, height = input$enrich_dl_h, dpi = input$enrich_dpi))
  output$download_enrich_down <- downloadHandler(
    filename = function() paste0("enrichment_DOWN_", input$deg_subset_value, ".png"),
    content  = function(file) ggsave(file, enrich_down_plot_obj(),
      width = input$enrich_dl_w, height = input$enrich_dl_h, dpi = input$enrich_dpi))
  output$download_enrich_up_csv <- downloadHandler(
    filename = function() paste0("enrichment_UP_",   input$deg_subset_value, ".csv"),
    content  = function(file) write.csv(enrich_results()$up,   file, row.names = FALSE))
  output$download_enrich_down_csv <- downloadHandler(
    filename = function() paste0("enrichment_DOWN_", input$deg_subset_value, ".csv"),
    content  = function(file) write.csv(enrich_results()$down, file, row.names = FALSE))

  # ================================================
  # TAB 3
  # ================================================

  full_network <- reactive({ GetTFNetwork(seurat_obj) })

  network_data <- eventReactive(input$update_network, {
    net           <- full_network() %>% filter(Gain >= input$min_gain)
    all_tfs_in_net <- unique(net$tf)
    target_mode   <- nchar(trimws(input$tf_target_search)) > 0

    if (target_mode) {
      tgt   <- trimws(input$tf_target_search)
      edges <- net %>% filter(gene == tgt)
      validate(need(nrow(edges) > 0,
                    paste0("No TFs found regulating '", tgt, "' at current Gain threshold.")))
    } else {
      tf <- input$tf_center
      validate(need(tf %in% net$tf, paste0("TF '", tf, "' not found in network.")))
      primary <- net %>% filter(tf == !!tf) %>% arrange(desc(Gain)) %>% head(input$n_top_targets)
      edges   <- primary
      if (input$target_level == 2) {
        secondary <- net %>% filter(tf %in% primary$gene) %>%
          arrange(desc(Gain)) %>% head(input$n_top_targets * 2)
        edges <- bind_rows(primary, secondary) %>% distinct()
      }
    }

    if (isTRUE(input$tf_tf_only)) {
      edges <- edges %>% filter(gene %in% all_tfs_in_net)
      validate(need(nrow(edges) > 0, "No TF–TF interactions found."))
    }

    all_nodes  <- unique(c(edges$tf, edges$gene))
    is_tf_node <- all_nodes %in% all_tfs_in_net
    center_node <- if (target_mode) trimws(input$tf_target_search) else input$tf_center

    nodes_df <- data.frame(
      id    = all_nodes,
      label = all_nodes,
      color = ifelse(all_nodes == center_node, "#e74c3c",
              ifelse(is_tf_node, "#2980b9", "#7f8c8d")),
      shape = ifelse(is_tf_node, "diamond", "dot"),
      size  = ifelse(all_nodes == center_node, 24,
              ifelse(is_tf_node, 18, 12)),
      title = paste0("<b>", all_nodes, "</b><br>",
                     ifelse(all_nodes == center_node, "⭐ Center",
                     ifelse(is_tf_node, "🔵 TF", "⚫ Target"))),
      font.color = "#ffffff",
      stringsAsFactors = FALSE
    )

    edges_df <- data.frame(
      from   = edges$tf,
      to     = edges$gene,
      value  = scales::rescale(edges$Gain, c(1, 5)),
      color  = ifelse(edges$Cor > 0, "#2980b9", "#c0392b"),
      title  = paste0("Gain: ", round(edges$Gain, 3), " | Cor: ", round(edges$Cor, 3)),
      arrows = "to",
      stringsAsFactors = FALSE
    )

    list(nodes = nodes_df, edges = edges_df, raw_edges = edges)
  })

  output$network_plot <- renderVisNetwork({
    req(network_data())
    nd <- network_data()
    physics_solver <- switch(input$network_layout,
      "circle" = "repulsion", "fr" = "forceAtlas2Based",
      "kk" = "forceAtlas2Based", "drl" = "forceAtlas2Based", "forceAtlas2Based")

    visNetwork(nd$nodes, nd$edges, height = "660px", width = "100%",
               background = "#0f1b2d") %>%
      visEdges(smooth = list(type = "curvedCW", roundness = 0.15),
               arrows = list(to = list(enabled = TRUE, scaleFactor = 0.6))) %>%
      visNodes(font = list(size = 13, color = "#ffffff")) %>%
      visOptions(highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE),
                 nodesIdSelection = list(enabled = TRUE, style = "background:#162032;color:#1cb5bf;border:1px solid #1cb5bf;border-radius:4px;padding:3px")) %>%
      visInteraction(navigationButtons = TRUE, zoomView = TRUE,
                     tooltipDelay = 100) %>%
      visPhysics(solver = physics_solver,
                 forceAtlas2Based = list(gravitationalConstant = -60),
                 repulsion = list(nodeDistance = 120),
                 stabilization = list(iterations = 200)) %>%
      visLayout(randomSeed = 42,
                improvedLayout = (input$network_layout %in% c("nicely", "kk")))
  })

  output$download_network_full <- downloadHandler(
    filename = "tf_network_full.csv",
    content  = function(file) write.csv(full_network(), file, row.names = FALSE))
  output$download_network_filtered <- downloadHandler(
    filename = function() paste0("tf_network_filtered_", input$tf_center, ".csv"),
    content  = function(file) {
      req(network_data()); write.csv(network_data()$raw_edges, file, row.names = FALSE)
    })

  # ================================================
  # TAB 4
  # ================================================

  coexp_data <- eventReactive(input$run_coexp, {
    hub_df <- GetHubGenes(seurat_obj, n_hubs = input$coexp_n_hubs) %>%
      arrange(module, desc(kME))
    all_genes_df <- GetModules(seurat_obj) %>%
      filter(module != "grey") %>%
      rowwise() %>%
      mutate(kME_own = get(paste0("kME_", module))) %>%
      ungroup() %>%
      select(gene_name, module, kME_own) %>%
      arrange(module, desc(kME_own))
    module_sizes <- all_genes_df %>% count(module, name = "n_genes") %>% arrange(desc(n_genes))
    list(hub_df = hub_df, all_genes_df = all_genes_df,
         module_sizes = module_sizes,
         module_names = sort(unique(hub_df$module)))
  })

  output$module_size_table <- DT::renderDataTable({
    req(coexp_data())
    datatable(coexp_data()$module_sizes, filter = "top",
              options = list(pageLength = 10, scrollX = TRUE))
  })

  gene_lookup_result <- eventReactive(input$run_gene_lookup, {
    req(coexp_data(), nchar(trimws(input$gene_lookup)) > 0)
    g   <- trimws(input$gene_lookup)
    res <- coexp_data()$all_genes_df %>% filter(gene_name == g)
    if (nrow(res) == 0) paste0("Gene '", g, "' not found in any module.")
    else paste0("Gene:   ", g, "\nModule: ", res$module[1],
                "\nkME:    ", round(res$kME_own[1], 4))
  })
  output$gene_lookup_result <- renderPrint({ req(gene_lookup_result()); gene_lookup_result() })

  output$coexp_module_selector <- renderUI({
    req(coexp_data())
    selectInput("coexp_selected_module", "Module:", choices = coexp_data()$module_names,
                selected = coexp_data()$module_names[1])
  })

  output$hub_genes_table <- DT::renderDataTable({
    req(coexp_data())
    datatable(coexp_data()$hub_df, filter = "top",
              options = list(pageLength = 15, scrollX = TRUE))
  })
  output$download_hub_genes <- downloadHandler(
    filename = "hub_genes_by_module.csv",
    content  = function(file) { req(coexp_data()); write.csv(coexp_data()$hub_df, file, row.names = FALSE) })

  coexp_enrich_result <- eventReactive(input$run_coexp_enrich, {
    req(coexp_data(), input$coexp_selected_module)
    genes_in_module <- coexp_data()$all_genes_df %>%
      filter(module == input$coexp_selected_module) %>% pull(gene_name) %>% unique()
    validate(need(length(genes_in_module) >= 5,
                  paste0("Too few genes in module '", input$coexp_selected_module, "'")))
    db  <- input$coexp_enrich_db
    Sys.sleep(1)
    res <- enrichR::enrichr(genes_in_module, db)[[db]]
    validate(need(!is.null(res) && nrow(res) > 0, "No enrichment results."))
    res %>% arrange(Adjusted.P.value) %>% head(input$coexp_n_enrich)
  })

  output$enrich_module_title <- renderText({
    req(input$coexp_selected_module)
    paste0("Enrichment — ", input$coexp_selected_module, "  ·  ", input$coexp_enrich_db)
  })

  coexp_enrich_plot_obj <- reactive({
    req(coexp_enrich_result())
    df        <- coexp_enrich_result()
    mod       <- input$coexp_selected_module
    pt        <- input$coexp_enrich_pt
    mod_color <- tryCatch({
      GetModules(seurat_obj) %>% filter(module == mod) %>% pull(color) %>% unique() %>% .[1]
    }, error = function(e) "#1cb5bf")
    if (is.na(mod_color) || length(mod_color) == 0) mod_color <- "#1cb5bf"
    ggplot(df, aes(x = reorder(Term, Combined.Score), y = Combined.Score,
                   fill = -log10(Adjusted.P.value + 1e-300))) +
      geom_col(color = "white", linewidth = 0.25) +
      scale_fill_gradient(low = "lightyellow", high = mod_color, name = "-log10(adj.p)") +
      coord_flip() + theme_minimal(base_size = pt) +
      labs(title = paste0("Module: ", mod),
           subtitle = paste0(input$coexp_enrich_db, "  |  n = ",
                             nrow(coexp_data()$all_genes_df %>% filter(module == mod))),
           x = NULL, y = "Combined Score") +
      theme(plot.title = element_text(face = "bold"),
            plot.subtitle = element_text(color = "grey40"))
  })

  output$coexp_enrich_container <- renderUI({
    req(input$coexp_enrich_width, input$coexp_enrich_height)
    plotOutput("coexp_enrich_plot",
               width  = paste0(input$coexp_enrich_width,  "px"),
               height = paste0(input$coexp_enrich_height, "px"))
  })
  output$coexp_enrich_plot  <- renderPlot(coexp_enrich_plot_obj())
  output$coexp_enrich_table <- DT::renderDataTable({
    req(coexp_enrich_result())
    datatable(
      coexp_enrich_result() %>%
        select(Term, Overlap, P.value, Adjusted.P.value, Combined.Score, Genes) %>%
        mutate(across(where(is.numeric), ~ signif(.x, 3))),
      filter = "top", options = list(pageLength = 8, scrollX = TRUE))
  })
  output$download_coexp_enrich <- downloadHandler(
    filename = function() paste0("enrichment_module_", input$coexp_selected_module, ".csv"),
    content  = function(file) { req(coexp_enrich_result()); write.csv(coexp_enrich_result(), file, row.names = FALSE) })
  output$download_coexp_enrich_plot <- downloadHandler(
    filename = function() paste0("enrichment_barplot_", input$coexp_selected_module, ".png"),
    content  = function(file)
      ggsave(file, coexp_enrich_plot_obj(),
             width = input$coexp_enrich_dl_w, height = input$coexp_enrich_dl_h, dpi = input$coexp_enrich_dpi))

# ================================================
  # TAB 5
# ================================================

  output$reg_module_selector <- renderUI({
    m <- tryCatch({
      mods <- unique(as.character(GetModules(seurat_obj)$module))
      c("All modules", sort(mods[mods != "grey"]))
    }, error = function(e) {
      "All modules"
    })
    
    selectInput(
      "reg_selected_module",
      "Filter TFs by module:",
      choices  = m,
      selected = "All modules"
    )
  })
  
  # ----------------------------------------------------------------
  # Cached AverageExpression on a downsampled subset.
  # Computing AverageExpression on full atlas is the heaviest single
  # operation in the app. Downsampling to 5000 cells per cluster keeps
  # the means statistically stable while cutting memory by an order
  # of magnitude. Cached as a reactive so it runs once per session.
  # ----------------------------------------------------------------
  avg_expr_cached <- reactive({
    validate(need("Population_level3" %in% colnames(seurat_obj@meta.data),
                  "Metadata column 'Population_level3' not found in the Seurat object."))
    
    cells_keep <- unlist(lapply(
      unique(seurat_obj$Population_level3),
      function(cl) {
        cells <- colnames(seurat_obj)[seurat_obj$Population_level3 == cl]
        if (length(cells) > 5000) sample(cells, 5000) else cells
      }
    ))
    
    am <- AverageExpression(
      subset(seurat_obj, cells = cells_keep),
      group.by = "Population_level3"
    )[[1]]
    
    as.matrix(am)
  })
  
  reg_heatmap_data <- eventReactive(input$run_reg_heatmap, {
    withProgress(message = "Building TF regulon heatmap...", value = 0, {
      
      incProgress(0.05, detail = "Loading TF network...")
      net <- GetTFNetwork(seurat_obj) %>% filter(Gain >= input$reg_min_gain)
      validate(need(nrow(net) > 0, "No edges after Gain filter."))
      
      incProgress(0.15, detail = "Average expression per cluster...")
      avg_mat      <- avg_expr_cached()
      clusters     <- colnames(avg_mat)
      genes_in_mat <- rownames(avg_mat)
      
      incProgress(0.30, detail = "Filtering TFs...")
      module_tfs <- character(0)
      if (!is.null(input$reg_selected_module) && input$reg_selected_module != "All modules") {
        tryCatch({
          module_tfs <- GetModules(seurat_obj) %>%
            filter(module == input$reg_selected_module) %>% pull(gene_name) %>% unique()
        }, error = function(e) {})
      }
      manual_tfs <- if (nchar(trimws(input$reg_tf_filter)) > 0)
        strsplit(input$reg_tf_filter, ",\\s*")[[1]] %>% trimws() else character(0)
      
      all_tfs <- unique(net$tf)
      if (length(module_tfs) > 0) {
        all_tfs <- intersect(all_tfs, module_tfs)
        validate(need(length(all_tfs) > 0,
                      paste0("No TFs from module '", input$reg_selected_module, "'.")))
      }
      tfs_valid <- all_tfs[all_tfs %in% genes_in_mat]
      min_sz    <- input$reg_min_regulon_size
      tfs_valid <- tfs_valid[vapply(tfs_valid, function(tf_name) {
        pos_n <- sum(net$tf == tf_name & net$Cor > 0 & net$gene %in% genes_in_mat)
        neg_n <- sum(net$tf == tf_name & net$Cor < 0 & net$gene %in% genes_in_mat)
        max(pos_n, neg_n) >= min_sz
      }, logical(1))]
      validate(need(length(tfs_valid) > 0, "No TFs pass filters."))
      
      if (length(manual_tfs) > 0) {
        tfs_use <- intersect(manual_tfs, tfs_valid)
        validate(need(length(tfs_use) > 0, "Specified TFs not in filtered network."))
      } else {
        reg_sizes <- vapply(tfs_valid, function(tf_name) sum(net$tf == tf_name), integer(1))
        tfs_use   <- tfs_valid[order(reg_sizes, decreasing = TRUE)][
          seq_len(min(input$reg_top_n_tfs, length(tfs_valid)))]
      }
      
      incProgress(0.50, detail = "Building expression matrix...")
      
      safe_regulon_mean <- function(genes) {
        if (length(genes) == 0) return(rep(NA_real_, length(clusters)))
        if (length(genes) == 1) return(as.numeric(avg_mat[genes, ]))
        as.numeric(colMeans(avg_mat[genes, , drop = FALSE], na.rm = TRUE))
      }
      
      rows <- lapply(tfs_use, function(tf_name) {
        pos_genes <- net %>% filter(tf == !!tf_name, Cor > 0) %>%
          pull(gene) %>% unique() %>% intersect(genes_in_mat)
        neg_genes <- net %>% filter(tf == !!tf_name, Cor < 0) %>%
          pull(gene) %>% unique() %>% intersect(genes_in_mat)
        data.frame(
          tf       = tf_name,
          cluster  = clusters,
          tf_expr  = as.numeric(avg_mat[tf_name, ]),
          pos_reg  = safe_regulon_mean(pos_genes),
          neg_reg  = safe_regulon_mean(neg_genes),
          n_pos    = length(pos_genes),
          n_neg    = length(neg_genes),
          stringsAsFactors = FALSE
        )
      })
      
      df_wide <- bind_rows(rows) %>% group_by(tf) %>%
        mutate(tf_expr_scaled = { rng <- range(tf_expr, na.rm = TRUE);
        if (diff(rng) == 0) rep(0.5, n()) else (tf_expr - rng[1]) / diff(rng) }) %>%
        ungroup()
      
      df_long <- df_wide %>%
        pivot_longer(cols = c(pos_reg, neg_reg),
                     names_to = "regulon_dir", values_to = "mean_expr") %>%
        mutate(regulon_dir = dplyr::case_when(
          regulon_dir == "pos_reg" ~ "Positive regulon  (Cor > 0)",
          regulon_dir == "neg_reg" ~ "Negative regulon  (Cor < 0)"),
          regulon_dir = factor(regulon_dir,
                               levels = c("Positive regulon  (Cor > 0)", "Negative regulon  (Cor < 0)")))
      
      incProgress(1, detail = "Done.")
      list(df_long = df_long, df_wide = df_wide, tfs_use = tfs_use,
           clusters = clusters, module = input$reg_selected_module, net = net)
    })  # close withProgress
  })    # close eventReactive

  output$reg_enrich_tf_selector <- renderUI({
    if (isTruthy(reg_heatmap_data())) {
      tfs <- reg_heatmap_data()$tfs_use
      if (length(tfs) > 0)
        selectInput("reg_enrich_tf", "TF for enrichment:", choices = tfs, selected = tfs[1])
      else helpText("Build the heatmap first.")
    } else helpText("Build the heatmap first.")
  })

  output$reg_status <- renderText({
    if (!isTruthy(reg_heatmap_data())) return("")
    d  <- reg_heatmap_data()
    mt <- if (!is.null(d$module) && d$module != "All modules") paste0("  |  Module: ", d$module) else ""
    paste0("Loaded: ", length(d$tfs_use), " TFs  ×  ", length(d$clusters), " clusters", mt)
  })

  reg_heatmap_plot_obj <- reactive({
    req(reg_heatmap_data())
    df   <- reg_heatmap_data()$df_long
    dmin <- input$reg_dot_range[1]; dmax <- input$reg_dot_range[2]; pt <- input$reg_pt

    tf_order <- df %>%
      filter(regulon_dir == "Positive regulon  (Cor > 0)") %>%
      group_by(tf) %>% summarise(total_pos = sum(mean_expr, na.rm = TRUE), .groups = "drop") %>%
      arrange(desc(total_pos)) %>% pull(tf)
    df <- df %>% mutate(tf = factor(tf, levels = rev(tf_order)))

    dot_theme <- theme_minimal(base_size = pt) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = pt - 1),
            axis.text.y = element_text(size = pt - 1),
            panel.grid.major = element_line(color = "grey92"),
            legend.position = "right",
            plot.title = element_text(face = "bold", size = pt + 1),
            legend.key.size = unit(0.45, "cm"))

    make_dot <- function(data, fill_high, fill_name, title_txt) {
      ggplot(data, aes(x = cluster, y = tf, size = tf_expr_scaled, fill = mean_expr)) +
        geom_point(shape = 21, color = "grey30", stroke = 0.3, alpha = 0.9) +
        scale_size_continuous(range = c(dmin, dmax), name = "TF expr\n(scaled 0–1)") +
        scale_fill_gradient(low = "white", high = fill_high, na.value = "grey85", name = fill_name) +
        labs(title = title_txt, x = NULL, y = "Transcription Factor") + dot_theme
    }

    p_pos <- make_dot(df %>% filter(regulon_dir == "Positive regulon  (Cor > 0)"),
                      "firebrick", "Mean expr\npos. regulon", "Positive regulon  (Cor > 0)")
    p_neg <- make_dot(df %>% filter(regulon_dir == "Negative regulon  (Cor < 0)"),
                      "steelblue4", "Mean expr\nneg. regulon", "Negative regulon  (Cor < 0)") +
      labs(x = "Cluster")

    title_suffix <- if (!is.null(reg_heatmap_data()$module) &&
                        reg_heatmap_data()$module != "All modules")
      paste0("  |  Module: ", reg_heatmap_data()$module) else ""

    p_pos / p_neg +
      plot_annotation(
        title = paste0("TF Regulon Heatmap", title_suffix,
                       "  |  dot size = TF expr (scaled)  |  dot color = mean regulon expr"),
        theme = theme(plot.title = element_text(size = pt, color = "grey40", hjust = 0.5)))
  })

  output$reg_heatmap_container <- renderUI({
    req(input$reg_width, input$reg_height)
    plotOutput("reg_heatmap_plot",
               width  = paste0(input$reg_width,  "px"),
               height = paste0(input$reg_height, "px"))
  })
  output$reg_heatmap_plot  <- renderPlot({ req(reg_heatmap_plot_obj()); reg_heatmap_plot_obj() })

  output$reg_summary_table <- DT::renderDataTable({
    req(reg_heatmap_data())
    datatable(
      reg_heatmap_data()$df_wide %>% select(tf, n_pos, n_neg) %>% distinct() %>%
        mutate(total_regulon = n_pos + n_neg) %>% arrange(desc(total_regulon)),
      filter = "top", options = list(pageLength = 15, scrollX = TRUE))
  })

  output$download_reg_plot <- downloadHandler(
    filename = function() {
      mod <- if (!is.null(input$reg_selected_module) && input$reg_selected_module != "All modules")
        paste0("_", input$reg_selected_module) else ""
      paste0("tf_regulon_heatmap", mod, ".png")
    },
    content = function(file) {
      req(reg_heatmap_plot_obj())
      n_tfs  <- length(reg_heatmap_data()$tfs_use)
      height <- max(input$reg_dl_h, n_tfs * 0.35 * 2)
      ggsave(file, reg_heatmap_plot_obj(), width = input$reg_dl_w,
             height = height, dpi = input$reg_dpi, limitsize = FALSE)
    })
  output$download_reg_csv <- downloadHandler(
    filename = function() {
      mod <- if (!is.null(input$reg_selected_module) && input$reg_selected_module != "All modules")
        paste0("_", input$reg_selected_module) else ""
      paste0("tf_regulon_heatmap_data", mod, ".csv")
    },
    content = function(file) { req(reg_heatmap_data()); write.csv(reg_heatmap_data()$df_wide, file, row.names = FALSE) })

  # EnrichR on regulon
  get_reg_genes <- function(sign = c("pos", "neg")) {
    sign <- match.arg(sign)
    net  <- if (isTruthy(reg_heatmap_data())) reg_heatmap_data()$net
            else GetTFNetwork(seurat_obj) %>% filter(Gain >= input$reg_min_gain)
    tf_sel <- input$reg_enrich_tf
    validate(need(!is.null(tf_sel) && tf_sel != "", "Build the heatmap first."))
    validate(need(tf_sel %in% net$tf, paste0("TF '", tf_sel, "' not in network.")))
    genes <- if (sign == "pos")
      net %>% filter(tf == tf_sel, Cor > 0) %>% pull(gene) %>% unique()
    else
      net %>% filter(tf == tf_sel, Cor < 0) %>% pull(gene) %>% unique()
    validate(need(length(genes) >= 5, paste0("Too few ", sign, " targets (n < 5).")))
    genes
  }

  run_enrichr_tidy <- function(genes, db) {
    Sys.sleep(1)
    res <- enrichR::enrichr(genes, databases = db)
    df  <- res[[db]]
    if (is.null(df) || nrow(df) == 0) return(NULL)
    df %>% arrange(Adjusted.P.value) %>%
      mutate(neg_log10_p = -log10(Adjusted.P.value + 1e-300),
             Term = str_trunc(Term, 65))
  }

  reg_enrich_results <- eventReactive(input$run_reg_enrich, {
    db <- input$reg_enrich_db
    list(pos = run_enrichr_tidy(get_reg_genes("pos"), db),
         neg = run_enrichr_tidy(get_reg_genes("neg"), db),
         tf  = input$reg_enrich_tf, db = db)
  })

  reg_enrich_pos_df <- reactive({ req(reg_enrich_results()$pos); reg_enrich_results()$pos %>% head(input$reg_enrich_n) })
  reg_enrich_neg_df <- reactive({ req(reg_enrich_results()$neg); reg_enrich_results()$neg %>% head(input$reg_enrich_n) })

  make_reg_enrich_plot <- function(df, low_col, high_col, tf_name, db, direction) {
    ggplot(df, aes(x = Combined.Score, y = reorder(Term, Combined.Score), fill = neg_log10_p)) +
      geom_col(color = "white", linewidth = 0.25) +
      scale_fill_gradient(low = low_col, high = high_col, name = "-log10\n(adj.p)") +
      labs(title = paste0(direction, " regulon — ", tf_name),
           subtitle = paste0(db, "  |  top ", nrow(df), " terms"),
           x = "Combined Score", y = NULL) +
      theme_minimal(base_size = input$reg_enrich_pt) +
      theme(plot.title = element_text(face = "bold"),
            plot.subtitle = element_text(color = "grey40"))
  }

  reg_enrich_pos_plot_obj <- reactive({
    req(reg_enrich_pos_df())
    make_reg_enrich_plot(reg_enrich_pos_df(), "#f4a582", "#ca0020",
                         reg_enrich_results()$tf, reg_enrich_results()$db, "Positive")
  })
  reg_enrich_neg_plot_obj <- reactive({
    req(reg_enrich_neg_df())
    make_reg_enrich_plot(reg_enrich_neg_df(), "#92c5de", "#0571b0",
                         reg_enrich_results()$tf, reg_enrich_results()$db, "Negative")
  })

  output$reg_enrich_pos_table <- DT::renderDataTable({
    req(reg_enrich_pos_df())
    datatable(reg_enrich_pos_df() %>%
        select(Term, Overlap, Adjusted.P.value, Combined.Score, Genes) %>%
        mutate(across(where(is.numeric), ~ signif(.x, 3))),
      filter = "top", options = list(pageLength = 8, scrollX = TRUE))
  })
  output$reg_enrich_neg_table <- DT::renderDataTable({
    req(reg_enrich_neg_df())
    datatable(reg_enrich_neg_df() %>%
        select(Term, Overlap, Adjusted.P.value, Combined.Score, Genes) %>%
        mutate(across(where(is.numeric), ~ signif(.x, 3))),
      filter = "top", options = list(pageLength = 8, scrollX = TRUE))
  })

  output$reg_enrich_pos_container <- renderUI({
    req(input$reg_enrich_width, input$reg_enrich_height)
    plotOutput("reg_enrich_pos_plot",
               width  = paste0(input$reg_enrich_width,  "px"),
               height = paste0(input$reg_enrich_height, "px"))
  })
  output$reg_enrich_neg_container <- renderUI({
    req(input$reg_enrich_width, input$reg_enrich_height)
    plotOutput("reg_enrich_neg_plot",
               width  = paste0(input$reg_enrich_width,  "px"),
               height = paste0(input$reg_enrich_height, "px"))
  })
  output$reg_enrich_pos_plot <- renderPlot(reg_enrich_pos_plot_obj())
  output$reg_enrich_neg_plot <- renderPlot(reg_enrich_neg_plot_obj())

  output$download_reg_enrich_pos_png <- downloadHandler(
    filename = function() paste0("POS_", input$reg_enrich_tf, "_", input$reg_enrich_db, ".png"),
    content  = function(file) ggsave(file, reg_enrich_pos_plot_obj(),
      width = input$reg_enrich_dl_w, height = input$reg_enrich_dl_h, dpi = input$reg_enrich_dpi))
  output$download_reg_enrich_neg_png <- downloadHandler(
    filename = function() paste0("NEG_", input$reg_enrich_tf, "_", input$reg_enrich_db, ".png"),
    content  = function(file) ggsave(file, reg_enrich_neg_plot_obj(),
      width = input$reg_enrich_dl_w, height = input$reg_enrich_dl_h, dpi = input$reg_enrich_dpi))
  output$download_reg_enrich_pos_csv <- downloadHandler(
    filename = function() paste0("POS_", input$reg_enrich_tf, "_", input$reg_enrich_db, ".csv"),
    content  = function(file) write.csv(reg_enrich_results()$pos, file, row.names = FALSE))
  output$download_reg_enrich_neg_csv <- downloadHandler(
    filename = function() paste0("NEG_", input$reg_enrich_tf, "_", input$reg_enrich_db, ".csv"),
    content  = function(file) write.csv(reg_enrich_results()$neg, file, row.names = FALSE))


  }# end server function
}

