# ============================================================
# server.R — Server builder for SeuratAtlasExplorer
# ============================================================

#' Build the Shiny server for the Atlas Explorer
#'
#' Returns a standard Shiny \code{server} function that closes over the
#' supplied \code{seurat_obj} and \code{metadata_choices}, so no global
#' variables are required.
#'
#' @param seurat_obj      A Seurat object.
#' @param metadata_choices Character vector from \code{get_metadata_choices()}.
#' @return A Shiny server function.
#' @keywords internal
atlas_server <- function(seurat_obj, metadata_choices) {

  function(input, output, session) {

    # ==========================================================
    # TAB 1 — UMAP + Expression
    # ==========================================================

    output$cluster_selector <- shiny::renderUI({
      shiny::req(input$umap_groupby)
      choices <- sort(unique(as.character(seurat_obj@meta.data[[input$umap_groupby]])))
      shiny::selectInput("cluster", paste0("Cluster (", input$umap_groupby, "):"),
                         choices = choices, selected = choices[1])
    })

    umap_plot_obj <- shiny::reactive({
      shiny::req(input$umap_groupby)
      Seurat::DimPlot(seurat_obj, reduction = "umap.harmony",
                      group.by = input$umap_groupby, label = TRUE,
                      label.size = input$umap_pt / 3) +
        ggplot2::theme_minimal(base_size = input$umap_pt)
    })

    output$umap_container <- shiny::renderUI({
      shiny::req(input$umap_width, input$umap_height)
      shiny::plotOutput("umap", width  = paste0(input$umap_width,  "px"),
                                height = paste0(input$umap_height, "px"))
    })
    output$umap <- shiny::renderPlot({ umap_plot_obj() })

    output$download_umap <- shiny::downloadHandler(
      filename = function() paste0("umap_", input$umap_groupby, ".png"),
      content  = function(file)
        ggplot2::ggsave(file, umap_plot_obj(),
                        width = input$umap_dl_w, height = input$umap_dl_h, dpi = input$umap_dpi)
    )

    markers_data <- shiny::eventReactive(input$markers, {
      shiny::req(input$cluster, input$umap_groupby)
      Seurat::Idents(seurat_obj) <- input$umap_groupby
      Seurat::FindMarkers(seurat_obj, ident.1 = input$cluster, only.pos = TRUE) |>
        tibble::rownames_to_column("gene") |>
        dplyr::arrange(dplyr::desc(.data$avg_log2FC)) |>
        utils::head(input$n_markers)
    })

    output$marker_table <- DT::renderDT({
      shiny::req(markers_data())
      DT::datatable(markers_data(), filter = "top",
                    options = list(pageLength = 10, scrollX = TRUE))
    })
    output$download_markers <- shiny::downloadHandler(
      filename = function() paste0("markers_", input$cluster, ".csv"),
      content  = function(file) utils::write.csv(markers_data(), file, row.names = FALSE)
    )

    expr_plot_obj <- shiny::reactive({
      shiny::req(input$umap_groupby)
      Seurat::Idents(seurat_obj) <- input$umap_groupby
      pt       <- input$expr_pt
      do_split <- isTRUE(input$use_split) && nchar(trimws(input$split_by)) > 0

      if (input$plot_type == "feature") {
        feats <- if (isTRUE(input$blend_mode)) c(input$feature, input$feature2) else input$feature
        args  <- list(seurat_obj, features = feats, reduction = "umap.harmony",
                      pt.size = 1, blend = isTRUE(input$blend_mode))
        if (do_split && !isTRUE(input$blend_mode)) args$split.by <- input$split_by
        do.call(Seurat::FeaturePlot, args)

      } else if (input$plot_type == "vln") {
        args <- list(seurat_obj, features = input$feature,
                     group.by = input$umap_groupby, pt.size = 0)
        if (do_split) args$split.by <- input$split_by
        p <- do.call(Seurat::VlnPlot, args)
        p & ggplot2::theme(axis.text.x = ggplot2::element_text(size = pt, angle = 45, hjust = 1))

      } else if (input$plot_type == "dot") {
        Seurat::DotPlot(seurat_obj,
                        features = unique(c(input$feature, input$feature2)),
                        group.by = input$umap_groupby) +
          ggplot2::theme_minimal(base_size = pt) +
          ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

      } else if (input$plot_type == "heatmap") {
        Seurat::DoHeatmap(seurat_obj,
                          features = unique(c(input$feature, input$feature2)),
                          group.by = input$umap_groupby, angle = 45, size = pt / 3)
      }
    })

    output$expr_container <- shiny::renderUI({
      shiny::req(input$expr_width, input$expr_height)
      shiny::plotOutput("expr_plot", width  = paste0(input$expr_width,  "px"),
                                     height = paste0(input$expr_height, "px"))
    })
    output$expr_plot <- shiny::renderPlot(expr_plot_obj())
    output$download_expr_plot <- shiny::downloadHandler(
      filename = function() paste0("expr_", input$feature, ".png"),
      content  = function(file)
        ggplot2::ggsave(file, expr_plot_obj(),
                        width = input$expr_dl_w, height = input$expr_dl_h, dpi = input$expr_dpi)
    )

    # ==========================================================
    # TAB 2 — DEGs + Volcano + Enrichment
    # ==========================================================

    output$deg_subset_value_selector <- shiny::renderUI({
      shiny::req(input$deg_subset_meta)
      choices <- sort(unique(as.character(seurat_obj@meta.data[[input$deg_subset_meta]])))
      shiny::selectInput("deg_subset_value", paste0("Subset to (", input$deg_subset_meta, "):"),
                         choices = choices, selected = choices[1])
    })

    deg_available_groups <- shiny::reactive({
      shiny::req(input$deg_subset_meta, input$deg_subset_value, input$deg_compare_meta)
      Seurat::Idents(seurat_obj) <- input$deg_subset_meta
      sub <- subset(seurat_obj, idents = input$deg_subset_value)
      sort(unique(as.character(sub@meta.data[[input$deg_compare_meta]])))
    })

    output$deg_group1_selector <- shiny::renderUI({
      shiny::req(deg_available_groups())
      groups <- deg_available_groups()
      shiny::selectInput("deg_group1", "Group 1 (ident.1):", choices = groups, selected = groups[1])
    })
    output$deg_group2_selector <- shiny::renderUI({
      shiny::req(deg_available_groups(), input$deg_group1)
      groups  <- deg_available_groups()
      default <- if (length(groups) >= 2) groups[groups != input$deg_group1][1] else groups[1]
      shiny::selectInput("deg_group2", "Group 2 (ident.2):", choices = groups, selected = default)
    })

    degs_data <- shiny::eventReactive(input$deg, {
      shiny::req(input$deg_subset_meta, input$deg_subset_value,
                 input$deg_compare_meta, input$deg_group1, input$deg_group2)
      shiny::validate(shiny::need(input$deg_group1 != input$deg_group2, "Groups must be different."))
      Seurat::Idents(seurat_obj) <- input$deg_subset_meta
      sub <- subset(seurat_obj, idents = input$deg_subset_value)
      Seurat::Idents(sub) <- input$deg_compare_meta
      Seurat::FindMarkers(sub, ident.1 = input$deg_group1, ident.2 = input$deg_group2,
                          fc.slot = "counts") |>
        tibble::rownames_to_column("gene") |>
        dplyr::mutate(
          direction = ifelse(.data$avg_log2FC > 0, "UP", "DOWN"),
          log_p     = -log10(.data$p_val_adj + 1e-300),
          abs_fc    = abs(.data$avg_log2FC)
        )
    })

    output$deg_table <- DT::renderDT({
      shiny::req(degs_data())
      df <- degs_data() |>
        dplyr::filter(.data$p_val_adj < input$pval_cut, .data$abs_fc >= input$lfc_cut) |>
        dplyr::arrange(.data$p_val_adj) |>
        utils::head(input$n_deg_table)
      DT::datatable(df, filter = "top", options = list(pageLength = 15, scrollX = TRUE)) |>
        DT::formatSignif(columns = c("p_val", "p_val_adj", "avg_log2FC"), digits = 3)
    })
    output$download_degs <- shiny::downloadHandler(
      filename = function() paste0("degs_", input$deg_subset_value, "_",
                                   input$deg_group1, "_vs_", input$deg_group2, ".csv"),
      content  = function(file) utils::write.csv(degs_data(), file, row.names = FALSE)
    )

    volcano_plot_obj <- shiny::reactive({
      shiny::req(degs_data())
      df        <- degs_data() |>
        dplyr::filter(.data$p_val_adj < input$pval_cut, .data$abs_fc >= input$lfc_cut)
      top_label <- df %>% dplyr::arrange(dplyr::desc(.data$abs_fc)) %>% utils::head(input$n_label)
      pt        <- input$volcano_pt
      ggplot2::ggplot(df, ggplot2::aes(
          x = .data$avg_log2FC, y = .data$log_p, color = .data$direction,
          text = paste0("Gene: ", .data$gene,
                        "<br>log2FC: ", round(.data$avg_log2FC, 3),
                        "<br>p.adj: ",  signif(.data$p_val_adj, 3)))) +
        ggplot2::geom_point(alpha = 0.55, size = 1.5) +
        ggplot2::scale_color_manual(values = c("UP" = "#c0392b", "DOWN" = "#2980b9")) +
        ggrepel::geom_text_repel(data = top_label,
                                 ggplot2::aes(label = .data$gene),
                                 size = pt / 3, max.overlaps = 30,
                                 segment.color = "grey60", show.legend = FALSE) +
        ggplot2::geom_hline(yintercept = -log10(input$pval_cut),
                            linetype = "dashed", color = "grey60", linewidth = 0.4) +
        ggplot2::geom_vline(xintercept = c(-input$lfc_cut, input$lfc_cut),
                            linetype = "dashed", color = "grey60", linewidth = 0.4) +
        ggplot2::labs(title = paste0(input$deg_group1, " vs ", input$deg_group2,
                                     "  |  ", input$deg_subset_value),
                      x = "log2 Fold Change", y = "-log10(adj. p-value)", color = NULL) +
        ggplot2::theme_minimal(base_size = pt) +
        ggplot2::theme(legend.position = "top",
                       plot.title = ggplot2::element_text(face = "bold", size = pt + 1))
    })

    output$volcano_container <- shiny::renderUI({
      shiny::req(input$volcano_width, input$volcano_height)
      plotly::plotlyOutput("volcano_plot",
                           width  = paste0(input$volcano_width,  "px"),
                           height = paste0(input$volcano_height, "px"))
    })
    output$volcano_plot <- plotly::renderPlotly({ plotly::ggplotly(volcano_plot_obj(), tooltip = "text") })
    output$download_volcano <- shiny::downloadHandler(
      filename = function() paste0("volcano_", input$deg_subset_value, ".png"),
      content  = function(file)
        ggplot2::ggsave(file, volcano_plot_obj(),
                        width = input$volcano_dl_w, height = input$volcano_dl_h, dpi = input$volcano_dpi)
    )

    enrich_results <- shiny::eventReactive(input$run_enrich, {
      shiny::req(degs_data())
      df         <- degs_data() |>
        dplyr::filter(.data$p_val_adj < input$pval_cut, .data$abs_fc >= input$lfc_cut)
      up_genes   <- df %>% dplyr::filter(.data$direction == "UP")   %>% dplyr::pull(.data$gene) %>% unique()
      down_genes <- df %>% dplyr::filter(.data$direction == "DOWN") %>% dplyr::pull(.data$gene) %>% unique()
      db         <- input$enrich_db
      up_res   <- if (length(up_genes)   >= 5) run_enrichr_tidy(up_genes,   db) else NULL
      down_res <- if (length(down_genes) >= 5) run_enrichr_tidy(down_genes, db) else NULL
      list(up = up_res, down = down_res)
    })

    enrich_up_plot_obj   <- shiny::reactive({
      shiny::req(enrich_results()$up)
      make_enrich_plot(enrich_results()$up,   "#c0392b",
                       paste0("UP in ",   input$deg_group1), input$n_enrich, input$enrich_pt)
    })
    enrich_down_plot_obj <- shiny::reactive({
      shiny::req(enrich_results()$down)
      make_enrich_plot(enrich_results()$down, "#2980b9",
                       paste0("DOWN in ", input$deg_group1), input$n_enrich, input$enrich_pt)
    })

    output$enrich_up_plot   <- shiny::renderPlot({ shiny::req(enrich_up_plot_obj());   enrich_up_plot_obj()   },
                                                  width = shiny::reactive(input$enrich_width),
                                                  height = shiny::reactive(input$enrich_height))
    output$enrich_down_plot <- shiny::renderPlot({ shiny::req(enrich_down_plot_obj()); enrich_down_plot_obj() },
                                                  width = shiny::reactive(input$enrich_width),
                                                  height = shiny::reactive(input$enrich_height))

    enrich_tidy <- function(df, n) {
      if (is.null(df)) return(NULL)
      utils::head(df, n) %>%
        dplyr::select("Term","Overlap","P.value","Adjusted.P.value","Combined.Score","Genes") %>%
        dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~ signif(.x, 3)))
    }
    output$enrich_up_table <- DT::renderDT({
      shiny::req(enrich_results()$up)
      DT::datatable(enrich_tidy(enrich_results()$up, input$n_enrich),
                    options = list(pageLength = 8, scrollX = TRUE), filter = "top")
    })
    output$enrich_down_table <- DT::renderDT({
      shiny::req(enrich_results()$down)
      DT::datatable(enrich_tidy(enrich_results()$down, input$n_enrich),
                    options = list(pageLength = 8, scrollX = TRUE), filter = "top")
    })
    output$download_enrich_up <- shiny::downloadHandler(
      filename = function() paste0("enrichment_UP_",   input$deg_subset_value, ".png"),
      content  = function(file) ggplot2::ggsave(file, enrich_up_plot_obj(),
        width = input$enrich_dl_w, height = input$enrich_dl_h, dpi = input$enrich_dpi))
    output$download_enrich_down <- shiny::downloadHandler(
      filename = function() paste0("enrichment_DOWN_", input$deg_subset_value, ".png"),
      content  = function(file) ggplot2::ggsave(file, enrich_down_plot_obj(),
        width = input$enrich_dl_w, height = input$enrich_dl_h, dpi = input$enrich_dpi))
    output$download_enrich_up_csv <- shiny::downloadHandler(
      filename = function() paste0("enrichment_UP_",   input$deg_subset_value, ".csv"),
      content  = function(file) utils::write.csv(enrich_results()$up,   file, row.names = FALSE))
    output$download_enrich_down_csv <- shiny::downloadHandler(
      filename = function() paste0("enrichment_DOWN_", input$deg_subset_value, ".csv"),
      content  = function(file) utils::write.csv(enrich_results()$down, file, row.names = FALSE))

    # ==========================================================
    # TAB 3 — TF Network
    # ==========================================================

    full_network <- shiny::reactive({ GetTFNetwork(seurat_obj) })

    network_data <- shiny::eventReactive(input$update_network, {
      net            <- full_network() %>% dplyr::filter(.data$Gain >= input$min_gain)
      all_tfs_in_net <- unique(net$tf)
      target_mode    <- nchar(trimws(input$tf_target_search)) > 0

      if (target_mode) {
        tgt   <- trimws(input$tf_target_search)
        edges <- net %>% dplyr::filter(.data$gene == tgt)
        shiny::validate(shiny::need(nrow(edges) > 0,
                                    paste0("No TFs regulating '", tgt, "' at current Gain.")))
      } else {
        tf <- input$tf_center
        shiny::validate(shiny::need(tf %in% net$tf, paste0("TF '", tf, "' not found.")))
        primary <- net %>% dplyr::filter(.data$tf == tf) %>%
          dplyr::arrange(dplyr::desc(.data$Gain)) %>% utils::head(input$n_top_targets)
        edges <- primary
        if (input$target_level == 2) {
          secondary <- net %>% dplyr::filter(.data$tf %in% primary$gene) %>%
            dplyr::arrange(dplyr::desc(.data$Gain)) %>% utils::head(input$n_top_targets * 2)
          edges <- dplyr::bind_rows(primary, secondary) %>% dplyr::distinct()
        }
      }

      if (isTRUE(input$tf_tf_only)) {
        edges <- edges %>% dplyr::filter(.data$gene %in% all_tfs_in_net)
        shiny::validate(shiny::need(nrow(edges) > 0, "No TF\u2013TF interactions found."))
      }

      all_nodes   <- unique(c(edges$tf, edges$gene))
      is_tf_node  <- all_nodes %in% all_tfs_in_net
      center_node <- if (target_mode) trimws(input$tf_target_search) else input$tf_center

      nodes_df <- data.frame(
        id    = all_nodes, label = all_nodes,
        color = ifelse(all_nodes == center_node, "#e74c3c",
                ifelse(is_tf_node, "#2980b9", "#7f8c8d")),
        shape = ifelse(is_tf_node, "diamond", "dot"),
        size  = ifelse(all_nodes == center_node, 24, ifelse(is_tf_node, 18, 12)),
        title = paste0("<b>", all_nodes, "</b><br>",
                       ifelse(all_nodes == center_node, "\u2B50 Center",
                       ifelse(is_tf_node, "\U0001F535 TF", "\u26AB Target"))),
        font.color = "#ffffff", stringsAsFactors = FALSE
      )
      edges_df <- data.frame(
        from = edges$tf, to = edges$gene,
        value  = scales::rescale(edges$Gain, c(1, 5)),
        color  = ifelse(edges$Cor > 0, "#2980b9", "#c0392b"),
        title  = paste0("Gain: ", round(edges$Gain, 3), " | Cor: ", round(edges$Cor, 3)),
        arrows = "to", stringsAsFactors = FALSE
      )
      list(nodes = nodes_df, edges = edges_df, raw_edges = edges)
    })

    output$network_plot <- visNetwork::renderVisNetwork({
      shiny::req(network_data())
      nd <- network_data()
      visNetwork::visNetwork(nd$nodes, nd$edges, height = "660px", width = "100%",
                             background = "#f2f6fa") %>%
        visNetwork::visEdges(
          smooth = list(type = "curvedCW", roundness = 0.15),
          arrows = list(to = list(enabled = TRUE, scaleFactor = 0.6))) %>%
        visNetwork::visNodes(font = list(size = 13, color = "#ffffff")) %>%
        visNetwork::visOptions(
          highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE),
          nodesIdSelection = list(enabled = TRUE)) %>%
        visNetwork::visInteraction(navigationButtons = TRUE, zoomView = TRUE) %>%
        visNetwork::visPhysics(
          solver = "forceAtlas2Based",
          forceAtlas2Based = list(gravitationalConstant = -60),
          stabilization = list(iterations = 200)) %>%
        visNetwork::visLayout(randomSeed = 42)
    })

    output$download_network_full <- shiny::downloadHandler(
      filename = "tf_network_full.csv",
      content  = function(file) utils::write.csv(full_network(), file, row.names = FALSE))
    output$download_network_filtered <- shiny::downloadHandler(
      filename = function() paste0("tf_network_filtered_", input$tf_center, ".csv"),
      content  = function(file) {
        shiny::req(network_data())
        utils::write.csv(network_data()$raw_edges, file, row.names = FALSE)
      })

    # ==========================================================
    # TAB 4 — Co-expression Modules
    # ==========================================================

    coexp_data <- shiny::eventReactive(input$run_coexp, {
      hub_df <- hdWGCNA::GetHubGenes(seurat_obj, n_hubs = input$coexp_n_hubs) %>%
        dplyr::arrange(.data$module, dplyr::desc(.data$kME))
      all_genes_df <- hdWGCNA::GetModules(seurat_obj) %>%
        dplyr::filter(.data$module != "grey") %>%
        dplyr::rowwise() %>%
        dplyr::mutate(kME_own = get(paste0("kME_", .data$module))) %>%
        dplyr::ungroup() %>%
        dplyr::select("gene_name", "module", "kME_own") %>%
        dplyr::arrange(.data$module, dplyr::desc(.data$kME_own))
      module_sizes <- all_genes_df %>%
        dplyr::count(.data$module, name = "n_genes") %>%
        dplyr::arrange(dplyr::desc(.data$n_genes))
      list(hub_df = hub_df, all_genes_df = all_genes_df,
           module_sizes = module_sizes,
           module_names = sort(unique(hub_df$module)))
    })

    output$module_size_table <- DT::renderDT({
      shiny::req(coexp_data())
      DT::datatable(coexp_data()$module_sizes, filter = "top",
                    options = list(pageLength = 10, scrollX = TRUE))
    })

    gene_lookup_result <- shiny::eventReactive(input$run_gene_lookup, {
      shiny::req(coexp_data(), nchar(trimws(input$gene_lookup)) > 0)
      g   <- trimws(input$gene_lookup)
      res <- coexp_data()$all_genes_df %>% dplyr::filter(.data$gene_name == g)
      if (nrow(res) == 0) paste0("Gene '", g, "' not found in any module.")
      else paste0("Gene:   ", g, "\nModule: ", res$module[1],
                  "\nkME:    ", round(res$kME_own[1], 4))
    })
    output$gene_lookup_result <- shiny::renderPrint({
      shiny::req(gene_lookup_result()); gene_lookup_result()
    })

    output$coexp_module_selector <- shiny::renderUI({
      shiny::req(coexp_data())
      shiny::selectInput("coexp_selected_module", "Module:",
                         choices  = coexp_data()$module_names,
                         selected = coexp_data()$module_names[1])
    })
    output$hub_genes_table <- DT::renderDT({
      shiny::req(coexp_data())
      DT::datatable(coexp_data()$hub_df, filter = "top",
                    options = list(pageLength = 15, scrollX = TRUE))
    })
    output$download_hub_genes <- shiny::downloadHandler(
      filename = "hub_genes_by_module.csv",
      content  = function(file) { shiny::req(coexp_data()); utils::write.csv(coexp_data()$hub_df, file, row.names = FALSE) })

    coexp_enrich_result <- shiny::eventReactive(input$run_coexp_enrich, {
      shiny::req(coexp_data(), input$coexp_selected_module)
      genes_in_module <- coexp_data()$all_genes_df %>%
        dplyr::filter(.data$module == input$coexp_selected_module) %>%
        dplyr::pull(.data$gene_name) %>% unique()
      shiny::validate(shiny::need(length(genes_in_module) >= 5,
                                  paste0("Too few genes in module '", input$coexp_selected_module, "'")))
      res <- run_enrichr_tidy(genes_in_module, input$coexp_enrich_db)
      shiny::validate(shiny::need(!is.null(res) && nrow(res) > 0, "No enrichment results."))
      utils::head(res, input$coexp_n_enrich)
    })

    output$enrich_module_title <- shiny::renderText({
      shiny::req(input$coexp_selected_module)
      paste0("Enrichment \u2014 ", input$coexp_selected_module, "  \u00b7  ", input$coexp_enrich_db)
    })

    coexp_enrich_plot_obj <- shiny::reactive({
      shiny::req(coexp_enrich_result())
      df        <- coexp_enrich_result()
      mod       <- input$coexp_selected_module
      pt        <- input$coexp_enrich_pt
      mod_color <- tryCatch({
        hdWGCNA::GetModules(seurat_obj) %>%
          dplyr::filter(.data$module == mod) %>%
          dplyr::pull(.data$color) %>% unique() %>% .[1]
      }, error = function(e) "#0fa3b1")
      if (is.na(mod_color) || length(mod_color) == 0) mod_color <- "#0fa3b1"
      ggplot2::ggplot(df, ggplot2::aes(
          x = reorder(.data$Term, .data$Combined.Score),
          y = .data$Combined.Score,
          fill = -log10(.data$Adjusted.P.value + 1e-300))) +
        ggplot2::geom_col(color = "white", linewidth = 0.25) +
        ggplot2::scale_fill_gradient(low = "lightyellow", high = mod_color, name = "-log10(adj.p)") +
        ggplot2::coord_flip() +
        ggplot2::theme_minimal(base_size = pt) +
        ggplot2::labs(title    = paste0("Module: ", mod),
                      subtitle = paste0(input$coexp_enrich_db, "  |  n genes = ",
                                        nrow(coexp_data()$all_genes_df %>% dplyr::filter(.data$module == mod))),
                      x = NULL, y = "Combined Score") +
        ggplot2::theme(plot.title    = ggplot2::element_text(face = "bold"),
                       plot.subtitle = ggplot2::element_text(color = "grey40"))
    })

    output$coexp_enrich_container <- shiny::renderUI({
      shiny::req(input$coexp_enrich_width, input$coexp_enrich_height)
      shiny::plotOutput("coexp_enrich_plot",
                        width  = paste0(input$coexp_enrich_width,  "px"),
                        height = paste0(input$coexp_enrich_height, "px"))
    })
    output$coexp_enrich_plot  <- shiny::renderPlot(coexp_enrich_plot_obj())
    output$coexp_enrich_table <- DT::renderDT({
      shiny::req(coexp_enrich_result())
      DT::datatable(
        coexp_enrich_result() %>%
          dplyr::select("Term","Overlap","P.value","Adjusted.P.value","Combined.Score","Genes") %>%
          dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~ signif(.x, 3))),
        filter = "top", options = list(pageLength = 8, scrollX = TRUE))
    })
    output$download_coexp_enrich <- shiny::downloadHandler(
      filename = function() paste0("enrichment_module_", input$coexp_selected_module, ".csv"),
      content  = function(file) { shiny::req(coexp_enrich_result()); utils::write.csv(coexp_enrich_result(), file, row.names = FALSE) })
    output$download_coexp_enrich_plot <- shiny::downloadHandler(
      filename = function() paste0("enrichment_barplot_", input$coexp_selected_module, ".png"),
      content  = function(file)
        ggplot2::ggsave(file, coexp_enrich_plot_obj(),
                        width = input$coexp_enrich_dl_w, height = input$coexp_enrich_dl_h,
                        dpi = input$coexp_enrich_dpi))

    # ==========================================================
    # TAB 5 — TF Regulon Heatmap + EnrichR
    # ==========================================================

    output$reg_module_selector <- shiny::renderUI({
      mods <- tryCatch({
        m <- unique(hdWGCNA::GetModules(seurat_obj)$module)
        c("All modules", sort(m[m != "grey"]))
      }, error = function(e) "All modules")
      shiny::selectInput("reg_selected_module", "Filter TFs by module:",
                         choices = mods, selected = "All modules")
    })

    reg_heatmap_data <- shiny::eventReactive(input$run_reg_heatmap, {
      shiny::withProgress(message = "Building TF regulon heatmap...", value = 0, {

        shiny::incProgress(0.05, detail = "Loading TF network...")
        net <- GetTFNetwork(seurat_obj) %>% dplyr::filter(.data$Gain >= input$reg_min_gain)
        shiny::validate(shiny::need(nrow(net) > 0, "No edges after Gain filter."))

        shiny::incProgress(0.15, detail = "Average expression per cluster...")
        avg_mat      <- Seurat::AverageExpression(seurat_obj, group.by = "Population_level3",
                                                  assays = Seurat::DefaultAssay(seurat_obj))[[1]]
        clusters     <- colnames(avg_mat)
        genes_in_mat <- rownames(avg_mat)

        shiny::incProgress(0.30, detail = "Filtering TFs...")
        module_tfs <- character(0)
        if (!is.null(input$reg_selected_module) && input$reg_selected_module != "All modules") {
          tryCatch({
            module_tfs <- hdWGCNA::GetModules(seurat_obj) %>%
              dplyr::filter(.data$module == input$reg_selected_module) %>%
              dplyr::pull(.data$gene_name) %>% unique()
          }, error = function(e) {})
        }
        manual_tfs <- if (nchar(trimws(input$reg_tf_filter)) > 0)
          trimws(strsplit(input$reg_tf_filter, ",\\s*")[[1]]) else character(0)

        all_tfs <- unique(net$tf)
        if (length(module_tfs) > 0) {
          all_tfs <- intersect(all_tfs, module_tfs)
          shiny::validate(shiny::need(length(all_tfs) > 0,
                                      paste0("No TFs from module '", input$reg_selected_module, "'.")))
        }
        tfs_valid <- all_tfs[all_tfs %in% genes_in_mat]
        min_sz    <- input$reg_min_regulon_size
        tfs_valid <- tfs_valid[vapply(tfs_valid, function(tf_name) {
          pos_n <- sum(net$tf == tf_name & net$Cor > 0 & net$gene %in% genes_in_mat)
          neg_n <- sum(net$tf == tf_name & net$Cor < 0 & net$gene %in% genes_in_mat)
          max(pos_n, neg_n) >= min_sz
        }, logical(1))]
        shiny::validate(shiny::need(length(tfs_valid) > 0, "No TFs pass filters."))

        if (length(manual_tfs) > 0) {
          tfs_use <- intersect(manual_tfs, tfs_valid)
          shiny::validate(shiny::need(length(tfs_use) > 0, "Specified TFs not in filtered network."))
        } else {
          reg_sizes <- vapply(tfs_valid, function(tf_name) sum(net$tf == tf_name), integer(1))
          tfs_use   <- tfs_valid[order(reg_sizes, decreasing = TRUE)][
            seq_len(min(input$reg_top_n_tfs, length(tfs_valid)))]
        }

        shiny::incProgress(0.50, detail = "Building expression matrix...")
        rows <- lapply(tfs_use, function(tf_name) {
          pos_genes <- net %>% dplyr::filter(.data$tf == tf_name, .data$Cor > 0) %>%
            dplyr::pull(.data$gene) %>% unique() %>% intersect(genes_in_mat)
          neg_genes <- net %>% dplyr::filter(.data$tf == tf_name, .data$Cor < 0) %>%
            dplyr::pull(.data$gene) %>% unique() %>% intersect(genes_in_mat)
          data.frame(
            tf      = tf_name, cluster = clusters,
            tf_expr = as.numeric(avg_mat[tf_name, ]),
            pos_reg = if (length(pos_genes) >= 1) as.numeric(colMeans(avg_mat[pos_genes, , drop=FALSE])) else rep(NA_real_, length(clusters)),
            neg_reg = if (length(neg_genes) >= 1) as.numeric(colMeans(avg_mat[neg_genes, , drop=FALSE])) else rep(NA_real_, length(clusters)),
            n_pos = length(pos_genes), n_neg = length(neg_genes), stringsAsFactors = FALSE
          )
        })
        df_wide <- dplyr::bind_rows(rows) %>% dplyr::group_by(.data$tf) %>%
          dplyr::mutate(tf_expr_scaled = {
            rng <- range(.data$tf_expr, na.rm = TRUE)
            if (diff(rng) == 0) rep(0.5, dplyr::n()) else (.data$tf_expr - rng[1]) / diff(rng)
          }) %>% dplyr::ungroup()
        df_long <- df_wide %>%
          tidyr::pivot_longer(cols = c("pos_reg","neg_reg"),
                              names_to = "regulon_dir", values_to = "mean_expr") %>%
          dplyr::mutate(
            regulon_dir = dplyr::case_when(
              .data$regulon_dir == "pos_reg" ~ "Positive regulon  (Cor > 0)",
              .data$regulon_dir == "neg_reg" ~ "Negative regulon  (Cor < 0)"),
            regulon_dir = factor(.data$regulon_dir,
              levels = c("Positive regulon  (Cor > 0)", "Negative regulon  (Cor < 0)"))
          )
        shiny::incProgress(1, detail = "Done.")
        list(df_long = df_long, df_wide = df_wide, tfs_use = tfs_use,
             clusters = clusters, module = input$reg_selected_module, net = net)
      })
    })

    output$reg_enrich_tf_selector <- shiny::renderUI({
      if (shiny::isTruthy(reg_heatmap_data())) {
        tfs <- reg_heatmap_data()$tfs_use
        if (length(tfs) > 0)
          shiny::selectInput("reg_enrich_tf", "TF for enrichment:", choices = tfs, selected = tfs[1])
        else shiny::helpText("Build the heatmap first.")
      } else shiny::helpText("Build the heatmap first.")
    })

    output$reg_status <- shiny::renderText({
      if (!shiny::isTruthy(reg_heatmap_data())) return("")
      d  <- reg_heatmap_data()
      mt <- if (!is.null(d$module) && d$module != "All modules") paste0("  |  Module: ", d$module) else ""
      paste0("Loaded: ", length(d$tfs_use), " TFs  \u00d7  ", length(d$clusters), " clusters", mt)
    })

    reg_heatmap_plot_obj <- shiny::reactive({
      shiny::req(reg_heatmap_data())
      df <- reg_heatmap_data()$df_long
      dmin <- input$reg_dot_range[1]; dmax <- input$reg_dot_range[2]; pt <- input$reg_pt
      tf_order <- df %>%
        dplyr::filter(.data$regulon_dir == "Positive regulon  (Cor > 0)") %>%
        dplyr::group_by(.data$tf) %>%
        dplyr::summarise(total_pos = sum(.data$mean_expr, na.rm=TRUE), .groups="drop") %>%
        dplyr::arrange(dplyr::desc(.data$total_pos)) %>% dplyr::pull(.data$tf)
      df <- df %>% dplyr::mutate(tf = factor(.data$tf, levels = rev(tf_order)))
      dot_theme <- ggplot2::theme_minimal(base_size = pt) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle=45, hjust=1, size=pt-1),
                       axis.text.y = ggplot2::element_text(size=pt-1),
                       panel.grid.major = ggplot2::element_line(color="grey92"),
                       legend.position = "right",
                       plot.title = ggplot2::element_text(face="bold", size=pt+1),
                       legend.key.size = ggplot2::unit(0.45,"cm"))
      make_dot <- function(data, fill_high, fill_name, title_txt) {
        ggplot2::ggplot(data, ggplot2::aes(x=.data$cluster, y=.data$tf,
                                            size=.data$tf_expr_scaled, fill=.data$mean_expr)) +
          ggplot2::geom_point(shape=21, color="grey30", stroke=0.3, alpha=0.9) +
          ggplot2::scale_size_continuous(range=c(dmin,dmax), name="TF expr\n(scaled 0\u20131)") +
          ggplot2::scale_fill_gradient(low="white", high=fill_high, na.value="grey85", name=fill_name) +
          ggplot2::labs(title=title_txt, x=NULL, y="Transcription Factor") + dot_theme
      }
      p_pos <- make_dot(df %>% dplyr::filter(.data$regulon_dir=="Positive regulon  (Cor > 0)"),
                        "firebrick","Mean expr\npos. regulon","Positive regulon  (Cor > 0)")
      p_neg <- make_dot(df %>% dplyr::filter(.data$regulon_dir=="Negative regulon  (Cor < 0)"),
                        "steelblue4","Mean expr\nneg. regulon","Negative regulon  (Cor < 0)") +
        ggplot2::labs(x="Cluster")
      title_suffix <- if (!is.null(reg_heatmap_data()$module) &&
                          reg_heatmap_data()$module != "All modules")
        paste0("  |  Module: ", reg_heatmap_data()$module) else ""
      p_pos / p_neg +
        patchwork::plot_annotation(
          title = paste0("TF Regulon Heatmap", title_suffix,
                         "  |  dot size = TF expr (scaled)  |  dot color = mean regulon expr"),
          theme = ggplot2::theme(plot.title = ggplot2::element_text(size=pt, color="grey40", hjust=0.5)))
    })

    output$reg_heatmap_container <- shiny::renderUI({
      shiny::req(input$reg_width, input$reg_height)
      shiny::plotOutput("reg_heatmap_plot",
                        width  = paste0(input$reg_width,  "px"),
                        height = paste0(input$reg_height, "px"))
    })
    output$reg_heatmap_plot  <- shiny::renderPlot({ shiny::req(reg_heatmap_plot_obj()); reg_heatmap_plot_obj() })
    output$reg_summary_table <- DT::renderDT({
      shiny::req(reg_heatmap_data())
      DT::datatable(
        reg_heatmap_data()$df_wide %>% dplyr::select("tf","n_pos","n_neg") %>%
          dplyr::distinct() %>% dplyr::mutate(total_regulon = .data$n_pos + .data$n_neg) %>%
          dplyr::arrange(dplyr::desc(.data$total_regulon)),
        filter = "top", options = list(pageLength = 15, scrollX = TRUE))
    })
    output$download_reg_plot <- shiny::downloadHandler(
      filename = function() {
        mod <- if (!is.null(input$reg_selected_module) && input$reg_selected_module != "All modules")
          paste0("_", input$reg_selected_module) else ""
        paste0("tf_regulon_heatmap", mod, ".png")
      },
      content = function(file) {
        shiny::req(reg_heatmap_plot_obj())
        n_tfs  <- length(reg_heatmap_data()$tfs_use)
        height <- max(input$reg_dl_h, n_tfs * 0.35 * 2)
        ggplot2::ggsave(file, reg_heatmap_plot_obj(), width=input$reg_dl_w,
                        height=height, dpi=input$reg_dpi, limitsize=FALSE)
      })
    output$download_reg_csv <- shiny::downloadHandler(
      filename = function() {
        mod <- if (!is.null(input$reg_selected_module) && input$reg_selected_module != "All modules")
          paste0("_", input$reg_selected_module) else ""
        paste0("tf_regulon_heatmap_data", mod, ".csv")
      },
      content = function(file) { shiny::req(reg_heatmap_data()); utils::write.csv(reg_heatmap_data()$df_wide, file, row.names=FALSE) })

    # Regulon EnrichR
    get_reg_genes <- function(sign = c("pos","neg")) {
      sign <- match.arg(sign)
      net <- if (shiny::isTruthy(reg_heatmap_data())) reg_heatmap_data()$net
             else GetTFNetwork(seurat_obj) %>% dplyr::filter(.data$Gain >= input$reg_min_gain)
      tf_sel <- input$reg_enrich_tf
      shiny::validate(shiny::need(!is.null(tf_sel) && tf_sel != "", "Build the heatmap first."))
      shiny::validate(shiny::need(tf_sel %in% net$tf, paste0("TF '", tf_sel, "' not in network.")))
      genes <- if (sign == "pos")
        net %>% dplyr::filter(.data$tf == tf_sel, .data$Cor > 0) %>% dplyr::pull(.data$gene) %>% unique()
      else
        net %>% dplyr::filter(.data$tf == tf_sel, .data$Cor < 0) %>% dplyr::pull(.data$gene) %>% unique()
      shiny::validate(shiny::need(length(genes) >= 5, paste0("Too few ", sign, " targets (n < 5).")))
      genes
    }

    reg_enrich_results <- shiny::eventReactive(input$run_reg_enrich, {
      db <- input$reg_enrich_db
      list(pos = run_enrichr_tidy(get_reg_genes("pos"), db),
           neg = run_enrichr_tidy(get_reg_genes("neg"), db),
           tf  = input$reg_enrich_tf, db = db)
    })

    reg_enrich_pos_df <- shiny::reactive({ shiny::req(reg_enrich_results()$pos); utils::head(reg_enrich_results()$pos, input$reg_enrich_n) })
    reg_enrich_neg_df <- shiny::reactive({ shiny::req(reg_enrich_results()$neg); utils::head(reg_enrich_results()$neg, input$reg_enrich_n) })

    make_reg_enrich_plot <- function(df, low_col, high_col, tf_name, db, direction, pt) {
      ggplot2::ggplot(df, ggplot2::aes(x=.data$Combined.Score, y=reorder(.data$Term,.data$Combined.Score), fill=.data$neg_log10_p)) +
        ggplot2::geom_col(color="white", linewidth=0.25) +
        ggplot2::scale_fill_gradient(low=low_col, high=high_col, name="-log10\n(adj.p)") +
        ggplot2::labs(title=paste0(direction," regulon \u2014 ",tf_name),
                      subtitle=paste0(db,"  |  top ",nrow(df)," terms"),
                      x="Combined Score", y=NULL) +
        ggplot2::theme_minimal(base_size=pt) +
        ggplot2::theme(plot.title=ggplot2::element_text(face="bold"),
                       plot.subtitle=ggplot2::element_text(color="grey40"))
    }

    reg_enrich_pos_plot_obj <- shiny::reactive({
      shiny::req(reg_enrich_pos_df())
      make_reg_enrich_plot(reg_enrich_pos_df(),"#f4a582","#ca0020",
                           reg_enrich_results()$tf, reg_enrich_results()$db,
                           "Positive", input$reg_enrich_pt)
    })
    reg_enrich_neg_plot_obj <- shiny::reactive({
      shiny::req(reg_enrich_neg_df())
      make_reg_enrich_plot(reg_enrich_neg_df(),"#92c5de","#0571b0",
                           reg_enrich_results()$tf, reg_enrich_results()$db,
                           "Negative", input$reg_enrich_pt)
    })

    enrich_dt_tidy <- function(df) {
      df %>% dplyr::select("Term","Overlap","Adjusted.P.value","Combined.Score","Genes") %>%
        dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~ signif(.x, 3)))
    }
    output$reg_enrich_pos_table <- DT::renderDT({ shiny::req(reg_enrich_pos_df());
      DT::datatable(enrich_dt_tidy(reg_enrich_pos_df()), filter="top", options=list(pageLength=8,scrollX=TRUE)) })
    output$reg_enrich_neg_table <- DT::renderDT({ shiny::req(reg_enrich_neg_df());
      DT::datatable(enrich_dt_tidy(reg_enrich_neg_df()), filter="top", options=list(pageLength=8,scrollX=TRUE)) })

    output$reg_enrich_pos_container <- shiny::renderUI({
      shiny::req(input$reg_enrich_width, input$reg_enrich_height)
      shiny::plotOutput("reg_enrich_pos_plot",
                        width=paste0(input$reg_enrich_width,"px"), height=paste0(input$reg_enrich_height,"px"))
    })
    output$reg_enrich_neg_container <- shiny::renderUI({
      shiny::req(input$reg_enrich_width, input$reg_enrich_height)
      shiny::plotOutput("reg_enrich_neg_plot",
                        width=paste0(input$reg_enrich_width,"px"), height=paste0(input$reg_enrich_height,"px"))
    })
    output$reg_enrich_pos_plot <- shiny::renderPlot(reg_enrich_pos_plot_obj())
    output$reg_enrich_neg_plot <- shiny::renderPlot(reg_enrich_neg_plot_obj())

    output$download_reg_enrich_pos_png <- shiny::downloadHandler(
      filename = function() paste0("POS_",input$reg_enrich_tf,"_",input$reg_enrich_db,".png"),
      content  = function(file) ggplot2::ggsave(file,reg_enrich_pos_plot_obj(),
        width=input$reg_enrich_dl_w,height=input$reg_enrich_dl_h,dpi=input$reg_enrich_dpi))
    output$download_reg_enrich_neg_png <- shiny::downloadHandler(
      filename = function() paste0("NEG_",input$reg_enrich_tf,"_",input$reg_enrich_db,".png"),
      content  = function(file) ggplot2::ggsave(file,reg_enrich_neg_plot_obj(),
        width=input$reg_enrich_dl_w,height=input$reg_enrich_dl_h,dpi=input$reg_enrich_dpi))
    output$download_reg_enrich_pos_csv <- shiny::downloadHandler(
      filename = function() paste0("POS_",input$reg_enrich_tf,"_",input$reg_enrich_db,".csv"),
      content  = function(file) utils::write.csv(reg_enrich_results()$pos, file, row.names=FALSE))
    output$download_reg_enrich_neg_csv <- shiny::downloadHandler(
      filename = function() paste0("NEG_",input$reg_enrich_tf,"_",input$reg_enrich_db,".csv"),
      content  = function(file) utils::write.csv(reg_enrich_results()$neg, file, row.names=FALSE))

  } # end server function
}
