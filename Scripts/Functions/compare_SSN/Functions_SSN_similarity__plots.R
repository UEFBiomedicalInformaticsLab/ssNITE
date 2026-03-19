# Function to plot the similarity as heatmaps

plot_SSN_similarity <- function(similarity_data, sample_info){
  
  require(ComplexHeatmap)
  
  sample_info <- sample_info %>% column_to_rownames("sample_id")
  
  # Check if all samples have sample info
  sample_id_list <- sort(unique(c(similarity_data$Network1, similarity_data$Network2))) %>% as.character()
  if(!all(sample_id_list %in% row.names(sample_info))){
    stop("All samples from similarity data not represented in sample info")
  }
  sample_info <- sample_info[sample_id_list, ]
  
  
  # Create annotation for the heatmap
  top_annot <- HeatmapAnnotation(which = "column",
                                 df = sample_info,  
                                 simple_anno_size = unit(0.1, "cm"),
                                 annotation_name_gp = gpar(fontsize = 3),
                                 annotation_legend_param = list(title_gp = gpar(fontsize = 3),
                                                                labels_gp = gpar(fontsize = 2),
                                                                grid_height = unit(0.1, "cm"),
                                                                grid_width = unit(0.1, "cm"))
  )
  
  
  heatmaps <- list()
  for(metric_select in c("Jaccard similarity (node)", "Jaccard similarity (edge)", "Jaccard similarity (top degree)", "Adj. Rand Index (communities)")){
    
    # Prepare similarity matrix
    plot_data <- similarity_data %>% 
      select(c("Network1", "Network2", !!metric_select)) %>% 
      rename(similarity = !!metric_select)
    
    do_cluster <- ifelse(any(is.na(plot_data$similarity)), FALSE, TRUE)
    
    plot_data <-  plot_data %>%
      bind_rows(plot_data %>% 
                  rename(Network2 = Network1,
                         Network1 = Network2))
    
    plot_data <- plot_data %>% 
      bind_rows(data.frame("Network1" = sample_id_list, 
                           "Network2" = sample_id_list, 
                           similarity = 1))
    
    plot_data <- plot_data %>% 
      pivot_wider(names_from = "Network2", 
                  values_from = "similarity") %>%
      column_to_rownames("Network1") %>% 
      as.matrix()
    
    # Rearrange as in sample info 
    plot_data <- plot_data[row.names(sample_info), row.names(sample_info)]
    
    # Generate the plot
    plot <- Heatmap(plot_data,
                    
                    col = circlize::colorRamp2(breaks = c(0, 1), colors = c("white", "blue")),
                    cluster_rows  = do_cluster,
                    cluster_columns = do_cluster,
                    
                    row_title = "Network 1",
                    column_title = "Network 2",
                    row_title_side = "left",
                    column_title_side = "bottom",
                    row_title_gp = gpar(fontsize = 4, face = "bold"),
                    column_title_gp = gpar(fontsize = 4, face = "bold"),
                    
                    show_row_dend = FALSE,
                    show_column_dend = FALSE,
                    
                    top_annotation = top_annot,
                    
                    show_row_names = FALSE,
                    show_column_names = FALSE,
                    row_names_gp = gpar(fontsize = 1),
                    column_names_gp = gpar(fontsize = 1, linebreak = TRUE),
                    column_names_rot = 45, 
                    
                    heatmap_legend_param = list(title = "Similarity",
                                                title_gp = gpar(fontsize = 3),
                                                labels_gp = gpar(fontsize = 2),
                                                legend_height = unit(1, "cm"),
                                                legend_width = unit(0.01, "cm")
                    ))
    
    heatmaps[[metric_select]] <- ggplotify::as.ggplot(plot) + 
      labs(title = metric_select) + 
      theme( plot.title = element_text(face = "bold", hjust = 0.5, size = 4) )
    
    dev.off()
    
  }
  
  plot <- ggpubr::ggarrange(plotlist = heatmaps, 
                            nrow = 1, 
                            common.legend = FALSE)
  
  return(plot)
  
}

