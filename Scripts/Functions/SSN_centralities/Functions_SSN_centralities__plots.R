# Function to cluster samples by centralities

cluster_SSN_by_node_centrality <- function(network_centralities, 
                                           sample_info, 
                                           ggplot_options = list( "shape" = NULL, 
                                                                  "color" = NULL,
                                                                  "manual_color" = NULL,
                                                                  "manual_shape" = NULL ), 
                                           IO_path = NULL){
  
  if(!dir.exists(paste0(IO_path, "pca_res"))){
    dir.create(paste0(IO_path, "pca_res"), recursive = TRUE)
  }
  
  plot_theme <- theme(panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.25, linetype = NULL),
                      panel.grid = element_blank(),
                      panel.spacing = unit(0.1, "cm"),
                      text = element_text(size = 4),
                      plot.title = element_text(face = "bold", hjust = 0.5, size = 4),
                      axis.text.x = element_text(angle = 0, vjust = 0, hjust = 0.5),
                      axis.ticks = element_line(colour = "black", linewidth = 0.2),
                      legend.position = "bottom",
                      legend.key = element_rect(fill = NA),
                      legend.key.size = unit(0.25, "cm"),
                      legend.key.spacing.y = unit(0.1, "cm"),
                      legend.title = element_text(size = 3, face = "bold", margin = margin(t = 1, r = 1.5, b = 1, l = 1.5)),
                      legend.text = element_text(size = 3, margin = margin(t = 1, r = 0, b = 1, l = 0.5)),
                      legend.margin = margin(1,1,1,1),
                      legend.spacing = unit(0, "cm")  )
  
  plot_list <- list()
  
  centrality_type_list <- network_centralities %>% pull(centrality_type, as_vector = TRUE) %>% levels()
  
  for(select_type in centrality_type_list){
    
    print(paste0(Sys.time(), " |   --- using: ", select_type))
    
    # plot_data <- network_centralities %>% filter(centrality_type %in% select_type) %>% collect() %>% as.data.table()
    # plot_data[, centrality_type := NULL]
    
    plot_data <- network_centralities %>% 
      filter(centrality_type %in% select_type) %>% 
      select(!c("centrality_type", "Node1", "Node2")) %>% 
      collect() %>% 
      as.matrix()
    
    
    if( (is_empty(plot_data)) | (nrow(plot_data) == 0) | (all(is.na(plot_data))) ){ #[, !c("Node1", "Node2")]
      
      message("'plot_data' is empty!!! Empty plot will be generated.")
      
      pca_res <- matrix( data = NA, nrow = length(colnames(plot_data)) , ncol = 10, dimnames = list(colnames(plot_data), paste0("PC", 1:10)) )
      write.csv( pca_res %>% as.data.table(keep.rownames = "sample_id") %>% mutate("centrality_type" = select_type), 
                 file = paste0(IO_path, "pca_res/PCA_", select_type, ".csv"), row.names = FALSE)
      rm(pca_res)
      
      plot <- ggplot() +
        labs(title = select_type, 
             x = paste0("PC1: ", "---", "% variance"),
             y = paste0("PC2: ", "---" , "% variance")) + 
        plot_theme
      
    }else{
      
      # plot_data <- t(plot_data[, !c("Node1", "Node2")])
      plot_data <- t(plot_data)
      plot_data[is.na(plot_data)] <- 0 # NA in closeness replaced by 0
      
      zero_var_col <- round(apply(plot_data, 2, var, na.rm = TRUE), 12)
      zero_var_col <- which(zero_var_col == 0)
      
      if(length(zero_var_col) == ncol(plot_data)){
        
        message("Whole data with zero variance. Empty plot will be generated.")
        
        pca_res <- matrix( data = NA, nrow = length(row.names(plot_data)) , ncol = 10, dimnames = list(row.names(plot_data), paste0("PC", 1:10)) )
        write.csv( pca_res %>% as.data.table(keep.rownames = "sample_id") %>% mutate("centrality_type" = select_type), 
                   file = paste0(IO_path, "pca_res/PCA_", select_type, ".csv"), row.names = FALSE)
        rm(pca_res)
        
        plot <- ggplot() +
          labs(title = select_type, 
               x = paste0("PC1: ", "---", "% variance"),
               y = paste0("PC2: ", "---" , "% variance")) + 
          plot_theme
        
      }else{
        
        if(length(zero_var_col) > 0){
          plot_data <- plot_data[,-zero_var_col] # Remove 0 variance columns
        }
        
        if(select_type %in% c("Edge betweenness", "NNEC")){
          pca_res <- irlba::prcomp_irlba(plot_data, n = 10, center = TRUE, scale. = TRUE) # using faster implementation 
          row.names(pca_res$x) <- row.names(plot_data)
        }else{
          pca_res <- stats::prcomp(plot_data, center = TRUE, scale. = TRUE) # using exact calculation
        }
        
        write.csv(pca_res$x[,1:10] %>% as.data.table(keep.rownames = "sample_id") %>% mutate("centrality_type" = select_type), 
                  file = paste0(IO_path, "pca_res/PCA_", select_type, ".csv"), row.names = FALSE)
        
        pca_eigen <- factoextra::get_eig(pca_res)
        
        pca_res <- as_tibble(pca_res$x, rownames = "sample_id") %>% 
          select(c("sample_id", "PC1", "PC2")) %>% 
          left_join(sample_info, by = "sample_id")
        
        aes_mapping <- aes(x = PC1, y = PC2)
        
        if (!is.null(ggplot_options$shape)) {
          aes_mapping <-  modifyList(aes_mapping, aes(shape = !!sym(ggplot_options$shape)))
        }
        
        if (!is.null(ggplot_options$color)) {
          aes_mapping <-  modifyList(aes_mapping, aes(color = !!sym(ggplot_options$color)))
        }
        
        
        plot <- ggplot() +
          geom_point(data = pca_res,
                     mapping = aes_mapping,
                     size = 1,
                     stroke = 0.25) +
          labs(title = select_type, 
               x = paste0("PC1: ", round(pca_eigen$variance.percent[1], 2), "% variance"),
               y = paste0("PC2: ", round(pca_eigen$variance.percent[2], 2) , "% variance"),
               shape = ggplot_options$shape,
               color = ggplot_options$color) + 
          plot_theme
        
        if(!is.null(ggplot_options$manual_color)){
          plot <- plot + scale_color_manual(values = ggplot_options$manual_color, na.value = "#A9A9A9")
        }
        
        if(!is.null(ggplot_options$manual_shape)){
          plot <- plot + scale_shape_manual(values = ggplot_options$manual_shape, na.value = 0)
        }
        
      }
      
    }
    
    plot_list[[select_type]] <- plot
    
    rm(list = c("plot", "plot_data"))
    gc()
    
  }
  
  pca_res <- open_csv_dataset(paste0(IO_path, "pca_res")) %>% 
    select(c("centrality_type", "sample_id", everything())) %>% 
    collect() %>% as.data.table()
  pca_res[, c("centrality_type", "sample_id") := lapply(.SD, as.factor), .SDcols = c("centrality_type", "sample_id")]
  saveRDS(pca_res, paste0(IO_path, "pca_res.rds"))
  rm(pca_res)
  unlink(paste0(IO_path, "pca_res"), recursive = TRUE)
  
  plot <- ggpubr::ggarrange(plotlist = plot_list, 
                            ncol = 4,
                            nrow = 2,
                            common.legend = TRUE, 
                            legend = "bottom")
  
  return(list("pca_res_file" = paste0(IO_path, "pca_res.rds"), 
              "plot" = plot))
  
}
