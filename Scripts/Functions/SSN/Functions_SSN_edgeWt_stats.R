# Function calculate edge weight stats


func_plot_edge_weight_stats <- function(network,
                                        as_directed = TRUE,
                                        IO_path = NULL,
                                        n_cores = 10){
  
  require(parallel)
  require(arrow)
  require(tidyverse)
  
  
  # Check IO directory
  if(!dir.exists(paste0(IO_path, "/tmp_edge_density/"))){
    dir.create(paste0(IO_path, "/tmp_edge_density/"), recursive = TRUE)
  }
  
  # Check if the input network dataframe is correct
  if (!all(c("regulators", "targets") %in% names(network))) {
    stop("Input network must contain columns: 'regulators' and 'targets'")
  }
  
  # Get list of samples
  sample_id_list <- network %>% select(!c("regulators", "targets")) %>% names()
  
  if (length(sample_id_list) == 0) {
    stop("Input network must contain at least one sample column besides 'regulators' and 'targets'")
  }
  
  
  #####
  
  
  # Calculate the edge weight statistics
  
  print(paste0(Sys.time(), " | -- Calculating stats"))
  
  edgeWt_stats <- tryCatch(expr = {
    
    mclapply(X = sample_id_list, 
             mc.cores = n_cores,
             FUN = function(sample_select){
               
               # Subset the network for single sample
               network_mat <- network %>%
                 select(c("regulators", "targets", all_of(sample_select))) %>% 
                 rename( "weight" = !!sample_select) %>%
                 collect() %>% as.data.table()
               
               # Remove loops
               network_mat <- network_mat[network_mat$regulators != network_mat$targets,]
               
               # If using undirected network, keep only one edge
               if(!as_directed){
                 network_mat[, c("u", "v") := .(pmin(regulators, targets), pmax(regulators, targets))]
                 network_mat <- network_mat[, .(weight = max(weight, na.rm = TRUE)), by = .(u, v)]
                 setnames(network_mat, c("u", "v"), c("regulators", "targets"))
               }
               
               # Calculate basic stats
               valid_edgeWt <- network_mat[is.finite(weight) & weight > 0, weight]
               qs <- quantile(valid_edgeWt, probs = c(0.01, 0.25, 0.5, 0.75, 0.99), na.rm = TRUE)
               stats <- data.frame( "sample_id" = sample_select, 
                                    "edge_count" = nrow(network_mat),
                                    "prop_0" = round(sum(network_mat$weight == 0, na.rm = TRUE) / nrow(network_mat), 3), 
                                    "prop_NA" = round(sum(is.na(network_mat$weight)) / nrow(network_mat), 3), 
                                    "mean" = mean(valid_edgeWt, na.rm = TRUE),
                                    "sd" = sd(valid_edgeWt, na.rm = TRUE),
                                    "q01" = qs[1],
                                    "q25" = qs[2],
                                    "median" = qs[3],
                                    "q75" = qs[4],
                                    "q99" = qs[5]
               )
               
               return(stats)
               
             })
    
  }, 
  # warning = function(w){ warning(w) },
  error = function(e){ stop(e) })
  
  edgeWt_stats <- rbindlist(edgeWt_stats)
  gc()
  
  
  #####
  
  
  # Plot the edge weight stats
  
  print(paste0(Sys.time(), " | -- Plotting edge weight stats"))
  
  plot1 <- ggplot(edgeWt_stats, aes(x = sample_id)) +
    geom_line(aes(y = mean, 
                  group = 1), 
              color = "steelblue", 
              linewidth = 0.3) +
    geom_line(aes(y = median, group = 1), 
              color = "tomato", 
              linewidth = 0.3) +
    geom_ribbon(aes(ymin = q25, ymax = q75, group = 1), 
                fill = "darkgreen", 
                alpha = 0.2) +
    labs(x = "Samples", 
         y = "Edge weight") +
    theme( panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.25, linetype = NULL),
           panel.grid = element_blank(),
           panel.spacing = unit(0.1, "cm"),
           text = element_text(size = 4),
           plot.title = element_text(face = "bold", hjust = 0.5, size = 4),
           axis.text.x = element_blank(),
           axis.ticks = element_line(colour = "black", linewidth = 0.1),
           legend.position = "bottom",
           legend.key = element_rect(fill = NA),
           legend.key.size = unit(0.25, "cm"),
           legend.title = element_text(size = 3.5, face = "bold", margin = margin(t = 0, r = 0, b = 0, l = 1)),
           legend.text = element_text(size = 3, margin = margin(t = 0, r = 0, b = 0, l = 1)),
           legend.margin = margin(1,1,1,1),
           legend.spacing = unit(0, "cm") )
  
  
  #####
  
  
  # Compute the densities
  
  print(paste0(Sys.time(), " | -- Calculating density"))
  
  if(min(edgeWt_stats$q01) == max(edgeWt_stats$q01)){
    eps <- max(max(edgeWt_stats$q01) * 0.01, 1e-6) 
    grid <- seq(min(edgeWt_stats$q01) - eps, max(edgeWt_stats$q99) + eps, length.out = 1024)
  }else{
    grid <- seq(min(edgeWt_stats$q01), max(edgeWt_stats$q99), length.out = 1024)
  }
  
  edgeWt_density <- tryCatch(expr = {
    
    mclapply(X = sample_id_list, 
             mc.cores = n_cores,
             FUN = function(sample_select){
               
               # Subset the network for single sample
               network_mat <- network %>%
                 select(c("regulators", "targets", all_of(sample_select))) %>% 
                 rename( "weight" = !!sample_select) %>%
                 collect() %>% as.data.table()
               
               # Remove loops
               network_mat <- network_mat[network_mat$regulators != network_mat$targets,]
               
               # If using undirected network, keep only one edge
               if(!as_directed){
                 network_mat[, c("u", "v") := .(pmin(regulators, targets), pmax(regulators, targets))]
                 network_mat <- network_mat[, .(weight = max(weight, na.rm = TRUE)), by = .(u, v)]
                 setnames(network_mat, c("u", "v"), c("regulators", "targets"))
               }
               
               # Calculate the density
               valid_edgeWt <- network_mat[is.finite(weight) & weight > 0, weight]
               if(length(unique(valid_edgeWt)) == 1){
                 valid_edgeWt <- jitter(valid_edgeWt, amount = 1e-6) # Using to enable density estimation
               }
               kd <- KernSmooth::bkde(valid_edgeWt, gridsize = 1024)
               kd  <- approx(x = kd$x, y = kd$y, xout = grid, rule = 2)$y
               kd <- data.frame( "sample_id" = sample_select,
                                 "x" = grid,
                                 "y" = kd )
               rm(valid_edgeWt)
               
               fwrite(kd, file = paste0(IO_path, "/tmp_edge_density/", sample_select, ".csv.gz") )
               
               return(paste0(IO_path, "/tmp_edge_density/", sample_select, ".csv.gz"))
               
             })
    
  }, 
  # warning = function(w){ warning(w) },
  error = function(e){ stop(e) })
  
  edgeWt_density <- unlist(edgeWt_density)
  gc()
  
  
  #####
  
  
  # plot density
  
  print(paste0(Sys.time(), " | -- Plotting edge weight density"))
  
  edgeWt_density <- open_csv_dataset(edgeWt_density) %>% collect() 
  
  plot2 <- ggplot(edgeWt_density, aes(x = x, y = y, group = sample_id)) +
    geom_line(linewidth = 0.1)  +
    theme( panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.25, linetype = NULL),
           panel.grid = element_blank(),
           panel.spacing = unit(0.1, "cm"),
           text = element_text(size = 4),
           plot.title = element_text(face = "bold", hjust = 0.5, size = 4),
           # axis.text.x = element_text(angle = 0, vjust = 0, hjust = 0.5),
           axis.ticks = element_line(colour = "black", linewidth = 0.1),
           legend.position = "bottom",
           legend.key = element_rect(fill = NA),
           legend.key.size = unit(0.25, "cm"),
           legend.title = element_text(size = 3.5, face = "bold", margin = margin(t = 0, r = 0, b = 0, l = 1)),
           legend.text = element_text(size = 3, margin = margin(t = 0, r = 0, b = 0, l = 1)),
           legend.margin = margin(1,1,1,1),
           legend.spacing = unit(0, "cm") )
  
  unlink(paste0(IO_path, "/tmp_edge_density/"), recursive = TRUE)
  
  return(list("edgeWt_stats" = edgeWt_stats, 
              "edgeWt_stats_plot" = plot1,
              "edgeWt_density_plot" = plot2))
  
}
