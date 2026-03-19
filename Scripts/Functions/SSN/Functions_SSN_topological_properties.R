# Function to calculate the topological properties of SSNs



func_compute_SSN_topological_properties <- function( network, 
                                                     as_directed = TRUE, 
                                                     n_cores = 1 ){
  
  require(arrow)
  require(igraph)
  require(data.table)
  require(tidyverse)
  require(parallel)
  
  # Get list of samples
  sample_id_list <- network %>% select(!c("regulators", "targets")) %>% names()
  
  if (length(sample_id_list) == 0) {
    stop("Input network must contain at least one sample column besides 'regulators' and 'targets'")
  }
  
  
  #####
  
  
  net_topol_stats <- tryCatch(expr = {
    
    mclapply(X = sample_id_list, 
             mc.cores = n_cores,
             FUN = function(sample_select){
               
               # Subset the network for single sample
               network_tmp <- network %>%
                 select(c("regulators", "targets", all_of(sample_select))) %>% 
                 rename( "weight" = !!sample_select) %>% 
                 collect()
               
               # Remove loops
               network_tmp <- network_tmp[network_tmp$regulators != network_tmp$targets,]
               
               # Convert to igraph object 
               network_tmp <- graph_from_data_frame(d = network_tmp, 
                                                    directed = as_directed)
               
               remove_edges <- E(network_tmp)[is.na(E(network_tmp)$weight)]
               network_tmp <- delete_edges(graph = network_tmp, edges = remove_edges)
               network_tmp <- delete_edge_attr(graph = network_tmp, name = "weight")
               
               if(!as_directed){
                 network_tmp <- igraph::simplify(graph = network_tmp, remove.multiple = TRUE, remove.loops = TRUE)
               }
               
               # Calculate statistics
               deg_mode <- ifelse(as_directed, "out", "all")
               tmp1 <- data.frame( "Sample ID" = sample_select, 
                                   "Density" = edge_density(network_tmp, loops = FALSE), 
                                   "Number of connected components" = count_components(network_tmp), 
                                   "Largest component size" = max(components(network_tmp)$csize), 
                                   "Number of isolates" = sum(degree(network_tmp, mode = "total") == 0), 
                                   "Number of multi_edge" = sum(count_multiple(graph = as_undirected(graph = network_tmp, mode = "each")) > 1),
                                   "Power law exponent (directed)" = fit_power_law(degree(network_tmp, mode = deg_mode), xmin = 1)$alpha,
                                   "Power law exponent (undirected)" = fit_power_law(degree(as_undirected(graph = network_tmp, mode = "collapse"), mode = "all"), xmin = 1)$alpha,
                                   "Girth" = girth(network_tmp, circle = FALSE)$girth,
                                   "Reciprocity" = reciprocity(network_tmp, ignore.loops = TRUE, mode = "default"),
                                   "Transitivity" = transitivity(network_tmp, type = "global"),
                                   "Assortativity" = assortativity_degree(network_tmp, directed = as_directed),
                                   "Mean distance" = mean_distance(network_tmp, directed = as_directed, unconnected = TRUE, details = FALSE),
                                   "Diameter" = diameter(network_tmp, directed = as_directed, unconnected = TRUE),
              "Gini coefficient (directed)" = gini_vec(degree(network, mode = deg_mode)), 
              "Gini coefficient (undirected)" = gini_vec(degree(as_undirected(graph = network, mode = "collapse"), mode = "all")),
              "Hill index (directed)" = hill_index(degree(network, mode = deg_mode)),
              "Hill index (undirected)" = hill_index(degree(as_undirected(graph = network, mode = "collapse"), mode = "all")),
              check.names = FALSE
               )
               
               return(tmp1)
               
             })
  }, 
  warning = function(w){ warnings(w) },
  error = function(e){ stop(e) })
  
  net_topol_stats <- rbindlist(net_topol_stats)
  
  print(paste0(Sys.time(), " |  --- done!!"))
  
  return(net_topol_stats)
  
}



#####


func_plot_SSN_topological_properties <- function(topology_data){
  
  require(data.table)
  require(tidyverse)
  
  plot_data <- melt(topology_data, 
                    id.vars = "sample_id", 
                    variable.name = "metric", 
                    value.name = "score")
  
  plot <- ggplot(plot_data, aes(x = score)) +
    geom_freqpoly(linewidth = 0.1) +
    facet_wrap(~ metric, scale = "free") +
    labs(x = "Score", y = "Number of samples") +
    theme(panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.25, linetype = NULL),
          panel.grid = element_blank(),
          panel.spacing = unit(0.1, "cm"),
          text = element_text(size = 8),
          plot.title = element_text(face = "bold", hjust = 0.5, size = 8),
          axis.text.x = element_text(angle = 0, vjust = 0, hjust = 0.5),
          axis.ticks = element_line(colour = "black", linewidth = 0.2),
          legend.position = "bottom",
          legend.key = element_rect(fill = NA),
          legend.key.size = unit(0.4, "cm"),
          legend.key.spacing.y = unit(0.1, "cm"),
          legend.title = element_text(size = 6, face = "bold", margin = margin(t = 1, r = 0, b = 1.5, l = 0)),
          legend.text = element_text(size = 6, margin = margin(t = 0.5, r = 0, b = 0.1, l = 1)),
          legend.margin = margin(1,1,1,1),
          legend.spacing = unit(0.2, "cm")  )
  
  return(plot)
  
}





