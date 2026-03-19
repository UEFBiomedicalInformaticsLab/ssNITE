# Function to get statistics of the prunned network

get_network_stats <- function(sample_select, network_file, prune_thresholds, as_directed){
  
  require(arrow)
  require(igraph)
  require(tidyverse)
  source("Scripts/Functions/Functions_diversity_evaluation.R")
  
  
  # Subset the network for single sample
  network <- open_dataset(network_file)
  network <- network %>%
    select(c("regulators", "targets", all_of(sample_select))) %>% 
    rename( "weight" = !!sample_select) %>% 
    collect()
  
  # Prune the network
  threshold <- as.numeric(prune_thresholds[sample_select])
  network$weight <- ifelse(network$weight > threshold, network$weight, NA)
  
  # Remove loops
  network <- network[network$regulators != network$targets,]
  
  # Convert to igraph object 
  network <- graph_from_data_frame(d = network, 
                                   directed = as_directed)
  
  remove_edges <- E(network)[is.na(E(network)$weight)]
  network <- delete_edges(graph = network, edges = remove_edges)
  network <- delete_edge_attr(graph = network, name = "weight")
  
  if(!as_directed){
    network <- igraph::simplify(graph = network, remove.multiple = TRUE, remove.loops = TRUE)
  }
  
  # Calculate statistics
  deg_mode <- ifelse(as_directed, "out", "all")
  
  degree_res <- degree(network, mode = deg_mode, normalized = FALSE)
  
  if(length(unique(degree_res)) > 1){
    fit1 <- fit_power_law(degree_res, xmin = NULL, force.continuous = FALSE, implementation = "plfit")
    fit2 <- WGCNA::scaleFreeFitIndex(degree_res)
  }else{
    fit1 <- list("alpha" = NA, "xmin" = NA)
    fit2 <- data.frame("Rsquared.SFT" = NA, "slope.SFT" = NA)
  }
  
  component_res <- components(network)
  largest_component <- names(component_res$membership[component_res$membership == which.max(component_res$csize)])
  largest_component <- induced_subgraph(graph = network, vids = V(network)[V(network)$name %in% largest_component])

  communities_res <- cluster_leiden(graph = as_undirected(network, 
                                                          mode = "collapse"),
                                    objective_function = "modularity", 
                                    resolution = 1, 
                                    n_iterations = 10)
  
  data.frame( "Sample ID" = sample_select, 
              
              "Number of nodes" = vcount(network),
              "Number of edges" = ecount(network),
              "Density" = edge_density(network, loops = FALSE), 
              "Connected component count" = component_res$no, 
              "Largest component size" = max(component_res$csize), 
              "Largest component density" = edge_density(largest_component), 
              "Number of isolates" = sum(degree(network, mode = "total") == 0), 
              "Number of multiedge" = sum(count_multiple(graph = as_undirected(graph = network, mode = "each")) > 1),
              "Mean degree" = mean(degree_res),
              
              "Girth" = girth(network, circle = FALSE)$girth,
              "Reciprocity" = reciprocity(network, ignore.loops = TRUE, mode = "default"),
              "Transitivity" = transitivity(network, type = "global"),
              "Assortativity" = assortativity_degree(network, directed = as_directed),
              "Mean distance" = mean_distance(network, directed = as_directed, unconnected = TRUE, details = FALSE),
              "Diameter" = diameter(network, directed = as_directed, unconnected = TRUE),
              
              "Modularity" = modularity(network, membership = membership(communities_res)),
              "Number of modules" = length(communities_res),
              
              "Power law exp." = fit1$alpha,
              "Power law xmin" = fit1$xmin,
              "Rsquared" = fit2[["Rsquared.SFT"]], 
              "Slope" = fit2[["slope.SFT"]], 
              
              "Gini coefficient" = gini_vec(degree_res), 
              "Hill index" = hill_index(degree_res),
              
              check.names = FALSE
  )
  
}