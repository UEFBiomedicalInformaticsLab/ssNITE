# Script containing functions to calculate edge centralities not included in the igraph R package



#####


# Nearest neighbour edge centrality
# Ref: A straightforward edge centrality concept derived from generalizing degree and strength (https://doi.org/10.1038/s41598-022-08254-5)
# Currently implementing only the unweighed version
# Edges connecting very high degree nodes receive high score


nearest_neighbour_edge_centrality <- function(graph){
  
  if(is_directed(graph)){
    message("Directed graph detected. Collapsing to undirected graph.")
    graph <- as_undirected(graph, mode = "collapse")
  }
  
  if(!is_named(graph)){stop("The input graph must be named, i.e., the nodes should have names")}
  
  graph_degree <- degree(graph = graph, normalized = FALSE, mode = "all")
  graph_edges <- as_edgelist(graph)
  
  C_a <- unname(graph_degree[graph_edges[,1]])
  C_b <- unname(graph_degree[graph_edges[,2]])
  nnec <- (C_a + C_b - 2) / (abs(C_a - C_b) + 1) # Eq. 6 from the paper
  
  return(nnec)
  
}


#####
