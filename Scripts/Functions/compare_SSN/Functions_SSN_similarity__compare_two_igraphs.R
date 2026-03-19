# Function to compute similarity between two igraph objects


compare_two_igraph_networks <- function(igraph1, igraph2){
  
  set.seed(5081)
  require(igraph)
  
  if(!is_igraph(igraph1)){stop("Input network1 is not an igraph object")}
  if(!is_igraph(igraph2)){stop("Input network2 is not an igraph object")}
  
  if(is_directed(igraph1) != is_directed(igraph2)){
    stop("Input networks are of different types")
  }
  as_directed <- is_directed(igraph1)
  
  common_network <- igraph::intersection(igraph1, igraph2, keep.all.vertices = FALSE)
  union_network <- igraph::union(igraph1, igraph2)
  
  jaccard_nodes <- vcount(common_network) / vcount(union_network)
  jaccard_edges <- ecount(common_network) / ecount(union_network)
  
  deg_mode <- ifelse(as_directed, "out", "all")
  deg1 <- degree(igraph1, mode = deg_mode, normalized = TRUE)
  deg2 <- degree(igraph2, mode = deg_mode, normalized = TRUE)
  
  igraph1_topDegree <- names(deg1[deg1 >= quantile(unique(deg1), 0.95)])
  igraph2_topDegree <- names(deg2[deg2 >= quantile(unique(deg2), 0.95)])
  jaccard_topDegree <- length(intersect(igraph1_topDegree, igraph2_topDegree)) / length(union(igraph1_topDegree, igraph2_topDegree))
  
  igraph1_modules <- cluster_leiden(igraph1, objective_function = "modularity", resolution = 1, n_iterations = 10)
  igraph2_modules <- cluster_leiden(igraph2, objective_function = "modularity", resolution = 1, n_iterations = 10)
  adjrand_modules <- compare(igraph1_modules, igraph2_modules, method = "adjusted.rand")
  
  rm(list = c("igraph1", "igraph2", "common_network", "union_network", "deg1", "deg2", "igraph1_topDegree", "igraph2_topDegree", "igraph1_modules", "igraph2_modules"))
  
  tmp1 <- data.frame("Jaccard similarity (node)" = jaccard_nodes, 
                     "Jaccard similarity (edge)" = jaccard_edges, 
                     "Jaccard similarity (top degree)" = jaccard_topDegree, 
                     "Adj. Rand Index (communities)" = adjrand_modules, 
                     check.names = FALSE)
  
  return(tmp1)
  
}
