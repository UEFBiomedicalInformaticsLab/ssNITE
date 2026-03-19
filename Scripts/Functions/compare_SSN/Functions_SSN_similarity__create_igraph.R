# Function to create igraph objects

create_igraph_network <- function(sample_select, network_file, as_directed, IO_path){
  
  network <- open_dataset(network_file) %>%
    select(c("regulators", "targets", all_of(sample_select))) %>%
    rename( "weight" = !!sample_select) %>%
    collect()
  
  network <- network[network$regulators != network$targets,]
  
  network <- graph_from_data_frame(d = network,
                                       directed = as_directed)
  
  remove_edges <- E(network)[is.na(E(network)$weight)]
  network <- delete_edges(graph = network, edges = remove_edges)
  network <- delete_edge_attr(graph = network, name = "weight")
  
  if(!as_directed){
    network <- simplify(graph = network, remove.multiple = TRUE, remove.loops = TRUE)
  }
  
  tmp_file <- paste0(IO_path, "/tmp_igraph_networks/", sample_select, ".rds")
  
  saveRDS(object = network, file = tmp_file)
  
  return(tmp_file)
  
}