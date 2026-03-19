# Function to compute centrality for one sample

compute_centralities_for_one_network <- function(sample_select, 
                                                 network_file, 
                                                 as_directed = FALSE,
                                                 IO_path = NULL){
  
  source("Scripts/Functions/SSN_centralities/Functions_SSN_centralities__additional_edge_centralities.R")
  
  message(paste0("Executing for sample: ", sample_select))
  
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
  rm(remove_edges)
  
  if(!as_directed){
    network <- simplify(graph = network, remove.multiple = TRUE, remove.loops = TRUE)
  }
  
  
  # Get the edge names
  if(as_directed){
    edge_names <- apply( as_edgelist(network), 1, function(x){ paste(x, collapse = "|") })
  }else{
    edge_names <- apply( as_edgelist(network), 1, function(x){ paste(sort(x), collapse = "|") })
  }
  
  centrality_res <- list()
  
  if(as_directed){
    centrality_res[["Degree (in)"]] <- degree(graph = network, normalized = TRUE, mode = "in")
    centrality_res[["Degree (out)"]] <- degree(graph = network, normalized = TRUE, mode = "out")
    centrality_res[["Degree"]] <- degree(graph = as_undirected(network, mode = "collapse"), normalized = TRUE, mode = "all")
    centrality_res[["Closeness (in)"]] <- closeness(graph = network, normalized = TRUE, mode = "in", weights = NA)
    centrality_res[["Closeness (out)"]] <- closeness(graph = network, normalized = TRUE, mode = "out", weights = NA)
    centrality_res[["Closeness"]] <- closeness(graph = as_undirected(network, mode = "collapse"), normalized = TRUE, mode = "all", weights = NA)
  }else{
    centrality_res[["Degree"]] <- degree(graph = network, normalized = TRUE, mode = "all")
    centrality_res[["Closeness"]] <- closeness(graph = network, normalized = TRUE, mode = "all", weights = NA)
  }
  
  centrality_res[["Betweenness"]] <- betweenness(graph = network, normalized = TRUE, directed = as_directed, weights = NA)
  centrality_res[["Clustering coef."]] <- transitivity(graph = network, type = "local", isolates = "zero", weights = NA)
  centrality_res[["Eigen centrality"]] <- eigen_centrality(network, directed = as_directed, weights = NA)$vector
  centrality_res[["Page rank"]] <- page_rank(graph = network, directed = as_directed, damping = 0.5, weights = NA)$vector
  
  centrality_res[["Edge betweenness"]] <- edge_betweenness(graph = network, directed = as_directed, weights = NA)
  centrality_res[["NNEC"]] <- nearest_neighbour_edge_centrality(graph = network)
  names(centrality_res[["Edge betweenness"]]) <-  names(centrality_res[["NNEC"]]) <- edge_names
  
  centrality_res <- lapply(centrality_res, function(x){
    x <- as.data.table(x, keep.rownames = "id")
    x[, sample_id := sample_select]
    setnames(x, "x", "centrality_score")
    return(x)
  })
  
  factor_levels <- names(centrality_res)
  centrality_res <- rbindlist(centrality_res, idcol = "centrality_type")
  centrality_res[, c("Node1", "Node2") := {
    x <- tstrsplit(id, "\\|", fill = NA_character_)
    if(length(x) == 1){ list(x[[1]], NA_character_) }else{ x }
  }]
  centrality_res[, id := NULL]
  
  factor_cols <- c("centrality_type", "sample_id", "Node1", "Node2")
  centrality_res[, (factor_cols) := lapply(.SD, as.factor), .SDcols = factor_cols]
  setcolorder(centrality_res, c("centrality_type", "sample_id", "Node1", "Node2", "centrality_score"))
  centrality_res[, centrality_type := factor(centrality_type, levels = factor_levels)]
  
  tmp_file <- paste0(IO_path, "tmp_SSN_centrality/centrality_scores_", sample_select, ".parquet")
  write_parquet(x = centrality_res, sink = tmp_file, compression = "lz4" )
  
  return(tmp_file)
  
}