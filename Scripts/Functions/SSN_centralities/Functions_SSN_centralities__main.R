# Functions to compare SSNs by network centralities

func_compute_SSN_centralities <- function(network_file, 
                                          as_directed = TRUE, 
                                          sample_info = NULL,
                                          ggplot_options = list( "shape" = NULL, 
                                                                 "color" = NULL,
                                                                 "manual_color" = NULL,
                                                                 "manual_shape" = NULL ),
                                          n_cores = 10, 
                                          IO_path = NULL,
                                          use_slurm = FALSE){
  
  require(arrow)
  require(igraph)
  require(data.table)
  require(tidyverse)
  
  # Check if network file exists
  if(!file.exists(network_file)){
    stop(paste0("File does not exist: ", network_file))
  }
  
  # Check if the input network dataframe is correct
  network <- open_dataset(network_file)
  if (!all(c("regulators", "targets") %in% names(network))) {
    stop("Input network must contain columns: 'regulators' and 'targets'")
  }  
  
  # Get the list of samples
  sample_id_list <- network %>% select(!c("regulators", "targets")) %>% names()
  
  # Check the sample info 
  if( !("sample_id" %in% colnames(sample_info)) ){
    stop("'sample_info' must contain the column: sample_id")
  }
  
  # Check the ggplot options
  if( !is_null(ggplot_options$shape) & !(ggplot_options$shape %in% colnames(sample_info)) ){
    stop(paste0("'sample_info' must contain the column: ", ggplot_options$shape))
  }
  if( !is_null(ggplot_options$color) & !(ggplot_options$color %in% colnames(sample_info)) ){
    stop(paste0("'sample_info' must contain the column: ", ggplot_options$color))
  }
  
  
  #####
  
  
  # Calculate network centralities
  
  print(paste0(Sys.time(), " |   - calculating centrality"))
  
  source("Scripts/Functions/SSN_centralities/Functions_SSN_centralities__compute_centralities.R")
  
  if(!dir.exists( paste0(IO_path, "tmp_SSN_centrality") )){
    dir.create( paste0(IO_path, "tmp_SSN_centrality"), recursive = TRUE )
  }
  
  
  if(use_slurm){
    
    source("Scripts/batch_tools/batchtools_lapply.R")
    
    file_list <- tryCatch(expr = {
      
      batch_lapply( FUN = compute_centralities_for_one_network,
                    FUN_args = list("sample_select" = sample_id_list),
                    FUN_more_args = list("network_file" = network_file,
                                         "as_directed" = as_directed,
                                         "IO_path" = IO_path),
                    IO_path = IO_path,
                    registry_name = "batchtools_registry__compute_SSN_centralities",
                    packages = c("tidyverse", "arrow", "igraph", "data.table"),
                    cluster_template = "Scripts/batch_tools/slurm_batchtools_config.tmpl",
                    slurm_resources = list(partition = "small",
                                           time = "0-00:30:00",
                                           memory = "20GB",
                                           ntasks = 1,
                                           exclude = "kudos9,kudos8"),
                    max_concurrent_jobs = n_cores )
      
    },
    # warning = function(w){ warning(w) },
    error = function(e){ stop(e) })
    
  }else{
    
    require(parallel)
    
    file_list <- tryCatch(expr = {
      
      mclapply(X = sample_id_list,
               mc.cores = n_cores,
               FUN = compute_centralities_for_one_network,
               network_file = network_file,
               as_directed = as_directed,
               IO_path = IO_path)
      
    },
    # warning = function(w){ warning(w) },
    error = function(e){ stop(e) })
    
  }
  
  network_centralities <- open_dataset(unlist(file_list))
  write_parquet(x = network_centralities,
                sink = paste0(IO_path, "SSN_centrality_scores.parquet"),
                compression = "lz4" )
  
  unlink(paste0(IO_path, "tmp_SSN_centrality"), recursive = TRUE, force = TRUE)
  rm(list = c("network_centralities", "file_list"))
  gc()
  
  
  #####
  
  
  # Collect the centralities as dataframes
  
  print(paste0(Sys.time(), " |   - compiling centralities"))
  
  source("Scripts/Functions/SSN_centralities/Functions_SSN_centralities__process_centrality_table.R")
  
  if(!dir.exists( paste0(IO_path, "tmp_process_table") )){
    dir.create( paste0(IO_path, "tmp_process_table"), recursive = TRUE )
  }
  
  network_centralities <- open_dataset(paste0(IO_path, "SSN_centrality_scores.parquet"))
  node_list <- network_centralities %>% pull(Node1, as_vector = TRUE) %>% levels()
  rm(list = c("network_centralities"))
  
  
  if(use_slurm){
    
    source("Scripts/batch_tools/batchtools_lapply.R")
    
    file_list <- tryCatch(expr = {
      
      batch_lapply( FUN = process_centrality_table,
                    FUN_args = list("node" = node_list),
                    FUN_more_args = list("centrality_table_file" = paste0(IO_path, "SSN_centrality_scores.parquet"),
                                         "IO_path" = IO_path),
                    IO_path = IO_path,
                    registry_name = "batchtools_registry__process_SSN_centralities",
                    packages = c("tidyverse", "arrow", "data.table"),
                    cluster_template = "Scripts/batch_tools/slurm_batchtools_config.tmpl",
                    slurm_resources = list(partition = "small",
                                           time = "0-00:30:00",
                                           memory = "20GB",
                                           ntasks = 1,
                                           exclude = "kudos9,kudos8"),
                    max_concurrent_jobs = n_cores )
      
    },
    # warning = function(w){ warning(w) },
    error = function(e){ stop(e) })
    
  }else{
    
    require(parallel)
    
    file_list <- tryCatch(expr = {
      
      mclapply(X = node_list,
               mc.cores = n_cores,
               FUN = process_centrality_table,
               centrality_table_file = paste0(IO_path, "SSN_centrality_scores.parquet"),
               IO_path = IO_path)
      
    },
    # warning = function(w){ warning(w) },
    error = function(e){ stop(e) })
    
  }
  
  
  network_centralities <- open_dataset(unlist(file_list))
  
  write_parquet(x = network_centralities,
                sink = paste0(IO_path, "SSN_centrality_scores.parquet"),
                compression = "lz4" )
  
  unlink(paste0(IO_path, "tmp_process_table"), recursive = TRUE, force = TRUE)
  rm(list = c("node_list", "network_centralities", "file_list"))
  gc()
  
  
  #####
  
  
  # Cluster SSNs based on the centralities
  
  print(paste0(Sys.time(), " |   - clustering samples"))
  
  source("Scripts/Functions/SSN_centralities/Functions_SSN_centralities__plots.R")
  
  network_centralities <- open_dataset(paste0(IO_path, "SSN_centrality_scores.parquet"))
  sample_cluster_plot <- cluster_SSN_by_node_centrality(network_centralities,
                                                        sample_info = sample_info,
                                                        ggplot_options = ggplot_options, 
                                                        IO_path = IO_path)
  
  
  #####
  
  
  print(paste0(Sys.time(), " |   - exporting results"))
  
  return( list( "centrality_tables" = paste0(IO_path, "SSN_centrality_scores.parquet"),
                "pca_res_file" = sample_cluster_plot$pca_res_file,
                "sample_cluster_plot" = sample_cluster_plot$plot ) )
  
  
}
