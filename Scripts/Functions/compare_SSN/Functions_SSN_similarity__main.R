# Function to calculate the similarity between networks based on node composition,
# edge composition, top 5% nodes and communities

func_compute_SSN_similarity <- function(network_file, 
                                        sample_info,
                                        as_directed = TRUE, 
                                        IO_path = NULL, 
                                        n_cores = 10, 
                                        n_jobs = 1,
                                        use_slurm = FALSE){
  
  require(arrow)
  require(igraph)
  require(data.table)
  require(tidyverse)
  require(parallel)
  
  
  # Check if network file exists
  if(!file.exists(network_file)){
    stop(paste0("File does not exist: ", network_file))
  }
  
  # Create temporary directory to store networks
  if(is.null(IO_path)){IO_path = tempdir()}
  if(grepl("/$", IO_path)) { IO_path <- sub("/$", "", IO_path) }
  
  if(!dir.exists(paste0(IO_path, "/tmp_igraph_networks"))){
    dir.create(paste0(IO_path, "/tmp_igraph_networks"), recursive = TRUE)
  }
  if(!dir.exists(paste0(IO_path, "/tmp_networks_similarity"))){
    dir.create(paste0(IO_path, "/tmp_networks_similarity"), recursive = TRUE)
  }
  
  # Check if the input network dataframe is correct
  network <- open_dataset(network_file)
  if (!all(c("regulators", "targets") %in% names(network))) {
    stop("Input network must contain columns: 'regulators' and 'targets'")
  }  
  
  # Get the list of samples
  sample_id_list <- network %>% select(!c("regulators", "targets")) %>% names()
  
  
  #####
  
  
  # Create igraph objects
  
  print(paste0(Sys.time(), " |   - Creating igraph objects for each network"))
  
  source("Scripts/Functions/compare_SSN/Functions_SSN_similarity__create_igraph.R")
  
  if(use_slurm){
    
    source("Scripts/batch_tools/batchtools_lapply.R")
    
    file_list <- tryCatch(expr = {
      
      batch_lapply( FUN = create_igraph_network, 
                    FUN_args = list("sample_select" = sample_id_list), 
                    FUN_more_args = list("network_file" = network_file, 
                                         "as_directed" = as_directed,
                                         "IO_path" = IO_path), 
                    IO_path = IO_path, 
                    registry_name = "batchtools_registry__create_igraph",
                    packages = c("tidyverse", "arrow", "data.table", "igraph"), 
                    cluster_template = "Scripts/batch_tools/slurm_batchtools_config.tmpl", 
                    slurm_resources = list(partition = "small",
                                           time = "0-00:30:00",
                                           memory = "10G",
                                           ntasks = 1,
                                           exclude = "kudos9,kudos8"), 
                    max_concurrent_jobs = n_cores )
      
    }, 
    warning = function(w){ warning(w) },
    error = function(e){ stop(e) })
    
  }else{
    
    file_list <- tryCatch(expr = {
      
      mclapply(X = sample_id_list,
               mc.cores = n_cores,
               FUN = create_igraph_network,
               network_file = network_file,
               as_directed = as_directed,
               IO_path = IO_path)
      
    },
    warning = function(w){ warning(w) },
    error = function(e){ stop(e) })
    
  }
  
  
  #####
  
  
  # Compare the networks
  
  print(paste0(Sys.time(), " |   - Comparing the networks"))
  
  source("Scripts/Functions/compare_SSN/Functions_SSN_similarity__compute_similarity.R")
  
  if(use_slurm){
    
    source("Scripts/batch_tools/batchtools_lapply.R")
    
    file_list <- tryCatch(expr = {
      
      batch_lapply( FUN = compute_similarity, 
                    FUN_args = list("sample_idx" = 1:(length(sample_id_list)-1)), 
                    FUN_more_args = list("sample_id_list" = sample_id_list, 
                                         "IO_path" = IO_path, 
                                         "n_jobs" = n_jobs), 
                    IO_path = IO_path, 
                    registry_name = "batchtools_registry__compute_SSN_similarity",
                    packages = c("tidyverse", "arrow", "data.table", "parallel", "igraph"), 
                    cluster_template = "Scripts/batch_tools/slurm_batchtools_config.tmpl", 
                    slurm_resources = list(partition = "small",
                                           time = "0-02:00:00",
                                           memory = "30G",
                                           ntasks = n_jobs,
                                           exclude = "kudos9,kudos8"), 
                    max_concurrent_jobs = n_cores )
      
    }, 
    warning = function(w){ warning(w) },
    error = function(e){ stop(e) })
    
  }else{
    
    require(parallel)
    
    file_list <- tryCatch(expr = {
      
      mclapply(X = 1:(length(sample_id_list)-1), 
               mc.cores = n_cores, 
               FUN = compute_similarity, 
               sample_id_list = sample_id_list, 
               IO_path = IO_path, 
               n_jobs = n_jobs)
      
    }, 
    warning = function(w){ warning(w) },
    error = function(e){ stop(e) })
    
  }
  
  
  # Compile the results
  compare_res <- lapply(file_list, fread)
  compare_res <- rbindlist(compare_res)
  
  system(paste0("rm -r ", IO_path, "/tmp_igraph_networks/"))
  system(paste0("rm -r ", IO_path, "/tmp_networks_similarity/"))
  
  
  #####
  
  
  # Plot the similarity between the SSN as heatmaps
  
  print(paste0(Sys.time(), " |   - Plotting the similarity as heatmap"))
  
  source("Scripts/Functions/compare_SSN/Functions_SSN_similarity__plots.R")
  
  compare_plot <- plot_SSN_similarity(compare_res, sample_info)
  
  print(paste0(Sys.time(), " |   --- done !!!"))
  
  return(list("similarity_data" = compare_res, 
              "similarity_plot" = compare_plot))
  
}
