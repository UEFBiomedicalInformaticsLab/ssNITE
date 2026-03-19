# Function to create single-sample network using SWEET


create_SWEET_net <- function( gene_expr, 
                              IO_path,
                              sweet_k = 0.1,
                              sweet_sparsify = FALSE,
                              n_cores = 10, 
                              n_jobs = 1, 
                              use_slurm = FALSE ){
  
  require(arrow)
  require(data.table)
  require(tidyverse)
  require(parallel)
  
  
  #####
  
  
  # Input checks
  if(!dir.exists(IO_path)){dir.create(IO_path, recursive = TRUE)}
  if(grepl("/$", IO_path)) { IO_path <- sub("/$", "", IO_path) }
  
  
  # Create the directory to store results
  if(!dir.exists(paste0(IO_path, "/SWEET_results_tmp/"))){dir.create(paste0(IO_path, "/SWEET_results_tmp/"), recursive = TRUE)}
  if(!dir.exists(paste0(IO_path, "/SWEET_edges_tmp/"))){dir.create(paste0(IO_path, "/SWEET_edges_tmp/"), recursive = TRUE)}
  
  
  # Get the list of sample IDs
  sample_id_list <- colnames(gene_expr)
  
  
  ######
  
  
  # Calculate sample weights
  print(paste0(Sys.time(), " |   - Calculating sample weights"))
  sample_cor <- stats::cor(as.matrix(gene_expr), method = "pearson")
  mean_cor <- apply(sample_cor, 2, function(x){ (sum(x) - 1) / (length(sample_id_list) - 1) })
  sample_weights <- (mean_cor - min(mean_cor) + 0.01 ) / ( max(mean_cor) - min(mean_cor) + 0.01 )
  sample_weights <- sample_weights * sweet_k * length(sample_id_list)
  rm(list = c("mean_cor", "sample_cor"))
  
  
  # Generate the reference network
  print(paste0(Sys.time(), " |   - Generating reference network"))
  ref_net_cor <- stats::cor(t(gene_expr), method = "pearson")
  
  
  #####
  
  
  # Calculate the edge scores for the individual networks
  
  print(paste0(Sys.time(), " |   - Calculating edge scores"))
  
  source("Scripts/Functions/SWEET/Functions_SWEET__calculate_edge_scores.R")
  
  file_list <- tryCatch(expr = {
    
    mclapply( X = sample_id_list, 
              mc.cores = n_cores, 
              FUN = calculate_edge_score, 
              gene_expr = gene_expr, 
              ref_net_cor = ref_net_cor,
              sample_weights = sample_weights, 
              IO_path = IO_path,
              n_jobs = n_jobs )
    
  },
  warning = function(w){ warning(w) },
  error = function(e){ stop(e) } )
  
  print(paste0(Sys.time(), " |   --- done!!"))
  gc()
  
  
  #####
  
  
  # Calculate mean and SD for z-score
  
  if(sweet_sparsify){
    print(paste0(Sys.time(), " |   - Calculating edge weight stats"))
    edge_stats <- open_dataset(unlist(file_list))
    edge_stats <- edge_stats %>% 
      summarise( mean_val = mean(raw_edge_score, na.rm = TRUE),
                 sd_val = sd(raw_edge_score, na.rm = TRUE) ) %>% 
      collect()
  }else{
    edge_stats <- data.frame("mean_val" = NA, "sd_val" = NA)
  }
  
  
  #####
  
  
  # Export edges
  
  print(paste0(Sys.time(), " |   - Exporting edges"))
  
  source("Scripts/Functions/SWEET/Functions_SWEET__export_edges.R")
  
  if(use_slurm){
    
    source("Scripts/batch_tools/batchtools_lapply.R")
    
    file_list <- tryCatch(expr = {
      
      batch_lapply( FUN = export_SWEET_edges, 
                    FUN_args = list("gene" = row.names(gene_expr)), 
                    FUN_more_args = list("sample_id_list" = sample_id_list, 
                                         "mean_edgeWt" = edge_stats$mean_val,
                                         "sd_edgeWt" = edge_stats$sd_val,  
                                         "IO_path" = IO_path, 
                                         "n_jobs" = n_jobs), 
                    IO_path = IO_path, 
                    registry_name = "batchtools_registry__export_edges",
                    packages = c("tidyverse", "arrow", "data.table", "parallel"), 
                    cluster_template = "Scripts/batch_tools/slurm_batchtools_config.tmpl", 
                    slurm_resources = list(partition = "small",
                                           time = "0-01:00:00",
                                           memory = "25G",
                                           ntasks = n_jobs,
                                           exclude = "kudos12,kudos9,kudos8"), 
                    max_concurrent_jobs = n_cores )
      
    }, 
    warning = function(w){ warning(w) },
    error = function(e){ stop(e) }
    )
    
  }else{
    
    require(parallel)
    
    file_list <- tryCatch(expr = {
      
      mclapply( X = row.names(gene_expr),
                mc.cores = n_cores,
                FUN = export_SWEET_edges,
                sample_id_list = sample_id_list, 
                mean_edgeWt = edge_stats$mean_val,
                sd_edgeWt = edge_stats$sd_val, 
                IO_path = IO_path,
                n_jobs = n_jobs )
      
    },
    warning = function(w){ warning(w) },
    error = function(e){ stop(e) }
    )
    
  }
  
  print(paste0(Sys.time(), " |   --- done!!"))
  gc()
  
  
  #####
  
  
  print(paste0(Sys.time(), " |   - Exporting final network"))
  network <- open_dataset(unlist(file_list)) 
  system(paste0("rm -r ", IO_path, "/SWEET_results_tmp/"))
  
  print(paste0(Sys.time(), " |   --- done!!"))
  
  return(network)
  
}


