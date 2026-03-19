# Function to create single-sample network using CSN

# Executes the csnet() function from https://github.com/wys8c764/CSN
# Uses MATLAB for implementation


create_CSN_net <- function( gene_expr,
                            csn_weighted = TRUE, 
                            IO_path,
                            n_cores = 10,
                            use_slurm = FALSE ){
  
  require(R.matlab)
  require(matlabr)
  require(arrow)
  
  
  #####
  
  
  # Check if matlab is installed
  if(!have_matlab()){
    stop("Requires MATLAB!!! Ensure MATLAB is in path.")
  }
  
  # Input checks
  if(!dir.exists(IO_path)){dir.create(IO_path, recursive = TRUE)}
  if(grepl("/$", IO_path)) { IO_path <- sub("/$", "", IO_path) }

  if(!is.logical(csn_weighted)){
    stop("csn_weighted should be logical")
  }
    
  # Create the directory to store results
  if(!dir.exists(paste0(IO_path, "/CSN_results_tmp/"))){dir.create(paste0(IO_path, "/CSN_results_tmp/"), recursive = TRUE)}
  if(!dir.exists(paste0(IO_path, "/CSN_edges_tmp/"))){dir.create(paste0(IO_path, "/CSN_edges_tmp/"), recursive = TRUE)}
  
  
  #####
  
  
  # Store the gene expression as matlab file
  print(paste0(Sys.time(), " |   - Writing input data in matlab format"))
  
  if(!is.matrix(gene_expr)){
    gene_expr <- as.matrix(gene_expr)
  }
  
  writeMat(paste0(IO_path, "/gene_expr.mat"),
           gene_expr = gene_expr,
           gene_ids = row.names(gene_expr),
           sample_ids = colnames(gene_expr))
  
  
  # Compute the edges
  
  print(paste0(Sys.time(), " |   - Computing edges"))
  
  csn_weighted <- ifelse(csn_weighted, 1, 0)
  n_cores_tmp <- min(n_cores, 12) # Max cores in the matlab profile
  matlab_commands <- c( "addpath('Scripts/Functions/CSN/');",
                        paste0("run_csn_parallel(", csn_weighted, ", '", IO_path, "', ", n_cores_tmp, ");") )
  
  run_matlab_code(matlab_commands)
  print(paste0(Sys.time(), " |   --- done!!"))
  gc()
  
  
  #####
  
  
  # Compile the logs
  file_list <- lapply(colnames(gene_expr), function(sample_select){
    paste0(IO_path, "/CSN_results_tmp/log_", sample_select, ".txt")
  })
  
  merged_log <- lapply(file_list, function(x){ readLines(x, warn = FALSE) })
  writeLines(unlist(merged_log), paste0(IO_path, "/merged_log.txt"))
  invisible(file.remove(unlist(file_list)))
  
  
  #####
  
  
  print(paste0(Sys.time(), " |   - Exporting edges"))
  
  source("Scripts/Functions/CSN/Functions_CSN__write_edges.R")
  
  if(use_slurm){
    
    source("Scripts/batch_tools/batchtools_lapply.R")
    
    file_list <- tryCatch(expr = {
      
      batch_lapply( FUN = export_CSN_edges, 
                    FUN_args = list("gene" = row.names(gene_expr)), 
                    FUN_more_args = list("gene_expr" = gene_expr, "IO_path" = IO_path), 
                    IO_path = IO_path, 
                    registry_name = "batchtools_registry__export_edges",
                    packages = c("tidyverse", "arrow", "data.table"), 
                    cluster_template = "Scripts/batch_tools/slurm_batchtools_config.tmpl", 
                    slurm_resources = list(partition = "small",
                                           time = "0-01:00:00",
                                           memory = "25G",
                                           ntasks = 1,
                                           exclude = "kudos12,kudos9,kudos8"), 
                    max_concurrent_jobs = n_cores )
      
    }, 
    warning = function(w){ stop(w) },
    error = function(e){ stop(e) }
    )
    
  }else{
    
    require(parallel)
    
    file_list <- tryCatch(expr = {
      
      mclapply( X = row.names(gene_expr),
                mc.cores = n_cores,
                FUN = export_CSN_edges,
                gene_expr = gene_expr, 
                IO_path = IO_path )
      
    },
    warning = function(w){ stop(w) },
    error = function(e){ stop(e) }
    )
    
  }
  
  print(paste0(Sys.time(), " |   --- done!!"))
  gc()
  
  
  #####
  
  
  print(paste0(Sys.time(), " |   - Exporting final network"))
  network <- open_dataset(unlist(file_list)) 
  
  file.remove(paste0(IO_path, "/gene_expr.mat"))  
  system(paste0("rm -r ", IO_path, "/CSN_results_tmp/"))
  
  print(paste0(Sys.time(), " |   --- done!!"))
  
  return(network)
  
}


