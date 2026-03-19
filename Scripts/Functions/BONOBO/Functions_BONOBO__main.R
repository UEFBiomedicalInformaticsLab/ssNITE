# Function to create single-sample network using BONOBO


create_BONOBO_net <- function( gene_expr, 
                               IO_path,
                               bonobo_sparsify = FALSE,
                               n_cores = 10 ){
  
  require(arrow)
  require(parallel)
  
  
  #####
  
  
  # Input checks
  if(!dir.exists(IO_path)){dir.create(IO_path, recursive = TRUE)}
  if(grepl("/$", IO_path)) { IO_path <- sub("/$", "", IO_path) }
  
  
  # Create the directory to store results
  if(!dir.exists(paste0(IO_path, "/BONOBO_results_tmp/"))){dir.create(paste0(IO_path, "/BONOBO_results_tmp/"), recursive = TRUE)}
  if(!dir.exists(paste0(IO_path, "/BONOBO_edges_tmp/"))){dir.create(paste0(IO_path, "/BONOBO_edges_tmp/"), recursive = TRUE)}
  
  
  #####
  
  
  if(!is.data.frame(gene_expr)){
    stop("gene_expr should be a datatable/dataframe with genes in rows and samples in columns.")
  }
  
  write.table(gene_expr, file = paste0(IO_path, "/gene_expr.tsv"), sep = "\t", row.names = TRUE)
  
  
  #####
  
  
  # Check for python virtual environment
  env_name <- paste0("Environment/python_env_kudos")
  source("Scripts/Functions/BONOBO/Function_BONOBO__load_py.R")
  load_python_env(env_name, create_missing_env = FALSE)
  
  
  # Run BONOBO using reticulate environment
  
  print(paste0(Sys.time(), " |   - Running BONOBO"))
  
  bonobo_obj <- nz$Bonobo(paste0(IO_path, "/gene_expr.tsv"))
  bonobo_obj$run_bonobo(keep_in_memory = FALSE,
                        output_fmt = ".h5",
                        sparsify = bonobo_sparsify,
                        output_folder = paste0(IO_path, "/BONOBO_results_tmp/"),
                        save_pvals = FALSE)
  
  print(paste0(Sys.time(), " |   --- done!!"))
  gc()
  
  
  #####
  
  
  print(paste0(Sys.time(), " |   - Exporting edges"))
  
  Sys.setenv(HDF5_USE_FILE_LOCKING = "FALSE", HDF5_DISABLE_VERSION_CHECK = "2")
  
  file_list <- tryCatch(expr = {
    
    mclapply(X = seq_along(row.names(gene_expr)),
             mc.cores = n_cores,
             FUN = function(gene_idx){
               
               gene <- row.names(gene_expr)[gene_idx]
               
               file_paths <- sapply(colnames(gene_expr), function(sample_select){
                 paste0(IO_path, "/BONOBO_results_tmp/bonobo_", sample_select, ".h5")
               }, USE.NAMES = FALSE)
               
               tmp1 <- lapply(file_paths,
                              function(file){
                                
                                attempt <- 0
                                success <- FALSE
                                
                                while((attempt < 3) & (!success)){
                                  
                                  Sys.sleep(runif(1,1, 3) * attempt) 
                                  
                                  tryCatch(expr = {
                                    
                                    cols <- rhdf5::h5read(file, "/bonobo/axis0")
                                    if(cols[gene_idx] != gene){stop("---")}
                                    df <- rhdf5::h5read(file, "/bonobo/block0_values", index = list(gene_idx, NULL))
                                    colnames(df) <- cols
                                    row.names(df) <- gene
                                    if(bonobo_sparsify){ df[abs(df) > 0] <- 1 } # Binarise the data
                                    df <- as.data.frame(df) %>% rownames_to_column("regulators")
                                    
                                    success <- TRUE
                                    
                                  }, 
                                  error = function(e){
                                    message(paste0("Attempt ", attempt + 1, " failed for file: ", file))
                                    message(paste0("Error: ", e$message))
                                  })
                                  
                                  attempt <- attempt + 1
                                }
                                
                                if(success){
                                  return(df)
                                }else{
                                  stop(paste("Critical Error: ", file, "failed after 3 attempts."))
                                }
                                
                              })
               
               names(tmp1) <- colnames(gene_expr)
               tmp1 <- rbindlist(tmp1, idcol = "sample_id")
               tmp1 <- melt(tmp1, id.vars = c("sample_id", "regulators"), variable.name = "targets", value.name = "score")
               tmp1 <- dcast(tmp1, regulators + targets ~ sample_id, value.var = "score")
               tmp1[, targets := as.character(targets)]
               
               write_parquet(tmp1, sink = paste0(IO_path, "/BONOBO_edges_tmp/BONOBO_edges_", gene, ".parquet"), compression = "lz4")
               rm(tmp1)
               return(paste0(IO_path, "/BONOBO_edges_tmp/BONOBO_edges_", gene, ".parquet"))
               
             })
    
  }, 
  warning = function(w){ warning(w) },
  error = function(e){ stop(e) }
  )  
  
  print(paste0(Sys.time(), " |   --- done!!"))
  gc()
  
  
  #####
  
  
  print(paste0(Sys.time(), " |   - Exporting final network"))
  network <- open_dataset(unlist(file_list)) 
  
  file.remove(paste0(IO_path, "/gene_expr.tsv"))  
  system(paste0("rm -r ", IO_path, "/BONOBO_results_tmp/"))
  
  print(paste0(Sys.time(), " |   --- done!!"))
  
  return(network)
  
}


