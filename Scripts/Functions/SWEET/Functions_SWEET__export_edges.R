# Function to export the edges in the desired format
# If the mean and SD values are provided, calculates the z-scores, and binarises the network


export_SWEET_edges <- function(gene, sample_id_list, mean_edgeWt, sd_edgeWt, IO_path, n_jobs = 1){
  
  require(arrow)
  require(data.table)
  require(tidyverse)
  setDTthreads(n_jobs)
  
  attempt <- 0
  success <- FALSE
  
  while((attempt < 3) & (!success)){
    
    Sys.sleep(runif(1,1, 3) * attempt)
    
    tryCatch(expr = {
      
      file_list <- paste0(IO_path, "/SWEET_results_tmp/SWEET_results_", sample_id_list, ".parquet")
      edges <- open_dataset(file_list) %>% mutate(sample_id = add_filename())
      edges <- edges %>% filter(regulators %in% gene) %>% collect() %>% as.data.table()
      edges[, sample_id := sub(".*SWEET_results_tmp/SWEET_results_(.*)\\.parquet", "\\1", sample_id)] 
      
      if( !(is.na(mean_edgeWt)) &  !(is.na(sd_edgeWt)) ){
        edges[, zscore := (raw_edge_score - mean_edgeWt) / sd_edgeWt]
        edges$raw_edge_score <- ifelse(abs(edges$zscore) > 1.960, 1, 0)
        edges[, zscore := NULL]
      }
      
      edges <- dcast(edges, regulators + targets ~ sample_id, value.var = "raw_edge_score")
      write_parquet(edges, sink = paste0(IO_path, "/SWEET_edges_tmp/SWEET_edges_", gene, ".parquet"), compression = "lz4")
      
      success <- TRUE
      
    },  
    error = function(e){
      message(paste0("Attempt ", attempt + 1, " failed for gene: ", gene))
      message(paste0("Error: ", e$message))
    })
    
    attempt <- attempt + 1
    
  }
  
  if(exists("edges")){rm(edges)}
  gc()
  
  
  if(success){
    message(paste0("Attempt ", attempt, " success for gene: ", gene))
    return(paste0(IO_path, "/SWEET_edges_tmp/SWEET_edges_", gene, ".parquet"))
  }else{
    stop(paste("Critical Error: Gene", gene, "failed after 3 attempts."))
  }
  
}

