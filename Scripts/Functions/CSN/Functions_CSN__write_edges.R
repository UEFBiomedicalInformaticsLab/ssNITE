# Function to compile the edges in the desired format

export_CSN_edges <- function(gene, gene_expr, IO_path){
  
  require(arrow)
  require(data.table)
  require(tidyverse)
  
  attempt <- 0
  success <- FALSE
  
  while((attempt < 3) & (!success)){
    
    Sys.sleep(runif(1,1, 3) * attempt)
    
    tryCatch(expr = {
      
      csn_result <- sapply(colnames(gene_expr),
                           function(sample_select){paste0(IO_path, "/CSN_results_tmp/CSN_result_", sample_select, ".parquet")},
                           USE.NAMES = FALSE)
      csn_result <- open_dataset(csn_result)
      
      csn_result <- csn_result %>% filter(Row %in% gene) %>% collect() %>% as.data.table()
      csn_result <- setnames(csn_result, "Row", "regulators")
      csn_result <- melt(csn_result, id.vars = c("regulators", "sample_id"), variable.name = "targets", value.name = "score")
      csn_result <- dcast(csn_result, regulators + targets ~ sample_id, value.var = "score")
      csn_result[, targets := as.character(targets)]
      
      write_parquet(csn_result, sink = paste0(IO_path, "/CSN_edges_tmp/CSN_edges_", gene, ".parquet"), compression = "lz4")
      
      success <- TRUE
      
    }, 
    error = function(e){
      message(paste0("Attempt ", attempt + 1, " failed for gene: ", gene))
      message(paste0("Error: ", e$message))
    })
    
    attempt <- attempt + 1
    
  }
  
  if(exists("csn_result")){rm(csn_result)}
  gc()
  
  if(success){
    message(paste0("Attempt ", attempt, " success for gene: ", gene))
    return(paste0(IO_path, "/CSN_edges_tmp/CSN_edges_", gene, ".parquet"))
  }else{
    stop(paste("Critical Error: Gene", gene, "failed after 3 attempts."))
  }
  
}