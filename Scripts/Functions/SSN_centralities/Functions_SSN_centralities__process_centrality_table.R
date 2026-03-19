# Function to process the centrality table

process_centrality_table <- function(node, centrality_table_file, IO_path){
  
  require(arrow)
  require(data.table)
  require(tidyverse)
  
  attempt <- 0
  success <- FALSE
  
  while((attempt < 3) & (!success)){
    
    Sys.sleep(runif(1,1, 3) * attempt)
    
    tryCatch(expr = {
      
      centrality_tables <- open_dataset(centrality_table_file)
      tmp1 <- centrality_tables %>% filter(Node1 == node) %>% collect() %>% as.data.table()
      tmp1 <- dcast(tmp1, centrality_type + Node1 + Node2 ~ sample_id, value.var = "centrality_score", fill = 0)
      
      tmp_file <- paste0(IO_path, "tmp_process_table/centrality_scores_", node, ".parquet")
      write_parquet(x = tmp1, sink = tmp_file, compression = "lz4" )
      success <- TRUE
    },  
    error = function(e){
      print(paste0("Attempt ", attempt + 1, " failed for node: ", node))
      print(paste0("Error: ", e$message))
    })
    
    attempt <- attempt + 1
    
  }
  
  if(success){
    message(paste0("Attempt ", attempt, " success for node: ", node))
    return(tmp_file)
  }else{
    stop(paste("Critical Error: Node", node, "failed after 3 attempts."))
  }
  
}



