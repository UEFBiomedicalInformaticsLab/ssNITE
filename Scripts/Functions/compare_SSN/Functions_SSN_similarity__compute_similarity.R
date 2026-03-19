# Function to compute similarity between one vs many samples

compute_similarity <- function(sample_idx, sample_id_list, IO_path, n_jobs = 1){
  
  source("Scripts/Functions/compare_SSN/Functions_SSN_similarity__compare_two_igraphs.R")
  
  sample_comb_grid <- expand.grid(sample_id_list[sample_idx], 
                                  sample_id_list[-c(1:sample_idx)], 
                                  stringsAsFactors = FALSE)
  colnames(sample_comb_grid) <- c("sample1", "sample2")
  sample_comb_grid$sample1_path <- paste0(IO_path, "/tmp_igraph_networks/", sample_comb_grid$sample1, ".rds")
  sample_comb_grid$sample2_path <- paste0(IO_path, "/tmp_igraph_networks/", sample_comb_grid$sample2, ".rds")
  
  res <- tryCatch(expr = {
    
    mclapply(X = 1:nrow(sample_comb_grid), 
             mc.cores = n_jobs, 
             FUN = function(idx){
               
               attempt <- 0
               success <- FALSE
               
               while((attempt < 3) & (!success)){
                 
                 Sys.sleep(runif(1,1, 3) * attempt)
                 
                 tryCatch(expr = {
                   
                   network1 <- sample_comb_grid[idx, ]$sample1_path
                   network2 <- sample_comb_grid[idx, ]$sample2_path
                   
                   if(!file.exists(network1)){ stop(paste0("File missing: ", network1)) }
                   if(!file.exists(network2)){ stop(paste0("File missing: ", network2)) }
                   
                   network1 <- readRDS(network1)
                   network2 <- readRDS(network2)
                   
                   df <- compare_two_igraph_networks(network1, network2)
                   df$Network1 <- sample_comb_grid[idx, ]$sample1
                   df$Network2 <- sample_comb_grid[idx, ]$sample2
                   
                   success <- TRUE
                   
                 }, 
                 error = function(e){
                   message(paste0("Attempt ", attempt + 1, " failed for grid row: ", idx))
                   message(paste0("Error: ", e$message))
                 })
                 
                 attempt <- attempt + 1
                 
               }
               
               if(success){
                 return(df)
               }else{
                 stop(paste("Critical Error: grid row", idx, "failed after 3 attempts."))
               }
               
             })
    
  }, 
  warning = function(w){ warning(w) },
  error = function(e){ stop(e) })
  
  
  res <- data.table::rbindlist(res) %>% select(c("Network1", "Network2", everything()))
  
  file_path <- paste0(IO_path, "/tmp_networks_similarity/similarity_", sample_id_list[sample_idx], ".tsv.gz")
  write.table(res, file = file_path, sep = "\t", row.names = FALSE)
  
  return(file_path)
  
}