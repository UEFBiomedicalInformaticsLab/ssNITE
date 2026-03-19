# Function to create single-sample network using LIONESS

# Using custom implementation based on https://github.com/mararie/lionessR/blob/master/R/lioness.R


create_LIONESS_net <- function( gene_expr, # format : [n_samples, n_features // first column is sample name]
                                IO_path, # If running several instances of the function, ensure different IO_path for each as it generates temporary files within the IO_path
                                n_cores = 10 ){
  
  require(parallel)
  
  #####
  
  
  # Input checks
  if(!dir.exists(IO_path)){dir.create(IO_path, recursive = TRUE)}
  
  if(grepl("/$", IO_path)) { IO_path <- sub("/$", "", IO_path) }
  
  
  if(!is.matrix(gene_expr)){
    stop("gene_expr should be a matrix with samples in columns and features in rows")
  }
  
  
  #####
  
  
  # Create the directory to store results
  if(!dir.exists(paste0(IO_path, "/LIONESS_edges_tmp/"))){dir.create(paste0(IO_path, "/LIONESS_edges_tmp/"), recursive = TRUE)}
  
  
  #####
  
  
  netFun <- function(x){ stats::cor(t(x), method = "pearson") }
  nrsamples <- ncol(gene_expr)
  ref_net_cor <- netFun(gene_expr)
  
  print(paste0(Sys.time(), " |   - Calculating correlations"))
  network <- mclapply(X = colnames(gene_expr),
                      mc.cores = n_cores,
                      FUN = function(sample_select){
                        pert_net_cor <- netFun(gene_expr[, !(colnames(gene_expr) %in% sample_select)])
                        tmp1 <- nrsamples * (ref_net_cor - pert_net_cor) + pert_net_cor
                        return(tmp1)
                      })
  
  names(network) <- colnames(gene_expr)
  
  
  #####
  
  
  print(paste0(Sys.time(), " |   - Exporting edges"))
  file_list <- mclapply(X = row.names(gene_expr), 
                        mc.cores = n_cores, 
                        FUN = function(gene){
                          
                          tmp1 <- lapply(network, function(x){
                            as.data.table(x[gene,, drop = FALSE], keep.rownames = "regulators")
                          })
                          tmp1 <- rbindlist(tmp1, idcol = "sample_id")
                          
                          tmp1 <- melt(tmp1, id.vars = c("sample_id", "regulators"), variable.name = "targets", value.name = "score")
                          tmp1[, score := round(score, 5)]
                          tmp1 <- dcast(tmp1, regulators + targets ~ sample_id, value.var = "score")
                          tmp1[, regulators := as.character(regulators)] 
                          tmp1[, targets := as.character(targets)] 
                          
                          write_parquet(tmp1, sink = paste0(IO_path, "/LIONESS_edges_tmp/LIONESS_edges_", gene, ".parquet"), compression = "lz4")
                          return(paste0(IO_path, "/LIONESS_edges_tmp/LIONESS_edges_", gene, ".parquet"))
                          
                        })
  
  
  print(paste0(Sys.time(), " |   - Exporting final network"))
  network <- open_dataset(unlist(file_list)) 
  
  return(network)
  
}


