# Function to calculate the SWEET raw edge score


calculate_edge_score <- function(sample_select, gene_expr, ref_net_cor, sample_weights, IO_path, n_jobs = 1){
  
  tryCatch(expr = {
    
    require(arrow)
    require(data.table)
    setDTthreads(n_jobs)
    
    gene_expr_tmp <- cbind(gene_expr, gene_expr[,sample_select, drop = FALSE])
    pert_net_cor <- stats::cor(t(gene_expr_tmp), method = "pearson")
    
    tmp1 <- sample_weights[[sample_select]] * (pert_net_cor - ref_net_cor) + ref_net_cor
    tmp1 <- as.data.table(tmp1, keep.rownames = "regulators")
    tmp1 <- melt(tmp1, id.vars = "regulators", variable.name = "targets", value.name = "raw_edge_score")
    tmp1[, targets := as.character(targets)]
    tmp1 <- tmp1[regulators != targets]
    
    write_parquet(tmp1, sink = paste0(IO_path, "/SWEET_results_tmp/SWEET_results_", sample_select, ".parquet"), compression = "lz4")
    
    return(paste0(IO_path, "/SWEET_results_tmp/SWEET_results_", sample_select, ".parquet"))
    
  },
  warning = function(w){ warning(w) },
  error = function(e){ stop(e) })

}