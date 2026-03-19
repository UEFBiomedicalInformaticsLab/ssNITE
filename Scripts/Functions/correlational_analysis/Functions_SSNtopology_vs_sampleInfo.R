# Function to calculate the correlation between SSN topologies and sample information



func_SSNtopol_vs_sampleInfo_cor <- function(SSN_topol, 
                                            sample_info){
  
  require(tidyverse)
  
  if(!( is.data.frame(SSN_topol) & ("sample_id" %in% colnames(SSN_topol)) & (ncol(SSN_topol) > 1) )){
    stop("'SSN_topol' must be a dataframe with the column 'sample_id' and should have at least one column containing variables with which to correlate")
  }
  
  if(!( is.data.frame(sample_info) & ("sample_id" %in% colnames(sample_info)) & (ncol(sample_info) > 1) )){
    stop("'sample_info' must be a dataframe with the column 'sample_id' and should have at least one column containing variables with which to correlate")
  }
  
  
  # Prepare the sample information
  sample_info <- sample_info %>% column_to_rownames("sample_id")
  
  factor_cols <- names(sample_info)[sapply(sample_info, function(x) is.factor(x) || is.character(x))]
  binary_cols <- names(sample_info)[sapply(sample_info, function(x){all(is.numeric(x)) & (length(na.exclude(unique(x))) == 2)})]
  continuous_cols <- names(sample_info)[sapply(sample_info, function(x) is.numeric(x) & length(na.exclude(unique(x))) > 2)]
  
  # Convert the factor columns to binary
  if(length(factor_cols) > 0){
    sample_info <- WGCNA::binarizeCategoricalColumns(data = sample_info, 
                                                     convertColumns = factor_cols, 
                                                     dropFirstLevelVsAll = FALSE,
                                                     minCount = 10, prefixSep = "::", 
                                                     nameForAll = "")
  }
  
  
  # Prepare the SSN topology data
  SSN_topol <- SSN_topol %>% column_to_rownames("sample_id") 
  
  
  # Remove uncommon samples and reorder
  common_samples <- intersect(row.names(SSN_topol), row.names(sample_info))
  SSN_topol <- SSN_topol[common_samples,]
  sample_info <- sample_info[common_samples,]
  rm(common_samples)
  
  if(!all(row.names(SSN_topol) == row.names(sample_info))){
    SSN_topol <- SSN_topol[row.names(sample_info), ]
  }
  
  
  # Calculate correlations
  cor_matrix <- suppressWarnings( cor(SSN_topol, sample_info, use = "pairwise.complete.obs", method = "spearman"), )
  
  return(cor_matrix)
  
}
