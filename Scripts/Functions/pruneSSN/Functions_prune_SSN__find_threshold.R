# Function to select the threshold to trim edges from the full network

find_threshold <- function(sample_select, network_file, number_thresholds = 100, n_jobs = 1){
  
  require(arrow)
  require(data.table)
  require(tidyverse)
  require(parallel)
  
  # Subset the network for single sample
  network <- open_dataset(network_file)
  network <- network %>%
    select(c("regulators", "targets", all_of(sample_select))) %>% 
    rename( "weight" = !!sample_select) %>%
    collect() %>% as.data.table()
  
  # Remove loops
  network <- network[network$regulators != network$targets,]
  
  # Get the thresholds to screen
  edge_wts <- network %>% pull(weight) %>% .[. != 0]
  screen_thresholds <- quantile(edge_wts, 
                                probs = seq(0.001, 0.999, length.out = number_thresholds - 2), 
                                na.rm = TRUE, 
                                names = FALSE)
  screen_thresholds <- c(0, min(edge_wts), screen_thresholds, max(edge_wts))
  
  # Prepare the adjacency
  network <- dcast(network, regulators ~ targets, value.var = "weight")
  network <- network %>%
    column_to_rownames("regulators") %>% 
    as.matrix()
  network <- network[sort(row.names(network)), sort(colnames(network))]
  
  # Fit linear model on the log-log node degree distribution
  scale_free_fit <- mclapply(X = screen_thresholds,
                             mc.cores = n_jobs,
                             FUN = function(cut){
                               
                               undirected_mat <- (network > cut) | t(network > cut)
                               node_degree <- rowSums(undirected_mat, na.rm = TRUE)
                               node_degree_filt <- node_degree[node_degree > 0]
                               
                               if(length(unique(node_degree_filt)) >= 5){
                                 
                                 degree_val <- sort(unique(node_degree_filt))
                                 ccdf <- sapply(degree_val, function(k) sum(node_degree_filt >= k)/length(node_degree_filt))
                                 
                                 df <- data.frame( "node_degree" = degree_val,
                                                   "ccdf" = ccdf )
                                 
                                 fit <- lm(log10(ccdf) ~ log10(node_degree), data = df)
                                 fit <- data.frame( "cut" = cut,
                                                    "Rsquared" = summary(fit)$r.squared,
                                                    "slope" = summary(fit)$coefficients[2,1],
                                                    "isolates" = length(node_degree[node_degree == 0]),
                                                    "max_degree" = max(node_degree) )
                                 
                               } else {
                                 fit <- data.frame( "cut" = cut,
                                                    "Rsquared" = NA, 
                                                    "slope" = NA, 
                                                    "isolates" = length(node_degree[node_degree == 0]),
                                                    "max_degree" = max(node_degree) )
                               }
                               
                               return(fit)
                               
                             })
  
  scale_free_fit <- rbindlist(scale_free_fit, fill = TRUE)
  scale_free_fit$approx_alpha <- 1 - scale_free_fit$slope
  
  scale_free_fit <- scale_free_fit %>% 
    filter(Rsquared >= 0.8) %>% 
    filter(isolates <= round(nrow(network) * 0.25)) %>% 
    filter((approx_alpha > 1.5) & (approx_alpha < 3.5))
  
  if(nrow(scale_free_fit) > 0){
    scale_free_fit <- scale_free_fit %>% arrange(cut) %>% slice(1) %>% pull(cut)
  }else{
    scale_free_fit <- NA_real_
  }
  
  rm(list = c("network"))
  
  return(scale_free_fit)
  
}