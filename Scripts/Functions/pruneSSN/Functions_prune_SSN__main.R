# Function to prune SSN

prune_SSN <- function( network_file, 
                       as_directed = TRUE, 
                       n_cores = 1,
                       n_jobs = 1,
                       number_thresholds = 100, 
                       use_common_threshold = TRUE,
                       IO_path = NULL,
                       use_slurm = FALSE ){
  
  require(arrow)
  require(tidyverse)
  
  # Check if network file exists
  if(!file.exists(network_file)){
    stop(paste0("File does not exist: ", network_file))
  }
  
  # Check if the input network dataframe is correct
  network <- open_dataset(network_file)
  if (!all(c("regulators", "targets") %in% names(network))) {
    stop("Input network must contain columns: 'regulators' and 'targets'")
  }
  
  # Get list of samples
  sample_id_list <- network %>% select(!c("regulators", "targets")) %>% names()
  
  if (length(sample_id_list) == 0) {
    stop("Input network must contain at least one sample column besides 'regulators' and 'targets'")
  }
  
  
  #####
  
  
  # Calculate the threshold based on scale-free topology
  print(paste0(Sys.time(), " |  - Calculating theresholds"))
  
  source("Scripts/Functions/pruneSSN/Functions_prune_SSN__find_threshold.R")
  
  if(use_slurm){
    
    source("Scripts/batch_tools/batchtools_lapply.R")
    
    select_thresholds <- tryCatch(expr = {
      
      batch_lapply( FUN = find_threshold, 
                    FUN_args = list("sample_select" = sample_id_list), 
                    FUN_more_args = list("network_file" = network_file, 
                                         "number_thresholds" = number_thresholds), 
                    IO_path = IO_path, 
                    registry_name = "batchtools_registry__prune_SSN",
                    packages = c("igraph", "arrow", "data.table", "tidyverse"), 
                    cluster_template = "Scripts/batch_tools/slurm_batchtools_config.tmpl", 
                    slurm_resources = list(partition = "small",
                                           time = "0-01:00:00",
                                           memory = "10G",
                                           ntasks = n_jobs,
                                           exclude = "kudos9,kudos8"), 
                    max_concurrent_jobs = n_cores )
      
    }, 
    warning = function(w){ stop(w) },
    error = function(e){ stop(e) }
    )
    
  }else{
    
    require(parallel)
    
    select_thresholds <- tryCatch(expr = {
      
      mclapply(X = sample_id_list,
               mc.cores = n_cores,
               FUN = find_threshold,
               network_file = network_file, 
               number_thresholds = number_thresholds, 
               n_jobs = n_jobs)
      
    }, 
    warning = function(w){ stop(w) },
    error = function(e){ stop(e) }
    )
    
  }
  
  names(select_thresholds) <- sample_id_list
  select_thresholds <- unlist(select_thresholds)
  
  
  # Check for failed threshold selection
  select_failed <- select_thresholds[is.na(select_thresholds)]
  if(length(select_failed) > 0){
    warning(paste0( "Scale-free based threshold selection criteria not reached for ", 
                    length(select_failed), "/" , length(select_thresholds), " (",
                    (length(select_failed)/length(select_thresholds)) * 100, "%) samples" ), 
            immediate. = TRUE)
  }
  
  
  # Set the threshold to use
  use_thresholds <- select_thresholds
  
  if(length(select_failed) == length(select_thresholds)){
    
    # What happens if all are failed?
    # Take the top 5%tile for each sample
    # If use_common_threshold is selected, the median of all the 5%tile values will be used
    
    alternate_thresholds <- sapply(sample_id_list, function(x){
      y <- network %>% pull(x, as_vector = TRUE)
      quantile(y, probs = 0.95, na.rm = TRUE, names = FALSE)
    })
    
    if(use_common_threshold){
      tmp1 <- median(alternate_thresholds, na.rm = TRUE)
      use_thresholds[] <- tmp1
      rm(tmp1)
    }else{
      use_thresholds <- alternate_thresholds
    }
    rm(alternate_thresholds)
    
  }else{
    
    # If threshold selection is not failure for all, 
    # the median value is used if use_common_threshold is selected
    # else the 5% value is used for the failed ones
    
    if(use_common_threshold){
      tmp1 <- median(select_thresholds, na.rm = TRUE)
      use_thresholds[] <- tmp1
      rm(tmp1)
    }else{
      
      if(length(select_failed) > 0){
        alternate_thresholds <- sapply(names(select_failed), function(x){
          y <- network %>% pull(x, as_vector = TRUE)
          quantile(y, probs = 0.95, na.rm = TRUE, names = FALSE)
        })
        use_thresholds <- select_thresholds[!is.na(select_thresholds)]
        use_thresholds <- c(use_thresholds, alternate_thresholds)
        use_thresholds <- use_thresholds[sample_id_list]
      }
      
    }
    
  }
  
  
  print(paste0(Sys.time(), " |  --- done!!"))
  
  gc()
  
  
  #####
  
  
  # Generate the pruned network
  # First creates a named list of expression/function to be applied
  # on each column and then used the mutate to execute
  
  print(paste0(Sys.time(), " |  - Pruning the networks"))
  
  expr_list <- list()
  for(column in sample_id_list){
    threshold <- use_thresholds[[column]]
    expr_list[[column]] <- expr(
      if_else(!!sym(column) > !!threshold, !!sym(column), NA)
    )
  }
  
  network <- network %>% mutate(!!!expr_list)
  print(paste0(Sys.time(), " |  --- done!!"))
  
  
  #####
  
  
  # Generate statistics about the pruned network
  print(paste0(Sys.time(), " |  - Calculating pruned network statistics"))
  
  source("Scripts/Functions/pruneSSN/Functions_prune_SSN__network_stats.R")
  
  if(use_slurm){
    
    source("Scripts/batch_tools/batchtools_lapply.R")
    
    pruned_net_stats <- tryCatch(expr = {
      
      batch_lapply( FUN = get_network_stats, 
                    FUN_args = list("sample_select" = sample_id_list), 
                    FUN_more_args = list("network_file" = network_file, 
                                         "prune_thresholds" = use_thresholds,
                                         "as_directed" = as_directed), 
                    IO_path = IO_path, 
                    registry_name = "batchtools_registry__pruned_SSN_stats",
                    packages = c("igraph", "arrow", "data.table", "tidyverse"), 
                    cluster_template = "Scripts/batch_tools/slurm_batchtools_config.tmpl", 
                    slurm_resources = list(partition = "small",
                                           time = "0-01:00:00",
                                           memory = "10G",
                                           ntasks = 1,
                                           exclude = "kudos9,kudos8"), 
                    max_concurrent_jobs = n_cores )
      
    }, 
    warning = function(w){ stop(w) },
    error = function(e){ stop(e) })
    
  }else{
    
    require(parallel)
    
    pruned_net_stats <- tryCatch(expr = {
      
      mclapply(X = sample_id_list, 
               mc.cores = n_cores, 
               FUN = get_network_stats, 
               network_file = network_file,
               prune_thresholds = use_thresholds,
               as_directed = as_directed)
      
    }, 
    warning = function(w){ stop(w) },
    error = function(e){ stop(e) })
    
  }
  
  pruned_net_stats <- pruned_net_stats %>% 
    bind_rows() %>% 
    left_join(select_thresholds %>% 
                as_tibble(rownames = "Sample ID") %>% rename("Selected threshold" = "value"),
              by = "Sample ID") %>% 
    left_join(use_thresholds %>% 
                as_tibble(rownames = "Sample ID") %>% rename("Used threshold" = "value"),
              by = "Sample ID") %>% 
    as.data.table()
  
  print(paste0(Sys.time(), " |  --- done!!"))
  
  gc()
  
  
  #####
  
  
  # Generate plot using the stats
  print(paste0(Sys.time(), " |  - Plotting the pruned network statistics"))
  
  plot <- melt(pruned_net_stats, 
               id.vars = "Sample ID", 
               variable.name = "Metric", 
               value.name = "Score")
  
  plot <- ggplot(plot, aes(x = Score)) +
    geom_freqpoly(linewidth = 0.1, na.rm = TRUE) +
    facet_wrap(~ Metric, scale = "free") +
    labs(x = "Value",
         y = "Number of samples") +
    theme( panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.25, linetype = NULL),
           panel.grid = element_blank(),
           panel.spacing = unit(0.1, "cm"),
           strip.background = element_rect(color = "black", linewidth = 0.25,),
           strip.text = element_text(size = 5, margin = margin(t = 1, b = 1, r = 1, l= 1)),
           text = element_text(size = 8),
           plot.title = element_text(face = "bold", hjust = 0.5, size = 8),
           axis.text.x = element_text(size = 4, angle = 0, vjust = 0, hjust = 0.5),
           axis.text.y = element_text(size = 4),
           axis.ticks = element_line(colour = "black", linewidth = 0.2),
           legend.position = "bottom",
           legend.key = element_rect(fill = NA),
           legend.key.size = unit(0.4, "cm"),
           legend.key.spacing.y = unit(0.1, "cm"),
           legend.title = element_text(size = 6, face = "bold", margin = margin(t = 1, r = 0, b = 1.5, l = 0)),
           legend.text = element_text(size = 6, margin = margin(t = 0.5, r = 0, b = 0.1, l = 1)),
           legend.margin = margin(1,1,1,1),
           legend.spacing = unit(0.2, "cm") )
  
  
  print(paste0(Sys.time(), " |  --- done!!"))
  
  gc()
  
  
  #####
  
  
  return(list("pruned_network" = network, 
              "pruned_network_stats" = pruned_net_stats, 
              "pruned_network_stat_plot" = plot))
  
}