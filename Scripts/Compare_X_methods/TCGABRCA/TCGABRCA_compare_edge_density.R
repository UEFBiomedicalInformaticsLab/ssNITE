set.seed(5081)


# Script to calculate the SSNs edge weight statistics for TCGA-BRCA samples across all methods (for publication)


# Load libraries
library(unixtools)
library(optparse)
library(arrow)
library(data.table)
library(tidyverse)
library(parallel)


# Set temporary directory
if(!dir.exists("tmp_dir/")){dir.create("tmp_dir/", recursive = TRUE)}
set.tempdir("tmp_dir/")

if(!dir.exists(paste0("IOFiles/Compare_X_methods/TCGABRCA_compare_edge_density/"))){
  dir.create({paste0("IOFiles/Compare_X_methods/TCGABRCA_compare_edge_density/")})
}


#####


# Define function to plots

func_run_plot <- function(network_type){
  
  require(arrow)
  require(tidyverse)
  
  source("Scripts/Functions/SSN/Functions_SSN_edgeWt_stats.R")
  
  network <- open_dataset(paste0("IOFiles/Create_SSN/TCGABRCA/", network_type, "/", network_type, "_SSN__DA_TCGA__DT_BreastCancer.parquet"))
  
  as_directed <- FALSE
  
  result <- func_plot_edge_weight_stats(network,
                                        as_directed = as_directed,
                                        IO_path = paste0("IOFiles/Compare_X_methods/TCGABRCA_compare_edge_density/", network_type),
                                        n_cores = 10)
  
  unlink(x = paste0("IOFiles/Compare_X_methods/TCGABRCA_compare_edge_density/", network_type), recursive = TRUE)
  
  return(result$edgeWt_density_plot)
  
}


#####


# Make the plots for all network types

source("Scripts/batch_tools/batchtools_lapply.R")

network_type_list <- c("ssNITEpy", "BONOBO", "CSN", "LIONESS", "SWEET")

plot_list <- tryCatch(expr = {
  
  batch_lapply( FUN = func_run_plot,
                FUN_args = list("network_type" = network_type_list),
                IO_path = paste0("IOFiles/Compare_X_methods/TCGABRCA_compare_edge_density/"),
                registry_name = paste0("batchtools_registry__edgeWt_compare"),
                packages = c("igraph", "arrow", "data.table", "tidyverse"),
                cluster_template = "Scripts/batch_tools/slurm_batchtools_config.tmpl",
                slurm_resources = list(partition = "small",
                                       time = "0-12:00:00",
                                       memory = "100G",
                                       ntasks = 10,
                                       exclude = "kudos9"),
                max_concurrent_jobs = 10 )
  
},
error = function(e){ stop(e) }
)

names(plot_list) <- network_type_list
saveRDS(plot_list, "IOFiles/Compare_X_methods/TCGABRCA_compare_edge_density/TCGABRCA_compare_edge_density_data.rds")


#####


# Compile plots

color_list <- c("ssNITEpy" = "#33A02C", "ssNITE" = "#B2DF8A", 
                "BONOBO" = "#1F78B4", "BONOBOsparse" = "#A6CEE3", 
                "CSN" = "#E72F31", "CSNsparse" = "#FB9A99", 
                "SWEET" = "#CC6600", "SWEETsparse" = "#FDBF6F", 
                "LIONESS" = "#6A3D9A")

for(network_type in network_type_list){
  
  tmp1 <- read.csv(paste0("IOFiles/Pruned_SSN_analysis/TCGABRCA/undir_", network_type, "/", network_type, "_prunedSSNstats__DA_TCGA__DT_BreastCancer.csv"),
                   check.names = FALSE)
  
  plot_title <- switch(network_type, "ssNITE" = "ssNITE+", "ssNITEpy" = "ssNITE", network_type)
  
  plot_list[[network_type]] <- plot_list[[network_type]] +  
    geom_line(color = color_list[[network_type]]) +
    # geom_vline(xintercept = unique(tmp1$`Used threshold`), 
    #            linetype = "dotted", 
    #            color = "darkgray", 
    #            linewidth = 0.2) +
    labs( title = plot_title, 
          x = NULL,
          y = NULL ) +
    theme( plot.margin = margin(l = 2, r = 2, t = 3, b = 3),
           text = element_text(size = 6),
           plot.title = element_text(face = "bold", hjust = 0.5, size = 6, margin = margin(b = 0)) )
  
}

plot <- ggpubr::ggarrange(plotlist = plot_list[c("ssNITEpy", "BONOBO", "CSN", "LIONESS", "SWEET")], 
                          ncol = length(plot_list), nrow = 1)

plot <- ggpubr::annotate_figure(plot, 
                                bottom = ggpubr::text_grob("Edge weight", size = 6, face = "bold"),
                                left = ggpubr::text_grob("Density", size = 6, face = "bold", rot = 90))


tiff(paste0("IOFiles/Compare_X_methods/TCGABRCA_compare_edge_density/TCGABRCA_compare_edge_density.tiff"),
     width = 20, height = 5,
     units = "cm", compression = "lzw", res = 1200)

plot

dev.off()


#####


print(warnings())