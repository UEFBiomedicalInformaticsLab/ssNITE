set.seed(5081)


# Script to calculate the SSNs edge weight statistics for TCGA-BRCA samples


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


#####


# Get arguments
option_list = list(
  make_option(c("--network_type"), type = "character", default = NULL,
              help = "Network type to be used. Possible values: ssNITE, ssNITEpy, BONOBO, CSN, LIONESS, SWEET, BONOBOsparse, CSNsparse, SWEETsparse", metavar = "character"),
  make_option(c("--n_cores"), type = "numeric", default = 1,
              help = "Number of jobs in level 1 parallelisation. Default: NULL", metavar = "numeric")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)


if(is.null(opt$network_type)){
  print_help(opt_parser)
  stop("--network_type argument needed", call. = FALSE)
}

if(!opt$network_type %in% c("ssNITE", "ssNITEpy", "BONOBO", "CSN", "LIONESS", "SWEET", "BONOBOsparse", "CSNsparse", "SWEETsparse")){
  print_help(opt_parser)
  stop("--network_type must be one of the following: ssNITE, ssNITEpy, BONOBO, CSN, LIONESS, SWEET, BONOBOsparse, CSNsparse, SWEETsparse", call. = FALSE)
}

if(!is.null(opt$n_cores)){
  if(!is.numeric(opt$n_cores) | (opt$n_cores %% 1 != 0)){
    print_help(opt_parser)
    stop("--n_cores should be be an integer.", call.=FALSE)
  }
}


# Define global options for this script
network_type <- opt$network_type
n_cores <- opt$n_cores


print(paste0(Sys.time(), " | Using the following parameters: "))
print(paste0(Sys.time(), " | Network type: ", network_type))
print(paste0(Sys.time(), " | n_cores: ", n_cores))


#####


print(paste0(Sys.time(), " | Plotting SSN edge weight stats"))

source("Scripts/Functions/SSN/Functions_SSN_edgeWt_stats.R")

network <- open_dataset(paste0("IOFiles/Create_SSN/TCGABRCA/", network_type, "/", network_type, "_SSN__DA_TCGA__DT_BreastCancer.parquet"))

# if(network_type %in% c("ssNITE", "ssNITEpy")){
#   as_directed <- TRUE
# }else{
#   as_directed <- FALSE
# }
as_directed <- FALSE

result <- func_plot_edge_weight_stats(network, 
                                      as_directed = as_directed,
                                      IO_path = paste0("IOFiles/Create_SSN/TCGABRCA/", network_type), 
                                      n_cores = n_cores)

write.csv(result$edgeWt_stats, paste0("IOFiles/Create_SSN/TCGABRCA/", network_type, "/", network_type, "_edgeWt_stats__DA_TCGA__DT_BreastCancer.csv"), row.names = FALSE)

tiff(paste0("IOFiles/Create_SSN/TCGABRCA/", network_type, "/", network_type, "_edgeWt_stats__DA_TCGA__DT_BreastCancer.tiff"),
     width = 20, height = 6,
     units = "cm", compression = "lzw", res = 1200)

print(result$edgeWt_stats_plot)

dev.off()


tiff(paste0("IOFiles/Create_SSN/TCGABRCA/", network_type, "/", network_type, "_edgeWt_density__DA_TCGA__DT_BreastCancer.tiff"),
     width = 6, height = 6,
     units = "cm", compression = "lzw", res = 1200)

print(result$edgeWt_density_plot)

dev.off()

rm(list = c("network", "result"))
gc()


#####


print(warnings())