set.seed(5081)


# Script to compute pruned SSNs centralities for TCGA-BRCA samples

# Load libraries
library(unixtools)
library(optparse)
library(data.table)
library(arrow)
library(tidyverse)


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


print(paste0(Sys.time(), " | Reading the sample information"))

TCGA_data <- readRDS("IOFiles/GDC_data/vstCounts__DA_TCGA__DT_BreastCancer.rds")
sample_info <- TCGA_data$sample_info %>% filter(condition != "SolidTissueNormal_NA") %>% droplevels()
sample_info <- sample_info %>%
  select(c("barcode", "paper_pathologic_stage", "paper_BRCA_Subtype_PAM50"))  %>% 
  rename(c("sample_id" = "barcode", "Stage" = "paper_pathologic_stage", "Subtype" = "paper_BRCA_Subtype_PAM50"))
row.names(sample_info) <- NULL

rm(TCGA_data)


#####


# Set if the network will be treated as directed
# if(network_type %in% c("ssNITE", "ssNITEpy")){
#   as_directed <- TRUE
# }else{
#   as_directed <- FALSE
# }
as_directed <- FALSE

if(as_directed){
  IO_path <- paste0("IOFiles/SSN_centrality/TCGABRCA/dir_", network_type, "/")
  network_file <- paste0("IOFiles/Pruned_SSN_analysis/TCGABRCA/dir_", network_type, "/", network_type, "_prunedSSN__DA_TCGA__DT_BreastCancer.parquet")
}else{
  IO_path <- paste0("IOFiles/SSN_centrality/TCGABRCA/undir_", network_type, "/")
  network_file <- paste0("IOFiles/Pruned_SSN_analysis/TCGABRCA/undir_", network_type, "/", network_type, "_prunedSSN__DA_TCGA__DT_BreastCancer.parquet")
}

if(!dir.exists(IO_path)){dir.create(IO_path, recursive = TRUE)}


######


# Cluster SSNs based on SSN node/edge centralities

source("Scripts/Functions/SSN_centralities/Functions_SSN_centralities__main.R")

print(paste0(Sys.time(), " | Calculating centrality for ", network_type, " SSN"))

centrality_res <- func_compute_SSN_centralities(network_file,
                                                as_directed = as_directed,
                                                sample_info = sample_info,
                                                ggplot_options = list( "shape" = colnames(sample_info)[2],
                                                                       "color" = colnames(sample_info)[3],
                                                                       "manual_color" = c("Basal" = "#A3A500", "Her2" = "#00BF7D",
                                                                                          "LumA" = "#00B0F6", "LumB" = "#F8766D",
                                                                                          "Normal" = "#E76BF3"),
                                                                       "manual_shape" = c("Stage_I" = 1, "Stage_II" = 2, "Stage_III" = 3, "Stage_IV" = 4) ),
                                                n_cores = n_cores,
                                                IO_path = IO_path,
                                                use_slurm = TRUE)

file.rename(from = centrality_res$centrality_tables, 
            to = paste0(IO_path, "Centrality_", network_type, "_prunedSSN__DA_TCGA__DT_BreastCancer.parquet"))

file.rename(from = centrality_res$pca_res_file, 
            to = paste0(IO_path, "CentralityPCAres_", network_type, "_prunedSSN__DA_TCGA__DT_BreastCancer.rds"))

tiff(paste0(IO_path, "SampleClusterByCentrality_", network_type, "_prunedSSN__DA_TCGA__DT_BreastCancer.tiff"),
     width = 20, height = 8,
     units = "cm", compression = "lzw", res = 1200)

print(centrality_res$sample_cluster_plot)

dev.off()

rm(list = c("centrality_res"))
gc()


#####


print(warnings())