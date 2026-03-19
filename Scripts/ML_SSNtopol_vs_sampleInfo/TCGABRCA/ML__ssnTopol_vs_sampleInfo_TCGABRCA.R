set.seed(5081)


# Script to train ML model to predict TCGA-BRCA sample info from network level topological properties

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


# Set if the network will be treated as directed and read the network topology data
# if(network_type %in% c("ssNITE", "ssNITEpy")){
#   as_directed <- TRUE
# }else{
#   as_directed <- FALSE
# }
as_directed <- FALSE

if(as_directed){
  IO_path <- paste0("IOFiles/ML_SSNtopol_vs_sampleInfo/TCGABRCA/dir_", network_type, "/")
  network_topol <- read.csv(paste0("IOFiles/Pruned_SSN_analysis/TCGABRCA/dir_", network_type, "/", network_type, "_prunedSSNstats__DA_TCGA__DT_BreastCancer.csv"), check.names = FALSE)
}else{
  IO_path <- paste0("IOFiles/ML_SSNtopol_vs_sampleInfo/TCGABRCA/undir_", network_type, "/")
  network_topol <- read.csv(paste0("IOFiles/Pruned_SSN_analysis/TCGABRCA/undir_", network_type, "/", network_type, "_prunedSSNstats__DA_TCGA__DT_BreastCancer.csv"), check.names = FALSE)
}

if(!dir.exists(IO_path)){dir.create(IO_path, recursive = TRUE)}

network_topol <- network_topol %>% 
  select(c("Sample ID", 
           "Density", "Connected component count", "Largest component size", 
           "Number of isolates", "Transitivity", "Assortativity", 
           "Mean distance", "Diameter", "Modularity", 
           "Number of modules", "Gini coefficient", "Hill index"))


#####


# Read the straified samples

print(paste0(Sys.time(), " | Reading stratified samples"))

stratified_samples <- readRDS("IOFiles/ML_SSNtopol_vs_sampleInfo/TCGABRCA/ML_stratifiedSamples__DA_TCGA_DT_BreastCancer.rds")


#####


# Get the accuracy for classification columns 

print(paste0(Sys.time(), " | Calculating classification accuracy"))
source("Scripts/Functions/machine_learning/Function_train_RF_model.R")

ml_res_1 <- list()

for(sel_col in stratified_samples$classification_cols){
  
  print(paste0(Sys.time(), " | -- running for: ", sel_col))
  
  ml_data <- network_topol %>% 
    left_join(stratified_samples$sample_info %>% 
                select(c("sample_id", all_of(sel_col))) %>% 
                rename("ml_target" = !!sel_col), 
              by = c("Sample ID" = "sample_id")) %>% 
    rename("sample_id" = "Sample ID")
  
  
  ml_res_1[[sel_col]] <- func_train_model( ml_data = ml_data, 
                                           ml_target_type = "classification", 
                                           ml_sample_folds = stratified_samples$stratifiedSamples[[sel_col]], 
                                           n_cores = n_cores )
  
}


ml_res_1 <- rbindlist(ml_res_1, idcol = "ml_target")


#####


# Get the accuracy for regression columns 

print(paste0(Sys.time(), " | Calculating regression accuracy"))
source("Scripts/Functions/machine_learning/Function_train_RF_model.R")

ml_res_2 <- list()

for(sel_col in stratified_samples$regression_cols){
  
  print(paste0(Sys.time(), " | -- running for: ", sel_col))
  
  ml_data <- network_topol %>% 
    left_join(stratified_samples$sample_info %>% 
                select(c("sample_id", all_of(sel_col))) %>% 
                rename("ml_target" = !!sel_col), 
              by = c("Sample ID" = "sample_id")) %>% 
    rename("sample_id" = "Sample ID")
  
  
  ml_res_2[[sel_col]] <- func_train_model( ml_data = ml_data, 
                                           ml_target_type = "regression", 
                                           ml_sample_folds = stratified_samples$stratifiedSamples$Subtype, 
                                           n_cores = n_cores )
  
}


ml_res_2 <- rbindlist(ml_res_2, idcol = "ml_target")


#####


ml_res <- rbindlist(list(ml_res_1, ml_res_2))
saveRDS(ml_res, paste0(IO_path, "MLres_ssnTopol_vs_sampleInfo__DA_TCGA__DT_BreastCancer.rds"))


#####


print(warnings())