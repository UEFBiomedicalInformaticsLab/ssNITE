set.seed(5081)


# Script compute similarity between pruned SSNs of  TCGA-BRCA samples

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
              help = "Number of jobs in level 1 parallelisation. Default: NULL", metavar = "numeric"),
  make_option(c("--n_jobs"), type = "numeric", default = 1, 
              help = "Number of jobs in level 2 parallelisation. Default: NULL", metavar = "numeric")
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

if(!is.null(opt$n_jobs)){
  if(!is.numeric(opt$n_jobs) | (opt$n_jobs %% 1 != 0)){
    print_help(opt_parser)
    stop("--n_jobs should be be an integer.", call.=FALSE)
  }
}


# Define global options for this script
network_type <- opt$network_type
n_cores <- opt$n_cores
n_jobs <- opt$n_jobs


print(paste0(Sys.time(), " | Using the following parameters: "))
print(paste0(Sys.time(), " | Network type: ", network_type))
print(paste0(Sys.time(), " | n_cores: ", n_cores))
print(paste0(Sys.time(), " | n_jobs: ", n_jobs))


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


# Compute the similarity between the SSNs

print(paste0(Sys.time(), " | Calculating similarity for ", network_type, " SSN"))

source("Scripts/Functions/compare_SSN/Functions_SSN_similarity__main.R")

# Set if the network will be treated as directed
# if(network_type %in% c("ssNITE", "ssNITEpy")){
#   as_directed <- TRUE
# }else{
#   as_directed <- FALSE
# }
as_directed <- FALSE

if(as_directed){
  IO_path <- paste0("IOFiles/SSN_similarity/TCGABRCA/dir_", network_type, "/")
  network_file <- paste0("IOFiles/Pruned_SSN_analysis/TCGABRCA/dir_", network_type, "/", network_type, "_prunedSSN__DA_TCGA__DT_BreastCancer.parquet")
}else{
  IO_path <- paste0("IOFiles/SSN_similarity/TCGABRCA/undir_", network_type, "/")
  network_file <- paste0("IOFiles/Pruned_SSN_analysis/TCGABRCA/undir_", network_type, "/", network_type, "_prunedSSN__DA_TCGA__DT_BreastCancer.parquet")
}

if(!dir.exists(IO_path)){dir.create(IO_path, recursive = TRUE)}


similarity <- func_compute_SSN_similarity(network_file,
                                          sample_info,
                                          as_directed = as_directed,
                                          IO_path = IO_path,
                                          n_cores = n_cores,
                                          n_jobs = n_jobs,
                                          use_slurm = TRUE)

saveRDS(similarity$similarity_data, paste0(IO_path, network_type, "_prunedSSNsimilarity__DA_TCGA__DT_BreastCancer.rds"))

tiff(paste0(IO_path, network_type, "_prunedSSNsimilarity__DA_TCGA__DT_BreastCancer.tiff"),
     width = 25, height = 6,
     units = "cm", compression = "lzw", res = 1200)

print(similarity$similarity_plot)

dev.off()


#####


print(warnings())
