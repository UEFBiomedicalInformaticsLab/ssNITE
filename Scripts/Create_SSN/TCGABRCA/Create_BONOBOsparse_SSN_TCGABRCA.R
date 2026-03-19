set.seed(5081)


# Script to prepare SSN using BONOBO (Bayesian Optimized sample-specific Networks Obtained By Omics data) for TCGA-BRCA samples [sparse mode]

# Load libraries
library(unixtools)
library(reticulate)
library(arrow)
library(data.table)
library(tidyverse)
library(parallel)


# Set temporary directory
if(!dir.exists("tmp_dir/")){dir.create("tmp_dir/", recursive = TRUE)}
set.tempdir("tmp_dir/")

if(!dir.exists("IOFiles/Create_SSN/TCGABRCA/BONOBOsparse/")){dir.create("IOFiles/Create_SSN/TCGABRCA/BONOBOsparse/", recursive = TRUE)}


#####


# Read the processed data

print(paste0(Sys.time(), " | Reading gene expression data"))

TCGA_data <- readRDS("IOFiles/GDC_data/vstCounts__DA_TCGA__DT_BreastCancer.rds")

rnaSeq_vst <- TCGA_data$rnaSeq_vst
sample_info <- TCGA_data$sample_info %>% filter(condition != "SolidTissueNormal_NA") %>% droplevels()
rm(TCGA_data)


# #####
# 
# 
# # Randomly select samples from the whole dataset
# selected_samples <- c()
# for(condition_select in as.vector(unique(sample_info$condition))){
#   tmp1 <- sample_info %>% filter(condition == condition_select) %>% pull(barcode)
#   tmp1 <- tmp1 %>% sample(ceiling(length(tmp1) * 0.1))
#   selected_samples <- c(selected_samples, tmp1)
# }
# sample_info <- sample_info %>% filter(barcode %in% selected_samples)
# rm(list = c("selected_samples", "tmp1"))
# 
# 
# #####


# Read the background gene list

print(paste0(Sys.time(), " | Preparing background gene list"))

gene_list <- read.csv("IOFiles/SSN_background_genes/SSN_background_genes.csv", header = TRUE)
gene_list <- unique(gene_list$ensembl_gene_id)
print(paste0(Sys.time(), " | Number of genes in background gene list: ", length(gene_list)))

gene_list <- gene_list[gene_list %in% row.names(rnaSeq_vst)]
print(paste0(Sys.time(), " |    Using genes: ", length(gene_list)))


#####


# Create BONOBO networks
print(paste0(Sys.time(), " | Creating networks for gene list"))

rnaSeq_vst <- rnaSeq_vst[row.names(rnaSeq_vst) %in% gene_list, colnames(rnaSeq_vst) %in% sample_info$barcode]

source("Scripts/Functions/BONOBO/Functions_BONOBO__main.R")

network <- create_BONOBO_net( gene_expr = rnaSeq_vst,
                              IO_path = "IOFiles/Create_SSN/TCGABRCA/BONOBOsparse/", 
                              bonobo_sparsify = TRUE, 
                              n_cores = 10 )

write_parquet( x = network,
               sink = "IOFiles/Create_SSN/TCGABRCA/BONOBOsparse/BONOBOsparse_SSN__DA_TCGA__DT_BreastCancer.parquet",
               compression = "lz4" )

system("rm -r IOFiles/Create_SSN/TCGABRCA/BONOBOsparse/BONOBO_edges_tmp/")
gc()


#####


print(warnings())