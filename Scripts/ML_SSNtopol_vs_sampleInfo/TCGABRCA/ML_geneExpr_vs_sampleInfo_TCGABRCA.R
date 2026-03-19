set.seed(5081)


# Script to train ML model to predict TCGA-BRCA sample info from gene expression

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


IO_path <- paste0("IOFiles/ML_SSNtopol_vs_sampleInfo/TCGABRCA/")


#####


# Read the processed data

print(paste0(Sys.time(), " | Reading gene expression data"))

TCGA_data <- readRDS("IOFiles/GDC_data/vstCounts__DA_TCGA__DT_BreastCancer.rds")
rnaSeq_vst <- TCGA_data$rnaSeq_vst
rm(TCGA_data)


# Read the background gene list

print(paste0(Sys.time(), " | Preparing background gene list"))

gene_list <- read.csv("IOFiles/SSN_background_genes/SSN_background_genes.csv", header = TRUE)
gene_list <- unique(gene_list$ensembl_gene_id)
print(paste0(Sys.time(), " | Number of genes in background gene list: ", length(gene_list)))

gene_list <- gene_list[gene_list %in% row.names(rnaSeq_vst)]
print(paste0(Sys.time(), " |    Using genes: ", length(gene_list)))


#####


# Read the stratified samples

print(paste0(Sys.time(), " | Reading stratified samples"))

stratified_samples <- readRDS("IOFiles/ML_SSNtopol_vs_sampleInfo/TCGABRCA/ML_stratifiedSamples__DA_TCGA_DT_BreastCancer.rds")


rnaSeq_vst <- rnaSeq_vst[row.names(rnaSeq_vst) %in% gene_list, colnames(rnaSeq_vst) %in% stratified_samples$sample_info$sample_id]
rnaSeq_vst <- rnaSeq_vst %>% t()
pca_res <- stats::prcomp(x = rnaSeq_vst, center = TRUE, scale. = FALSE)
pca_res <- pca_res$x[,1:10]


#####


# Get the accuracy for classification columns 

print(paste0(Sys.time(), " | Calculating classification accuracy"))
source("Scripts/Functions/machine_learning/Function_train_RF_model.R")

ml_res_1 <- list()

for(sel_col in stratified_samples$classification_cols){
  
  print(paste0(Sys.time(), " | -- running for: ", sel_col))
  
  ml_data <- pca_res %>% as.data.frame() %>% rownames_to_column("sample_id") %>% 
    left_join(stratified_samples$sample_info %>% 
                select(c("sample_id", all_of(sel_col))) %>% 
                rename("ml_target" = !!sel_col), 
              by = "sample_id")
  
  
  ml_res_1[[sel_col]] <- func_train_model( ml_data = ml_data, 
                                           ml_target_type = "classification", 
                                           ml_sample_folds = stratified_samples$stratifiedSamples[[sel_col]], 
                                           n_cores = 10 )
  
}

ml_res_1 <- rbindlist(ml_res_1, idcol = "ml_target")


#####


# Get the accuracy for regression columns 

print(paste0(Sys.time(), " | Calculating regression accuracy"))
source("Scripts/Functions/machine_learning/Function_train_RF_model.R")

ml_res_2 <- list()

for(sel_col in stratified_samples$regression_cols){
  
  print(paste0(Sys.time(), " | -- running for: ", sel_col))
  
  ml_data <- pca_res %>% as.data.frame() %>% rownames_to_column("sample_id") %>% 
    left_join(stratified_samples$sample_info %>% 
                select(c("sample_id", all_of(sel_col))) %>% 
                rename("ml_target" = !!sel_col), 
              by = "sample_id")
  
  
  ml_res_2[[sel_col]] <- func_train_model( ml_data = ml_data, 
                                           ml_target_type = "regression", 
                                           ml_sample_folds = stratified_samples$stratifiedSamples$Subtype, 
                                           n_cores = 10 )
  
}

ml_res_2 <- rbindlist(ml_res_2, idcol = "ml_target")


#####


ml_res <- rbindlist(list(ml_res_1, ml_res_2))
saveRDS(ml_res, paste0(IO_path, "MLres_geneExpr_vs_sampleInfo__DA_TCGA__DT_BreastCancer.rds"))


#####


print(warnings())