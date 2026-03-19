set.seed(5081)


# Script to check the clustering of the samples using expression of SSN background genes in TCGA-BRCA processed RNA-Seq data (only tumour samples)


# Load libraries
library(unixtools)
library(UpSetR)
library(tidyverse)


# Set temporary directory
if(!dir.exists("tmp_dir/")){dir.create("tmp_dir/", recursive = TRUE)}
set.tempdir("tmp_dir/")


#####


# Read the processed data

print("Reading gene expression data")

TCGA_data <- readRDS("IOFiles/GDC_data/vstCounts__DA_TCGA__DT_BreastCancer.rds")

rnaSeq_vst <- TCGA_data$rnaSeq_vst
sample_info <- TCGA_data$sample_info %>% filter(condition != "SolidTissueNormal_NA") %>% droplevels()
rm(TCGA_data)


#####


# Read the background gene list

print("Preparing background gene list")

gene_list <- read.csv("IOFiles/SSN_background_genes/SSN_background_genes.csv", header = TRUE)

print(paste0("Number of genes in background gene list: ", length( unique(gene_list$ensembl_gene_id) )))

gene_list <- gene_list[gene_list$ensembl_gene_id %in% row.names(rnaSeq_vst), ]
print(paste0("    Using genes: ", length( unique(gene_list$ensembl_gene_id) )))


# Print the overlap between the gene sets after filtering
plot_data <- gene_list %>% select(c(ensembl_gene_id, source)) %>% separate_rows(source, sep = ";") 
plot_data <- split(plot_data$ensembl_gene_id, plot_data$source)
plot_data <- lapply(plot_data, unique)


tiff(paste0("IOFiles/SSN_background_genes/SSN_backgroundGenes_overlap__DA_TCGA__DT_BreastCancer.tiff"),
     width = 35, 
     height = 21,
     units = "cm", compression = "lzw", res = 1200)

upset( data = fromList(plot_data), 
       nsets = length(plot_data), 
       nintersects = NA, 
       keep.order = FALSE, 
       order.by = "freq", 
       text.scale = c(1, 1, 1, 1, 1, 0.75), 
       number.angles = 0 )

dev.off()


tiff(paste0("IOFiles/SSN_background_genes/SSN_backgroundGenes_top25overlap__DA_TCGA__DT_BreastCancer.tiff"),
     width = 20,
     height = 12,
     units = "cm", compression = "lzw", res = 1200)

upset( data = fromList(plot_data),
       nsets = length(plot_data),
       nintersects = 25,
       keep.order = FALSE,
       order.by = "freq",
       text.scale = c(1, 1, 0.9, 0.9, 0.8, 0.8), 
       number.angles = 0 )

dev.off()


gene_list <- unique(gene_list$ensembl_gene_id)


#####


# Plot PCA with all samples
rnaSeq_vst_filt <- rnaSeq_vst[row.names(rnaSeq_vst) %in% gene_list, colnames(rnaSeq_vst) %in% sample_info$barcode]

pca_res <- prcomp(t(rnaSeq_vst_filt))
pca_eigen <- factoextra::get_eig(pca_res)

pca_data <- as_tibble(pca_res$x, rownames = "sample_id")

pca_data <- pca_data %>%
  left_join(sample_info %>% dplyr::select(c("barcode", "condition", "paper_pathologic_stage")),
            by = c("sample_id" = "barcode"))  %>%
  separate(condition,
           into = c("sample_type", "cancer_subtype"),
           sep = "_") %>%
  rename("stage" = "paper_pathologic_stage")
pca_data$stage[pca_data$stage == "NA"] <- NA

plot_1 <- ggplot() +
  geom_point(data = pca_data, 
             mapping = aes(x = PC1, y = PC2, shape = stage, color = cancer_subtype),
             size = 2, 
             stroke = 0.5) +
  labs(title = "all samples", 
       x = paste0("PC1: ", round(pca_eigen$variance.percent[1], 2), "% variance"),
       y = paste0("PC2: ", round(pca_eigen$variance.percent[2], 2) , "% variance"), 
       shape = "Stage: ",
       color = "Subtype: ") +
  scale_color_manual(values = c("Basal" = "#A3A500", "Her2" = "#00BF7D", 
                                "LumA" = "#00B0F6", "LumB" = "#F8766D", 
                                "Normal" = "#E76BF3"), na.value = "#A9A9A9") + 
  scale_shape_manual(values = c("Stage_I" = 0, "Stage_II" = 1, "Stage_III" = 2, "Stage_IV" = 3), na.value = 4) + 
  theme(panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.25, linetype = NULL),
        panel.grid = element_blank(),
        panel.spacing = unit(0.1, "cm"),
        text = element_text(size = 8),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 8),
        axis.text.x = element_text(angle = 0, vjust = 0, hjust = 0.5),
        axis.ticks = element_line(colour = "black", linewidth = 0.2),
        legend.position = "bottom",
        legend.key = element_rect(fill = NA), 
        legend.key.size = unit(0.5, "cm"),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.margin = margin(1,1,1,1),
        legend.spacing = unit(0.2, "cm")  ) 

tiff(paste0("IOFiles/SSN_background_genes/plot_PCA_SSN_bkgGene__DA_TCGA__DT_BreastCancer.tiff"),
     width = 18, height = 15,
     units = "cm", compression = "lzw", res = 1200)

print(plot_1)

dev.off()


#####


# Read TNBC samples list

TNBC_samples <- readRDS("IOFiles/GDC_data/sample_groups__DA_TCGA__DT_BreastCancer.rds")
TNBC_samples <- TNBC_samples$TNBC_status
TNBC_samples <- TNBC_samples %>% filter(TNBC_status == "TRUE")
sample_info <- sample_info %>% filter(patient %in% TNBC_samples$bcr_patient_barcode)


#####


# Plot PCA with TNBC samples
rnaSeq_vst_filt <- rnaSeq_vst[row.names(rnaSeq_vst) %in% gene_list, colnames(rnaSeq_vst) %in% sample_info$barcode]

pca_res <- prcomp(t(rnaSeq_vst_filt))
pca_eigen <- factoextra::get_eig(pca_res)

pca_data <- as_tibble(pca_res$x, rownames = "sample_id")

pca_data <- pca_data %>%
  left_join(sample_info %>% dplyr::select(c("barcode", "condition", "paper_pathologic_stage")),
            by = c("sample_id" = "barcode"))  %>%
  separate(condition,
           into = c("sample_type", "cancer_subtype"),
           sep = "_") %>%
  rename("stage" = "paper_pathologic_stage")
pca_data$stage[pca_data$stage == "NA"] <- NA

plot_2 <- ggplot() +
  geom_point(data = pca_data, 
             mapping = aes(x = PC1, y = PC2, shape = stage, color = cancer_subtype),
             size = 2, 
             stroke = 0.5) +
  labs(title = "TNBC samples", 
       x = paste0("PC1: ", round(pca_eigen$variance.percent[1], 2), "% variance"),
       y = paste0("PC2: ", round(pca_eigen$variance.percent[2], 2) , "% variance"), 
       shape = "Sample type: ",
       color = "Cancer subtype: ") +
  scale_color_manual(values = c("Basal" = "#A3A500", "Her2" = "#00BF7D", 
                                "LumA" = "#00B0F6", "LumB" = "#F8766D", 
                                "Normal" = "#E76BF3"), na.value = "#A9A9A9") + 
  scale_shape_manual(values = c("Stage_I" = 0, "Stage_II" = 1, "Stage_III" = 2, "Stage_IV" = 3), na.value = 4) + 
  theme(panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.25, linetype = NULL),
        panel.grid = element_blank(),
        panel.spacing = unit(0.1, "cm"),
        text = element_text(size = 8),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 8),
        axis.text.x = element_text(angle = 0, vjust = 0, hjust = 0.5),
        axis.ticks = element_line(colour = "black", linewidth = 0.2),
        legend.position = "bottom",
        legend.key = element_rect(fill = NA), 
        legend.key.size = unit(0.5, "cm"),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.margin = margin(1,1,1,1),
        legend.spacing = unit(0, "cm")  ) 

tiff(paste0("IOFiles/SSN_background_genes/plot_PCA_SSN_bkgGene__DA_TCGA__DT_BreastCancer__DS_TNBC.tiff"),
     width = 18, height = 15,
     units = "cm", compression = "lzw", res = 1200)

print(plot_2)

dev.off()


#####


print(warnings())