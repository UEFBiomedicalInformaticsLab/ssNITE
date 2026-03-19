set.seed(5081)


# Script to process the RNA-Seq data

# Load libraries
library(unixtools)
library(DESeq2)
library(edgeR)
library(tidyverse)
library(DDoutlier)


# Set temporary directory
if(!dir.exists("tmp_dir/")){dir.create("tmp_dir/", recursive = TRUE)}
set.tempdir("tmp_dir/")


#####


# Read the TCGA data
TCGA_data <- readRDS("IOFiles/GDC_data/TCGABRCA_data.rds")

rnaSeq_counts <- TCGA_data$rnaSeq_counts
sample_info <- TCGA_data$rnaSeq_counts_sampleInfo %>% column_to_rownames("sample_id")
gene_info <- TCGA_data$rnaSeq_counts_geneInfo


# Create a condition by merging sample type and BRCA sub-type info
# Needed as samples are marked with NA in sub-type
sample_info$condition <- paste(gsub(" ", "", sample_info$sample_type),
                               sample_info$paper_BRCA_Subtype_PAM50, sep = "_")
sample_info$condition <- as.factor(sample_info$condition)


# Filter to keep only protein coding gene counts
gene_info <- gene_info %>% filter(gene_type == "protein_coding")
rnaSeq_counts <- rnaSeq_counts %>% 
  filter(ensembl_gene_id %in% gene_info$ensembl_gene_id) %>%
  column_to_rownames("ensembl_gene_id")


# Filter low count genes
# Filtering genes with CPM less than 10 in n samples
# n is the number of samples in the smallest group

keep <- row.names(rnaSeq_counts)[ rowSums(cpm(rnaSeq_counts)>10) > min(table(sample_info$condition)) ] 
rnaSeq_counts <- rnaSeq_counts[row.names(rnaSeq_counts) %in% keep, ]


# Rearrange count matrix as in sample info 
if(!all(sample_info$barcode == colnames(rnaSeq_counts))){
  rnaSeq_counts <- rnaSeq_counts[, sample_info$barcode]
}


#####


# Import data to DESeq2 
dds <- DESeqDataSetFromMatrix(countData = rnaSeq_counts, 
                              colData = sample_info, 
                              design = ~ condition)

dds$condition <- relevel(dds$condition, ref = "SolidTissueNormal_NA")

dds <- DESeq(dds)


# Transform data data using VST
vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
rnaSeq_vst <- assay(vsd)


#####


# Check for batch effect
pca_res <- prcomp(t(rnaSeq_vst))
pca_eigen <- factoextra::get_eig(pca_res)
pca_data <- pca_res$x %>% 
  as.data.frame() %>% 
  rownames_to_column("sample_id") %>% 
  left_join(sample_info %>% select(c("barcode", "sample_type", "paper_BRCA_Subtype_PAM50")), 
            by = c("sample_id" = "barcode")) %>%
  rename("cancer_subtype" = "paper_BRCA_Subtype_PAM50") %>% 
  column_to_rownames("sample_id")


# Identify outliers
LOF_scores <- LOF(pca_res$x[, 1:2], 
                  k = max(5, round(sqrt(length(sample_info$barcode)))) )
outlier_samples <- row.names(pca_res$x)[which(( LOF_scores - mean(LOF_scores) ) > sd(LOF_scores) * 5)]


# Plot

plot_list <- list()


plot_list[["A"]] <- ggplot() +
  geom_point(data = pca_data, 
             mapping = aes(x = PC1, y = PC2, shape = sample_type, color = cancer_subtype),
             size = 0.5, 
             stroke = 0.1) +
  geom_point(data = pca_data[row.names(pca_data) %in% outlier_samples, ], 
             mapping = aes(x = PC1, y = PC2),
             size = 1, 
             shape = 5,
             stroke = 0.1,
             color = "red") +
  labs(title = "All samples", 
       x = paste0("PC1: ", round(pca_eigen$variance.percent[1], 2), "% variance"),
       y = paste0("PC2: ", round(pca_eigen$variance.percent[2], 2) , "% variance"), 
       shape = "Sample type: ",
       color = "Cancer subtype:")+
  scale_color_manual(values = c("Basal" = "#A3A500", "Her2" = "#00BF7D", 
                                "LumA" = "#00B0F6", "LumB" = "#F8766D", 
                                "Normal" = "#E76BF3"), na.value = "#A9A9A9") + 
  theme(panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.25, linetype = NULL),
        panel.grid = element_blank(),
        panel.spacing = unit(0.1, "cm"),
        text = element_text(size = 4),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 4),
        axis.text.x = element_text(angle = 0, vjust = 0, hjust = 0.5),
        axis.ticks = element_line(colour = "black", linewidth = 0.2),
        legend.position = "bottom",
        legend.key = element_rect(fill = NA), 
        legend.key.size = unit(0.25, "cm"),
        legend.title = element_text(size = 2.5),
        legend.text = element_text(size = 2),
        legend.margin = margin(1,1,1,1),
        legend.spacing = unit(0, "cm")  ) 


#####


# Reprocess the data after filtering the samples


# Re-read the data
rnaSeq_counts_filt <- TCGA_data$rnaSeq_counts
sample_info_filt <- TCGA_data$rnaSeq_counts_sampleInfo %>% column_to_rownames("sample_id")
gene_info_filt <- TCGA_data$rnaSeq_counts_geneInfo


# Create a condition by merging sample type and BRCA sub-type info
# Needed as samples are marked with NA in sub-type
sample_info_filt$condition <- paste(gsub(" ", "", sample_info_filt$sample_type),
                                    sample_info_filt$paper_BRCA_Subtype_PAM50, sep = "_")
sample_info_filt$condition <- as.factor(sample_info_filt$condition)


# Filter outliers
sample_info_filt <- sample_info_filt[!sample_info_filt$barcode %in% outlier_samples, ]


# Remove primary tumour samples without labelled sub-groups and normal-like subtype
sample_info_filt <- sample_info_filt[!sample_info_filt$condition %in% c("PrimaryTumor_NA", "PrimaryTumor_Normal"), ]
sample_info_filt$condition <- droplevels(sample_info_filt$condition)


# Extract the raw counts
rnaSeq_counts_filt <- rnaSeq_counts_filt %>% 
  select(c("ensembl_gene_id", sample_info_filt$barcode)) %>%
  filter(ensembl_gene_id %in% gene_info$ensembl_gene_id) %>%
  column_to_rownames("ensembl_gene_id")


# Filter low count genes
# Filtering genes with CPM less than 1 in n samples
# n is the number of samples in the smallest group

keep <- row.names(rnaSeq_counts_filt)[ rowSums(cpm(rnaSeq_counts_filt)>10) > min(table(sample_info_filt$condition)) ] 
rnaSeq_counts_filt <- rnaSeq_counts_filt[row.names(rnaSeq_counts_filt) %in% keep, ]


# Rearrange count matrix as in sample info 
if(!all(sample_info_filt$barcode == colnames(rnaSeq_counts_filt))){
  rnaSeq_counts_filt <- rnaSeq_counts_filt[, sample_info_filt$barcode]
}


# Import data to DESeq2 
dds <- DESeqDataSetFromMatrix(countData = rnaSeq_counts_filt, 
                              colData = sample_info_filt, 
                              design = ~ condition)

dds$condition <- relevel(dds$condition, ref = "SolidTissueNormal_NA")

dds <- DESeq(dds)


# Transform data data using VST
vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
rnaSeq_vst <- assay(vsd)


# Check for batch effect
pca_res <- prcomp(t(rnaSeq_vst))
pca_eigen <- factoextra::get_eig(pca_res)
pca_data <- pca_res$x %>% 
  as.data.frame() %>% 
  rownames_to_column("sample_id") %>% 
  left_join(sample_info_filt %>% select(c("barcode", "sample_type", "paper_BRCA_Subtype_PAM50")), 
            by = c("sample_id" = "barcode")) %>%
  rename("cancer_subtype" = "paper_BRCA_Subtype_PAM50") %>% 
  column_to_rownames("sample_id")


# # Identify outliers
# LOF_scores <- LOF(pca_res$x[, 1:2], 
#                   k = max(5, round(sqrt(length(sample_info$barcode)))) )
# outlier_samples <- row.names(pca_res$x)[which(( LOF_scores - mean(LOF_scores) ) > sd(LOF_scores) * 5)]


# Plot
plot_list[["B"]] <- ggplot() +
  geom_point(data = pca_data, 
             mapping = aes(x = PC1, y = PC2, shape = sample_type, color = cancer_subtype),
             size = 0.5, 
             stroke = 0.1) +
  geom_point(data = pca_data[row.names(pca_data) %in% outlier_samples, ], 
             mapping = aes(x = PC1, y = PC2),
             size = 1, 
             shape = 5,
             stroke = 0.1,
             color = "red") +
  labs(title = "Filtered samples", 
       x = paste0("PC1: ", round(pca_eigen$variance.percent[1], 2), "% variance"),
       y = paste0("PC2: ", round(pca_eigen$variance.percent[2], 2) , "% variance"), 
       shape = "Sample type: ",
       color = "Cancer subtype: ") +
  scale_color_manual(values = c("Basal" = "#A3A500", "Her2" = "#00BF7D", 
                                "LumA" = "#00B0F6", "LumB" = "#F8766D", 
                                "Normal" = "#E76BF3"), na.value = "#A9A9A9") + 
  theme(panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.25, linetype = NULL),
        panel.grid = element_blank(),
        panel.spacing = unit(0.1, "cm"),
        text = element_text(size = 4),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 4),
        axis.text.x = element_text(angle = 0, vjust = 0, hjust = 0.5),
        axis.ticks = element_line(colour = "black", linewidth = 0.2),
        legend.position = "bottom",
        legend.key = element_rect(fill = NA), 
        legend.key.size = unit(0.25, "cm"),
        legend.title = element_text(size = 2.5),
        legend.text = element_text(size = 2),
        legend.margin = margin(1,1,1,1),
        legend.spacing = unit(0, "cm")  ) 


#####


# Save plot to file
if(!dir.exists("IOFiles/DE_analysis/")){
  dir.create("IOFiles/DE_analysis/", recursive = TRUE)
}

tiff("IOFiles/DE_analysis/plot_PCA__DA_TCGA__DT_BreastCancer.tiff",
     width = 15, height = 8,
     units = "cm", compression = "lzw", res = 1200)

ggpubr::ggarrange(plotlist = plot_list, ncol = 2, common.legend = TRUE, legend = "bottom")

dev.off()


#####


# Save the normalised data

if(!dir.exists("IOFiles/GDC_data/")){dir.create("IOFiles/GDC_data/", recursive = TRUE)}

saveRDS(list("rnaSeq_vst" = as.data.frame(rnaSeq_vst), 
             "sample_info" = sample_info_filt, 
             "gene_info" = gene_info_filt), 
        file = "IOFiles/GDC_data/vstCounts__DA_TCGA__DT_BreastCancer.rds")


#####


# Perform DE analysis
DE_res <- list()
DE_res[["LumA"]] <- results(dds, contrast = c("condition", "PrimaryTumor_LumA", "SolidTissueNormal_NA"))
DE_res[["LumB"]] <- results(dds, contrast = c("condition", "PrimaryTumor_LumB", "SolidTissueNormal_NA"))
DE_res[["Basal"]] <- results(dds, contrast = c("condition", "PrimaryTumor_Basal", "SolidTissueNormal_NA"))
DE_res[["Her2"]] <- results(dds, contrast = c("condition", "PrimaryTumor_Her2", "SolidTissueNormal_NA"))

# lapply(DE_res, function(x){ nrow(x[abs(x$log2FoldChange) > 2 & x$padj < 0.001, ]) })

if(!dir.exists("IOFiles/DE_analysis/")){dir.create("IOFiles/DE_analysis/", recursive = TRUE)}
saveRDS(DE_res, "IOFiles/DE_analysis/DEA__DA_TCGA__DT_BreastCancer.rds")


#####


print(warnings())