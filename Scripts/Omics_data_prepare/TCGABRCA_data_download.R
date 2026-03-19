set.seed(5081)


# Script to download TCGA data


# Load libraries
library(unixtools)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(biomaRt)
library(tidyverse)

# Set temporary directory
if(!dir.exists("tmp_dir/")){dir.create("tmp_dir/", recursive = TRUE)}
set.tempdir("tmp_dir/")

if(!dir.exists("Databases/Genomic_Data_Commons/")){dir.create("Databases/Genomic_Data_Commons/", recursive = TRUE)}


#####


# Download the RNA-Seq data

# Create query for downloading the data
query <- GDCquery(project = "TCGA-BRCA",
                  data.category = "Transcriptome Profiling", 
                  data.type = "Gene Expression Quantification",
                  sample.type = c("Primary Tumor", "Solid Tissue Normal"), 
                  experimental.strategy = "RNA-Seq", 
                  workflow.type = "STAR - Counts")

# Download the data
GDCdownload(query = query, 
            directory = "Databases/Genomic_Data_Commons/", 
            method = "api", 
            files.per.chunk = 10)

# Prepare the download files 
data <- GDCprepare(query = query, 
                   save = FALSE, 
                   directory = "Databases/Genomic_Data_Commons/")

rnaSeq_counts <- as_tibble(assay(data, i = "unstranded"), rownames = "gene_id")
rnaSeq_counts_sampleInfo <- as_tibble(colData(data), rownames = "sample_id")
rnaSeq_counts_geneInfo <- as_tibble(rowRanges(data)) %>% select(c("gene_id", everything()))


#####


# Add Ensembl Gene IDs (without version) to the RNA-Seq data

# Extract the gene IDs without the version numbers
# TCGA uses GENCODE version 36 which is equivalent to Ensembl version 102
# Ref: https://gdc.cancer.gov/about-data/gdc-data-processing/gdc-reference-files
# Ref: https://www.gencodegenes.org/human/releases.html

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", version = "102")
ensemblGeneIDver_2_ensemblGeneID <- getBM(attributes = c("ensembl_gene_id", "ensembl_gene_id_version"), 
                                          mart = ensembl, 
                                          filters = "ensembl_gene_id_version", 
                                          values = rnaSeq_counts_geneInfo$gene_id)

# 44 PAR genes not mapped

# Filter to keep only the genes with valid ID
rnaSeq_counts_geneInfo <- rnaSeq_counts_geneInfo %>% 
  left_join(ensemblGeneIDver_2_ensemblGeneID, by = c("gene_id" = "ensembl_gene_id_version")) %>% 
  filter(!is.na(ensembl_gene_id)) %>% 
  dplyr::select(-c("gene_id"))


rnaSeq_counts <- rnaSeq_counts %>% 
  left_join(ensemblGeneIDver_2_ensemblGeneID, by = c("gene_id" = "ensembl_gene_id_version")) %>% 
  filter(!is.na(ensembl_gene_id)) %>% 
  dplyr::select(c("gene_id", "ensembl_gene_id", everything())) %>% 
  dplyr::select(-c("gene_id"))


#####


# Download the clinical information 

query <- GDCquery(project = "TCGA-BRCA", 
                  data.category = "Clinical", 
                  data.type = "Clinical Supplement",
                  data.format = "bcr xml")

GDCdownload(query = query, 
            directory = "Databases/Genomic_Data_Commons/", 
            method = "api", 
            files.per.chunk = 10)

clinicalInfo_drug <- GDCprepare_clinic(query = query, 
                                       clinical.info = c("drug"), 
                                       directory = "Databases/Genomic_Data_Commons/")

clinicalInfo_followUp <- GDCprepare_clinic(query = query, 
                                           clinical.info = c("follow_up"), 
                                           directory = "Databases/Genomic_Data_Commons/")

clinicalInfo_radiation <- GDCprepare_clinic(query = query, 
                                            clinical.info = c("radiation"), 
                                            directory = "Databases/Genomic_Data_Commons/")

clinicalInfo_patient <- GDCprepare_clinic(query = query, 
                                          clinical.info = c("patient"), 
                                          directory = "Databases/Genomic_Data_Commons/")


clinicalInfo_stageEvent <- GDCprepare_clinic(query = query, 
                                             clinical.info = c("stage_event"), 
                                             directory = "Databases/Genomic_Data_Commons/")

clinicalInfo_newTumorEvent <- GDCprepare_clinic(query = query, 
                                                clinical.info = c("new_tumor_event"), 
                                                directory = "Databases/Genomic_Data_Commons/")


#####


# Download the proteomic data

# Create query for downloading the data
query <- GDCquery(project = "TCGA-BRCA",
                  data.category = "Proteome Profiling", 
                  data.type = "Protein Expression Quantification")

# Download the data
GDCdownload(query = query, 
            directory = "Databases/Genomic_Data_Commons/", 
            method = "api", 
            files.per.chunk = 10)

# Prepare the download files 
proteome_profile <- GDCprepare(query = query, 
                               directory = "Databases/Genomic_Data_Commons/")


#####


# Add annotations to the proteomic data 

# Download the annotation file from NCI website // GDC Reference Files
if(!file.exists("Databases/Genomic_Data_Commons/TCGA_antibodies_descriptions.gencode.v36.tsv")){
  download.file(url = "https://api.gdc.cancer.gov/v0/data/62647302-b4d3-4a81-a7c0-d141f5dbd300",
                destfile = "Databases/Genomic_Data_Commons/TCGA_antibodies_descriptions.gencode.v36.tsv", 
                mode = "wb", method = "wget")
}

TCGA_antibodies <- read.table("Databases/Genomic_Data_Commons/TCGA_antibodies_descriptions.gencode.v36.tsv", 
                              header = TRUE, row.names = NULL, sep = "\t")

TCGA_antibodies <- TCGA_antibodies %>% 
  separate_rows(gene_name, gene_id, sep = "/") %>% 
  filter(str_detect(gene_id, "^ENSG"))


# Add Ensembl Gene IDs (without version) to the peptide annotations

# Extract the gene IDs without the version numbers
# TCGA uses GENCODE version 36 which is equivalent to Ensembl version 102
# Ref: https://gdc.cancer.gov/about-data/gdc-data-processing/gdc-reference-files
# Ref: https://www.gencodegenes.org/human/releases.html

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", version = "102")
ensemblGeneIDver_2_ensemblGeneID <- getBM(attributes = c("ensembl_gene_id", "ensembl_gene_id_version"), 
                                          mart = ensembl, 
                                          filters = "ensembl_gene_id_version", 
                                          values = unique(TCGA_antibodies$gene_id))

# Filter to keep only the genes with valid ID
TCGA_antibodies <- TCGA_antibodies %>% 
  left_join(ensemblGeneIDver_2_ensemblGeneID, by = c("gene_id" = "ensembl_gene_id_version")) %>% 
  filter(!is.na(ensembl_gene_id)) 

proteome_profile <- proteome_profile %>% 
  left_join(TCGA_antibodies %>% dplyr::select(c("AGID", "ensembl_gene_id")), 
            by = "AGID") %>% 
  dplyr::select(c("AGID", "ensembl_gene_id", everything())) %>% 
  filter(!is.na(ensembl_gene_id)) 


#####


# Download the mutation profile

# Create query for downloading the data
query <- GDCquery(project = "TCGA-BRCA",
                  data.category = "Simple Nucleotide Variation", 
                  data.type = "Masked Somatic Mutation", 
                  access = "open")

# Download the data
GDCdownload(query = query, 
            directory = "Databases/Genomic_Data_Commons/", 
            method = "api", 
            files.per.chunk = 10)

# Prepare the download files 
SNV_profile <- GDCprepare(query = query, 
                          directory = "Databases/Genomic_Data_Commons/")


#####


# Download the copy number variation data (gene level)

# Create query for downloading the data
query <- GDCquery(project = "TCGA-BRCA",
                  data.category = "Copy Number Variation", 
                  data.type = "Gene Level Copy Number", 
                  workflow.type = "ABSOLUTE LiftOver",
                  sample.type = c("Primary Tumor", "Solid Tissue Normal"),
                  access = "open")

# Download the data
GDCdownload(query = query, 
            directory = "Databases/Genomic_Data_Commons/", 
            method = "api", 
            files.per.chunk = 10)

# Prepare the download files 
data <- GDCprepare(query = query, 
                   directory = "Databases/Genomic_Data_Commons/")

CNV_geneLevel <- assays(data)
CNV_geneLevel_geneInfo <- as_tibble(rowData(data))
CNV_geneLevel_sampleInfo <- as_tibble(colData(data))


#####


# Add Ensembl Gene IDs (without version) to CNV data

# Extract the gene IDs without the version numbers
# TCGA uses GENCODE version 36 which is equivalent to Ensembl version 102
# Ref: https://gdc.cancer.gov/about-data/gdc-data-processing/gdc-reference-files
# Ref: https://www.gencodegenes.org/human/releases.html

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", version = "102")
ensemblGeneIDver_2_ensemblGeneID <- getBM(attributes = c("ensembl_gene_id", "ensembl_gene_id_version"), 
                                          mart = ensembl, 
                                          filters = "ensembl_gene_id_version", 
                                          values = unique(CNV_geneLevel_geneInfo$gene_id))

# Filter to keep only the genes with valid ID
CNV_geneLevel_geneInfo <- CNV_geneLevel_geneInfo %>% 
  left_join(ensemblGeneIDver_2_ensemblGeneID, by = c("gene_id" = "ensembl_gene_id_version")) %>% 
  filter(!is.na(ensembl_gene_id)) 

# Add IDs to the data
CNV_geneLevel <- lapply(CNV_geneLevel, function(x){  
  x %>% as_tibble(rownames = "gene_id") %>% 
    left_join(ensemblGeneIDver_2_ensemblGeneID, 
              by = c("gene_id" = "ensembl_gene_id_version")) %>% 
    select(c("ensembl_gene_id", everything())) %>% 
    filter(!is.na(ensembl_gene_id)) 
})


#####


# Download the methylation data

# Create query for downloading the data
query <- GDCquery(project = "TCGA-BRCA",
                  data.category = "DNA Methylation", 
                  data.type = "Methylation Beta Value", 
                  platform = "Illumina Human Methylation 450",
                  sample.type = c("Primary Tumor", "Solid Tissue Normal"),
                  access = "open")

# Download the data
GDCdownload(query = query, 
            directory = "Databases/Genomic_Data_Commons/", 
            method = "api", 
            files.per.chunk = 10)

# Prepare the download files 
data <- GDCprepare(query = query, 
                   directory = "Databases/Genomic_Data_Commons/")

methylation_betaVals <- as_tibble(assay(data), rownames = "cpg_site")
methylation_betaVals_geneInfo <- as_tibble(rowData(data), rownames = "cpg_site")
methylation_betaVals_sampleInfo <- as_tibble(colData(data))


#####


# Get annotations to the methylation data 

# Download the annotation file from NCI website // GDC Reference Files
if(!file.exists("Databases/Genomic_Data_Commons/HM450.hg38.manifest.gencode.v36.tsv.gz")){
  download.file(url = "https://api.gdc.cancer.gov/v0/data/021a2330-951d-474f-af24-1acd77e7664f",
                destfile = "Databases/Genomic_Data_Commons/HM450.hg38.manifest.gencode.v36.tsv.gz", 
                mode = "wb", method = "wget")
}

methylation_array_annot <- read.table("Databases/Genomic_Data_Commons/HM450.hg38.manifest.gencode.v36.tsv.gz", 
                                      header = TRUE, row.names = NULL, sep = "\t")


#####


# Save the data
if(!dir.exists("IOFiles/GDC_data/")){dir.create("IOFiles/GDC_data/", recursive = TRUE)}


tmp1 <- list("rnaSeq_counts" = rnaSeq_counts, 
             "rnaSeq_counts_geneInfo" = rnaSeq_counts_geneInfo, 
             "rnaSeq_counts_sampleInfo" = rnaSeq_counts_sampleInfo, 
             
             "clinicalInfo_drug" = clinicalInfo_drug, 
             "clinicalInfo_followUp" = clinicalInfo_followUp, 
             "clinicalInfo_radiation" = clinicalInfo_radiation, 
             "clinicalInfo_patient" = clinicalInfo_patient, 
             "clinicalInfo_stageEvent" = clinicalInfo_stageEvent, 
             "clinicalInfo_newTumorEvent" = clinicalInfo_newTumorEvent, 
             
             "proteome_profile" = proteome_profile, 
             
             "SNV_profile" = SNV_profile, 
             
             "CNV_geneLevel" = CNV_geneLevel, 
             "CNV_geneLevel_geneInfo" = CNV_geneLevel_geneInfo, 
             "CNV_geneLevel_sampleInfo" = CNV_geneLevel_sampleInfo, 
             
             "methylation_betaVals" = methylation_betaVals, 
             "methylation_betaVals_geneInfo" = methylation_betaVals_geneInfo, 
             "methylation_betaVals_sampleInfo" = methylation_betaVals_sampleInfo, 
             "methylation_array_annot" = methylation_array_annot)


saveRDS(tmp1, "IOFiles/GDC_data/TCGABRCA_data.rds")


#####


print(warnings())