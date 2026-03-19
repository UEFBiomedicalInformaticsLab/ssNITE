set.seed(5081)


# Script to prepare gene-disease associations from different databases
# Notes:
# (a) Disease-gene associations are stored as R named list 
# (b) Only diseases with more than 5 genes retained


# Load libraries 
library(unixtools)
library(org.Hs.eg.db)
library(OmnipathR)
library(tidyverse)
library(arrow)
library(openxlsx)

# source("Scripts/Functions/Functions_data_manipulation.R")
source("Scripts/Functions/Functions_ID_Conversion.R")


# Set temporary directory
if(!dir.exists("tmp_dir/")){dir.create("tmp_dir/", recursive = TRUE)}
set.tempdir("tmp_dir/")


#####


# OpenTargets  

# Different types explained: https://platform-docs.opentargets.org/evidence

print("Processing OpenTargets")

# Download data using wget in terminal
if(!dir.exists("Databases/OpenTargets/")){dir.create("Databases/OpenTargets/", recursive = TRUE)}

if(!dir.exists("Databases/OpenTargets/associationByDatatypeDirect")){
  system("wget --recursive --no-parent --no-host-directories -P Databases/OpenTargets/associationByDatatypeDirect/ --cut-dirs 8 ftp://ftp.ebi.ac.uk/pub/databases/opentargets/platform/25.06/output/association_by_datatype_direct")
}

if(!dir.exists("Databases/OpenTargets/diseases")){
  system("wget --recursive --no-parent --no-host-directories -P Databases/OpenTargets/diseases/ --cut-dirs 8 ftp://ftp.ebi.ac.uk/pub/databases/opentargets/platform/25.06/output/disease")
}


# Read the disease info
OpenTargets_Diseases <- open_dataset("Databases/OpenTargets/diseases/")
OpenTargets_Diseases <- OpenTargets_Diseases %>% select(id, name, description)
OpenTargets_Diseases <- OpenTargets_Diseases %>% collect()


# Read the drug-gene associations
OpenTargets_Target_Disease <- open_dataset("Databases/OpenTargets/associationByDatatypeDirect/")


## Disease gene enrichment library based genetic association (GA) 
OpenTargets_Target_Disease_GA <- OpenTargets_Target_Disease %>%
  filter(datatypeId == "genetic_association")

OpenTargets_Target_Disease_GA <- OpenTargets_Target_Disease_GA %>% collect()  
OpenTargets_Target_Disease_GA <- OpenTargets_Target_Disease_GA[OpenTargets_Target_Disease_GA$score > quantile(unique(OpenTargets_Target_Disease_GA$score), 0.75),]

OpenTargets_Target_Disease_GA <- OpenTargets_Target_Disease_GA %>% 
  left_join(OpenTargets_Diseases %>% select(c("id", "name")), 
            by = c("diseaseId" = "id")) %>% 
  rename(c("diseaseLabel" = "name"))

OpenTargets_Target_Disease_GA$diseaseLabel <- paste0(OpenTargets_Target_Disease_GA$diseaseLabel, " (", OpenTargets_Target_Disease_GA$diseaseId, ")")
OpenTargets_Disease2Gene_GA_lib <- split(OpenTargets_Target_Disease_GA$targetId, OpenTargets_Target_Disease_GA$diseaseLabel)
OpenTargets_Disease2Gene_GA_lib <- lapply(OpenTargets_Disease2Gene_GA_lib, unique)

OpenTargets_Disease2Gene_GA_lib <- OpenTargets_Disease2Gene_GA_lib[lengths(OpenTargets_Disease2Gene_GA_lib) >= 5]

if(!dir.exists("IOFiles/Geneset_libraries/")){dir.create("IOFiles/Geneset_libraries/", recursive = TRUE)}
saveRDS(OpenTargets_Disease2Gene_GA_lib, "IOFiles/Geneset_libraries/OpenTargets_Disease2Gene_GA_lib.rds")


## Disease gene enrichment library based RNA expression score 
OpenTargets_Target_Disease_RNA <- OpenTargets_Target_Disease %>%
  filter(datatypeId == "rna_expression") 

OpenTargets_Target_Disease_RNA <- OpenTargets_Target_Disease_RNA %>% collect()  
OpenTargets_Target_Disease_RNA <- OpenTargets_Target_Disease_RNA[OpenTargets_Target_Disease_RNA$score > quantile(unique(OpenTargets_Target_Disease_RNA$score), 0.75),]

OpenTargets_Target_Disease_RNA <- OpenTargets_Target_Disease_RNA %>% 
  left_join(OpenTargets_Diseases %>% select(c("id", "name")), 
            by = c("diseaseId" = "id")) %>% 
  rename(c("diseaseLabel" = "name"))

OpenTargets_Target_Disease_RNA$diseaseLabel <- paste0(OpenTargets_Target_Disease_RNA$diseaseLabel, " (", OpenTargets_Target_Disease_RNA$diseaseId, ")")
OpenTargets_Disease2Gene_RNA_lib <- split(OpenTargets_Target_Disease_RNA$targetId, OpenTargets_Target_Disease_RNA$diseaseLabel)
OpenTargets_Disease2Gene_RNA_lib <- lapply(OpenTargets_Disease2Gene_RNA_lib, unique)

OpenTargets_Disease2Gene_RNA_lib <- OpenTargets_Disease2Gene_RNA_lib[lengths(OpenTargets_Disease2Gene_RNA_lib) >= 5]

if(!dir.exists("IOFiles/Geneset_libraries/")){dir.create("IOFiles/Geneset_libraries/", recursive = TRUE)}
saveRDS(OpenTargets_Disease2Gene_RNA_lib, "IOFiles/Geneset_libraries/OpenTargets_Disease2Gene_RNA_lib.rds")


## Disease gene enrichment library based literature 
OpenTargets_Target_Disease_lit <- OpenTargets_Target_Disease %>%
  filter(datatypeId == "literature") 

OpenTargets_Target_Disease_lit <- OpenTargets_Target_Disease_lit %>% collect()  
OpenTargets_Target_Disease_lit <- OpenTargets_Target_Disease_lit[OpenTargets_Target_Disease_lit$score > quantile(unique(OpenTargets_Target_Disease_lit$score), 0.75),]

OpenTargets_Target_Disease_lit <- OpenTargets_Target_Disease_lit %>% 
  left_join(OpenTargets_Diseases %>% select(c("id", "name")), 
            by = c("diseaseId" = "id")) %>% 
  rename(c("diseaseLabel" = "name"))

OpenTargets_Target_Disease_lit$diseaseLabel <- paste0(OpenTargets_Target_Disease_lit$diseaseLabel, " (", OpenTargets_Target_Disease_lit$diseaseId, ")")
OpenTargets_Disease2Gene_lit_lib <- split(OpenTargets_Target_Disease_lit$targetId, OpenTargets_Target_Disease_lit$diseaseLabel)
OpenTargets_Disease2Gene_lit_lib <- lapply(OpenTargets_Disease2Gene_lit_lib, unique)

OpenTargets_Disease2Gene_lit_lib <- OpenTargets_Disease2Gene_lit_lib[lengths(OpenTargets_Disease2Gene_lit_lib) >= 5]

if(!dir.exists("IOFiles/Geneset_libraries/")){dir.create("IOFiles/Geneset_libraries/", recursive = TRUE)}
saveRDS(OpenTargets_Disease2Gene_lit_lib, "IOFiles/Geneset_libraries/OpenTargets_Disease2Gene_lit_lib.rds")


## Disease gene enrichment library based on somatic mutation 
OpenTargets_Target_Disease_SM <- OpenTargets_Target_Disease %>%
  filter(datatypeId == "somatic_mutation") 

OpenTargets_Target_Disease_SM <- OpenTargets_Target_Disease_SM %>% collect()  
OpenTargets_Target_Disease_SM <- OpenTargets_Target_Disease_SM[OpenTargets_Target_Disease_SM$score > quantile(unique(OpenTargets_Target_Disease_SM$score), 0.75),]

OpenTargets_Target_Disease_SM <- OpenTargets_Target_Disease_SM %>% 
  left_join(OpenTargets_Diseases %>% select(c("id", "name")), 
            by = c("diseaseId" = "id")) %>% 
  rename(c("diseaseLabel" = "name"))

OpenTargets_Target_Disease_SM$diseaseLabel <- paste0(OpenTargets_Target_Disease_SM$diseaseLabel, " (", OpenTargets_Target_Disease_SM$diseaseId, ")")
OpenTargets_Disease2Gene_SM_lib <- split(OpenTargets_Target_Disease_SM$targetId, OpenTargets_Target_Disease_SM$diseaseLabel)
OpenTargets_Disease2Gene_SM_lib <- lapply(OpenTargets_Disease2Gene_SM_lib, unique)

OpenTargets_Disease2Gene_SM_lib <- OpenTargets_Disease2Gene_SM_lib[lengths(OpenTargets_Disease2Gene_SM_lib) >= 5]

if(!dir.exists("IOFiles/Geneset_libraries/")){dir.create("IOFiles/Geneset_libraries/", recursive = TRUE)}
saveRDS(OpenTargets_Disease2Gene_SM_lib, "IOFiles/Geneset_libraries/OpenTargets_Disease2Gene_SM_lib.rds")


#####


# IntOGen

# Retrieve cancer driver gene list from https://intogen.org/download
# Download and unzip files from: 
# (a) https://intogen.org/download?file=IntOGen-Drivers-20200201.zip
# (b) https://intogen.org/download?file=IntOGen-Cohorts-20200201.zip

print("Processing IntOGen")

if(!dir.exists("Databases/IntOGen")){dir.create("Databases/IntOGen", recursive = TRUE)}
if(!file.exists("Databases/IntOGen/Compendium_Cancer_Genes.tsv")){
  download.file(url = "https://www.intogen.org/download?file=IntOGen-Drivers-20240920.zip", 
                destfile = "Databases/IntOGen/IntOGen-Drivers-20240920.zip", method = "wget")
  unzip("Databases/IntOGen/IntOGen-Drivers-20240920.zip", files = "2024-06-18_IntOGen-Drivers/Compendium_Cancer_Genes.tsv", 
        exdir = "Databases/IntOGen/", junkpaths = TRUE)
}

if(!file.exists("Databases/IntOGen/cohorts.tsv")){
  download.file(url = "https://www.intogen.org/download?file=IntOGen-Cohorts-20240920.zip",
                destfile = "Databases/IntOGen/IntOGen-Cohorts-20240920.zip", method = "wget")
  unzip("Databases/IntOGen/IntOGen-Cohorts-20240920.zip", files = "2024-06-18_IntOGen-Cohorts/cohorts.tsv", 
        exdir = "Databases/IntOGen/", junkpaths = TRUE)
}


# Read the cancer types in IntOGen
IntOGen_cancers <- read.table("Databases/IntOGen/cohorts.tsv", 
                              sep = "\t", 
                              header = TRUE, 
                              check.names = FALSE, 
                              quote = "")
IntOGen_cancers <- unique(IntOGen_cancers[, c("CANCER", "CANCER_NAME")])


# Read disease gene association data
IntOGen_data <- read.table("Databases/IntOGen/Compendium_Cancer_Genes.tsv", 
                           sep = "\t", 
                           header = TRUE, 
                           check.names = FALSE)


# Add ensembl gene IDs
geneSymbol_2_ensemblId <- AnnotationDbi::select(org.Hs.eg.db, 
                                                keys = unique(IntOGen_data$SYMBOL), 
                                                columns = "ENSEMBL", 
                                                keytype = "SYMBOL")

IntOGen_data <- IntOGen_data %>% 
  left_join(geneSymbol_2_ensemblId, 
            by = "SYMBOL", 
            relationship = "many-to-many") %>% 
  filter(!is.na(ENSEMBL))
rm(geneSymbol_2_ensemblId)


# Add cancer names
IntOGen_data <- IntOGen_data %>% left_join(IntOGen_cancers, by = c("CANCER_TYPE" = "CANCER"))

IntOGen_data$CANCER_TYPE <- paste0(IntOGen_data$CANCER_NAME, " (", IntOGen_data$CANCER_TYPE, ")")

IntOGen_Disease2Gene_lib <- split(IntOGen_data$ENSEMBL, f = IntOGen_data$CANCER_TYPE)
IntOGen_Disease2Gene_lib <- lapply(IntOGen_Disease2Gene_lib, unique)

IntOGen_Disease2Gene_lib <- IntOGen_Disease2Gene_lib[lengths(IntOGen_Disease2Gene_lib) >= 5]

if(!dir.exists("IOFiles/Geneset_libraries/")){dir.create("IOFiles/Geneset_libraries/", recursive = TRUE)}
saveRDS(IntOGen_Disease2Gene_lib, "IOFiles/Geneset_libraries/IntOGen_Disease2Gene_lib.rds")


#####


# PharmGKB 

print("Processing PharmGKB")

if(!dir.exists("Databases/PharmGKB/")){dir.create("Databases/PharmGKB/", recursive = TRUE)}
if(!file.exists("Databases/PharmGKB/relationships.tsv")){
  download.file(url = "https://api.pharmgkb.org/v1/download/file/data/relationships.zip",
                destfile = "Databases/PharmGKB/relationships.zip", method = "wget")
  unzip("Databases/PharmGKB/relationships.zip", files = "relationships.tsv", exdir = "Databases/PharmGKB/")
}

PharmGKB_data <- read.table("Databases/PharmGKB/relationships.tsv", 
                            sep = "\t", 
                            fill = TRUE, 
                            header = TRUE, 
                            quote = "")
PharmGKB_Gene_Disease <- PharmGKB_data[(PharmGKB_data$Entity1_type == "Gene" & PharmGKB_data$Entity2_type == "Disease" & PharmGKB_data$Association == "associated"),]

# Map PharmGKB gene ID to Ensembl gene ID
tmp1 <- func_pharmgkb_mapping(PharmGKB_Gene_Disease$Entity1_id, map_to = "ensembl_gene_id")

PharmGKB_Gene_Disease <- PharmGKB_Gene_Disease %>% 
  left_join(tmp1, 
            by = c("Entity1_id" = "Entity_id"), 
            relationship = "many-to-many") %>% 
  rename("Entity1_ensembl_gene_id" = "ensembl_gene_id")
rm(tmp1)


PharmGKB_Gene_Disease <- PharmGKB_Gene_Disease[,c("Entity1_ensembl_gene_id", "Entity2_name")]
PharmGKB_Gene_Disease[PharmGKB_Gene_Disease == ""] <- NA
PharmGKB_Gene_Disease <- na.exclude(PharmGKB_Gene_Disease)

PharmGKB_Disease2Gene_lib <- split(x = PharmGKB_Gene_Disease$Entity1_ensembl_gene_id, f = PharmGKB_Gene_Disease$Entity2_name)
PharmGKB_Disease2Gene_lib <- lapply(PharmGKB_Disease2Gene_lib, unique)

PharmGKB_Disease2Gene_lib <- PharmGKB_Disease2Gene_lib[lengths(PharmGKB_Disease2Gene_lib) >= 5]

if(!dir.exists("IOFiles/Geneset_libraries/")){dir.create("IOFiles/Geneset_libraries/", recursive = TRUE)}
saveRDS(PharmGKB_Disease2Gene_lib, "IOFiles/Geneset_libraries/PharmGKB_Disease2Gene_lib.rds")


#####


# Comparative Toxicogenimics Database (CTD) 

print("Processing Comparative Toxicogenimics Database (CTD)")

# Download and read the gene disease association file from CTD
if(!dir.exists("Databases/ComparativeToxicogenomicsDatabase/")){dir.create("Databases/ComparativeToxicogenomicsDatabase/", recursive = TRUE)}
if(!file.exists("Databases/ComparativeToxicogenomicsDatabase/CTD_genes_diseases.csv.gz")){
  download.file(url = "http://ctdbase.org/reports/CTD_genes_diseases.csv.gz",
                destfile = "Databases/ComparativeToxicogenomicsDatabase/CTD_genes_diseases.csv.gz", method = "wget")
}

CTD_Gene_Disease <- read.csv(gzfile("Databases/ComparativeToxicogenomicsDatabase/CTD_genes_diseases.csv.gz"), 
                             header = FALSE, 
                             fill = TRUE, 
                             skip = 29)

colnames(CTD_Gene_Disease) <- c("GeneSymbol", "GeneID", "DiseaseName", 
                                "DiseaseID", "DirectEvidence", "InferenceChemicalName", 
                                "InferenceScore", "OmimIDs", "PubMedIDs")

CTD_Gene_Disease$GeneID <- as.character(CTD_Gene_Disease$GeneID)


# Filter to keep only direct associations
CTD_Gene_Disease <- CTD_Gene_Disease[CTD_Gene_Disease$DirectEvidence %in% c("marker/mechanism", "therapeutic", "marker/mechanism|therapeutic"), ]


# Map entrez gene IDs to ensembl IDs

entrezId_2_ensemblId <- AnnotationDbi::select(org.Hs.eg.db, 
                                              keys = unique(CTD_Gene_Disease$GeneID), 
                                              columns = "ENSEMBL", 
                                              keytype = "ENTREZID")

CTD_Gene_Disease <- CTD_Gene_Disease %>% 
  left_join(entrezId_2_ensemblId, 
            by = c("GeneID" = "ENTREZID"), 
            relationship = "many-to-many") %>% 
  filter(!is.na(ENSEMBL))
rm(entrezId_2_ensemblId)


CTD_Gene_Disease$DiseaseName <- paste0(CTD_Gene_Disease$DiseaseName, " (", CTD_Gene_Disease$DiseaseID, ")")

CTD_Disease2Gene_lib <- split(x = CTD_Gene_Disease$ENSEMBL, f = CTD_Gene_Disease$DiseaseName)
CTD_Disease2Gene_lib <- lapply(CTD_Disease2Gene_lib, unique)

CTD_Disease2Gene_lib <- CTD_Disease2Gene_lib[lengths(CTD_Disease2Gene_lib) >= 5]

if(!dir.exists("IOFiles/Geneset_libraries/")){dir.create("IOFiles/Geneset_libraries/", recursive = TRUE)}
saveRDS(CTD_Disease2Gene_lib, "IOFiles/Geneset_libraries/CTD_Disease2Gene_lib.rds")


#####


# Enrichr 

# Disease gene enrichment library based on Disease_Signatures_from_GEO_2014 (Enrichr)

print("Processing Enrichr")

# Download libraries from Enrichr
if(!dir.exists("Databases/Enrichr")){dir.create("Databases/Enrichr", recursive = TRUE)}
if(!file.exists("Databases/Enrichr/Disease_Signatures_from_GEO_down_2014.txt")){
  download.file(url = "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=Disease_Signatures_from_GEO_down_2014",
                destfile = "Databases/Enrichr/Disease_Signatures_from_GEO_down_2014.txt", method = "wget")
  download.file(url = "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=Disease_Signatures_from_GEO_up_2014",
                destfile = "Databases/Enrichr/Disease_Signatures_from_GEO_up_2014.txt", method = "wget")
}


Enrichr_GeoDiseaseSignatures_Up <- read.table("Databases/Enrichr/Disease_Signatures_from_GEO_up_2014.txt", 
                                              sep = "\t", 
                                              quote = "", 
                                              fill = TRUE)
result <- list()
for(i in 1:nrow(Enrichr_GeoDiseaseSignatures_Up)){
  tmp1 <- Enrichr_GeoDiseaseSignatures_Up[i, -c(1:2)]
  tmp2 <- as.character(apply(tmp1, 1, function(x)gsub(pattern = ",[0-9]+.[0-9]+", replacement = "", x = x)))
  tmp2 <- tmp2[!is.na(tmp2)]
  tmp3 <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, 
                                                 keys = unique(tmp2),
                                                 columns = "ENSEMBL", 
                                                 keytype = "SYMBOL"))
  tmp3 <- tmp3[!is.na(tmp3$ENSEMBL), ]
  result[[i]] <- unique(tmp3$ENSEMBL)
  names(result)[i] <- Enrichr_GeoDiseaseSignatures_Up[i,1]
}
Enrichr_GeoDiseaseSignatures_Up <- lapply(result, unique)
Enrichr_GeoDiseaseSignatures_Up <- Enrichr_GeoDiseaseSignatures_Up[lengths(Enrichr_GeoDiseaseSignatures_Up) >= 5]
names(Enrichr_GeoDiseaseSignatures_Up) <- paste0(names(Enrichr_GeoDiseaseSignatures_Up), " (upreg)")
rm(list = c("tmp1", "tmp2", "tmp3", "result"))


Enrichr_GeoDiseaseSignatures_Down <- read.table("Databases/Enrichr/Disease_Signatures_from_GEO_down_2014.txt", sep = "\t", quote = "", fill = TRUE)
result <- list()
for(i in 1:nrow(Enrichr_GeoDiseaseSignatures_Down)){
  tmp1 <- Enrichr_GeoDiseaseSignatures_Down[i, -c(1:2)]
  tmp2 <- as.character(apply(tmp1, 1, function(x)gsub(pattern = ",[0-9]+.[0-9]+", replacement = "", x = x)))
  tmp2 <- tmp2[!is.na(tmp2)]
  tmp3 <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, 
                                                 keys = unique(tmp2),
                                                 columns = "ENSEMBL", 
                                                 keytype = "SYMBOL"))
  tmp3 <- tmp3[!is.na(tmp3$ENSEMBL), ]
  result[[i]] <- unique(tmp3$ENSEMBL)
  names(result)[i] <- Enrichr_GeoDiseaseSignatures_Down[i,1]
}
Enrichr_GeoDiseaseSignatures_Down <- lapply(result, unique)
Enrichr_GeoDiseaseSignatures_Down <- Enrichr_GeoDiseaseSignatures_Down[lengths(Enrichr_GeoDiseaseSignatures_Down) >= 5]
names(Enrichr_GeoDiseaseSignatures_Down) <- paste0(names(Enrichr_GeoDiseaseSignatures_Down), " (downreg)")
rm(list = c("tmp1", "tmp2", "tmp3", "result"))

Enrichr_Disease2Gene_GeoDiseaseSig_lib <- c(Enrichr_GeoDiseaseSignatures_Up, Enrichr_GeoDiseaseSignatures_Down)

if(!dir.exists("IOFiles/Geneset_libraries/")){dir.create("IOFiles/Geneset_libraries/", recursive = TRUE)}
saveRDS(Enrichr_Disease2Gene_GeoDiseaseSig_lib, "IOFiles/Geneset_libraries/Enrichr_Disease2Gene_GeoDiseaseSig_lib.rds")


#####


# ADReCS 

# Refer for ID info: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4383906/
# The ADR ID can be explained as hierchy. The data contains level 3 and 4 ADR IDs.


print("Processing ADReCS")

if(!dir.exists("Databases/ADReCS/")){dir.create("Databases/ADReCS/", recursive = TRUE)}

if(!file.exists("Databases/ADReCS/ADRAlert2GENE2ID.xlsx")){
  download.file(url = "http://www.bio-add.org/ADReCS-Target/files/download/ADRAlert2GENE2ID.xlsx",
                destfile = "Databases/ADReCS/ADRAlert2GENE2ID.xlsx", method = "wget")
}

if(!file.exists("Databases/ADReCS/P_D_A.xlsx")){
  download.file(url = "https://www.bio-add.org/ADReCS-Target/files/download/P_D_A.xlsx",
                destfile = "Databases/ADReCS/P_D_A.xlsx", method = "wget")
}


# Process the ADR-protein associations
ADReCS_Protein_ADR <- read.xlsx("Databases/ADReCS/P_D_A.xlsx", 
                                check.names = FALSE)
ADReCS_Protein_ADR$ADR.Term <- paste0(ADReCS_Protein_ADR$ADR.Term, " (", ADReCS_Protein_ADR$ADReCS.ID, ")")

uniprotId_2_ensemblId <- AnnotationDbi::select(org.Hs.eg.db, 
                                               keys = unique(ADReCS_Protein_ADR$Uniprot.AC), 
                                               columns = "ENSEMBL", 
                                               keytype = "UNIPROT")
ADReCS_Protein_ADR <- ADReCS_Protein_ADR %>% 
  left_join(uniprotId_2_ensemblId, 
            by = c("Uniprot.AC" = "UNIPROT"), 
            relationship = "many-to-many") %>% 
  filter(!is.na(ENSEMBL))%>%
  rename("ADR.ID" = "ADReCS.ID")
rm(uniprotId_2_ensemblId)


# Process the ADR-gene associations
ADReCS_Gene_ADR <- read.xlsx("Databases/ADReCS/ADRAlert2GENE2ID.xlsx", 
                             check.names = FALSE)
ADReCS_Gene_ADR$ADR.Term <- paste0(ADReCS_Gene_ADR$ADR.Term, " (", ADReCS_Gene_ADR$ADR.ID, ")")
ADReCS_Gene_ADR$GeneID <- as.character(ADReCS_Gene_ADR$GeneID)

entrezId_2_ensemblId <- AnnotationDbi::select(org.Hs.eg.db, 
                                              keys = unique(ADReCS_Gene_ADR$GeneID), 
                                              columns = "ENSEMBL", 
                                              keytype = "ENTREZID")
ADReCS_Gene_ADR <- ADReCS_Gene_ADR %>% 
  left_join(entrezId_2_ensemblId, 
            by = c("GeneID" = "ENTREZID"), 
            relationship = "many-to-many") %>% 
  filter(!is.na(ENSEMBL)) 
rm(entrezId_2_ensemblId)


# Merge the protein and gene level info together
ADReCS_Gene_ADR <- rbind(ADReCS_Protein_ADR[,c("ADR.ID", "ADR.Term", "ENSEMBL")], 
                         ADReCS_Gene_ADR[,c("ADR.ID", "ADR.Term", "ENSEMBL")])

### Level 4 ADR IDs library
ADReCS_Gene_ADR_level4 <- ADReCS_Gene_ADR[lengths(strsplit(ADReCS_Gene_ADR$ADR.ID, split = "\\."))==4,]

ADReCS_ADR2Gene_level4_lib <- split(x = ADReCS_Gene_ADR_level4$ENSEMBL, f = ADReCS_Gene_ADR_level4$ADR.Term)
ADReCS_ADR2Gene_level4_lib <- lapply(ADReCS_ADR2Gene_level4_lib, unique)
ADReCS_ADR2Gene_level4_lib <- ADReCS_ADR2Gene_level4_lib[lengths(ADReCS_ADR2Gene_level4_lib) >= 5]

if(!dir.exists("IOFiles/Geneset_libraries/")){dir.create("IOFiles/Geneset_libraries/", recursive = TRUE)}
saveRDS(ADReCS_ADR2Gene_level4_lib, "IOFiles/Geneset_libraries/ADReCS_ADR2Gene_level4_lib.rds")


### Level 3 ADR IDs library
ADReCS_Gene_ADR_level3 <- ADReCS_Gene_ADR[lengths(strsplit(ADReCS_Gene_ADR$ADR.ID, split = "\\."))==3,]
ADReCS_ADR2Gene_level3_lib <- split(x = ADReCS_Gene_ADR_level3$ENSEMBL, f = ADReCS_Gene_ADR_level3$ADR.Term)
ADReCS_ADR2Gene_level3_lib <- lapply(ADReCS_ADR2Gene_level3_lib, unique)
ADReCS_ADR2Gene_level3_lib <- ADReCS_ADR2Gene_level3_lib[lengths(ADReCS_ADR2Gene_level3_lib) >= 5]

if(!dir.exists("IOFiles/Geneset_libraries/")){dir.create("IOFiles/Geneset_libraries/", recursive = TRUE)}
saveRDS(ADReCS_ADR2Gene_level3_lib, "IOFiles/Geneset_libraries/ADReCS_ADR2Gene_level3_lib.rds")


#####


# DiSignAtlas

print("Processing DiSignAtlas")

if(!dir.exists("Databases/DiSignAtlas/")){dir.create("Databases/DiSignAtlas/", recursive = TRUE)}
if(!file.exists("Databases/DiSignAtlas/Disease_information_DEGs.gmt")){
  download.file(url = "http://www.inbirg.com/disignatlas/download/dis_info_degs",
                destfile = "Databases/DiSignAtlas/Disease_information_DEGs.gmt", method = "wget")
}

if(!file.exists("Databases/DiSignAtlas/Disease_information_Datasets.csv")){
  download.file(url = "http://www.inbirg.com/disignatlas/download/dis_info_datasets",
                destfile = "Databases/DiSignAtlas/Disease_information_Datasets.csv", method = "wget")
}


# Read the disease list
DiSignAtlas_diseases <- read_csv(file = "Databases/DiSignAtlas/Disease_information_Datasets.csv")
DiSignAtlas_diseases <- DiSignAtlas_diseases %>% 
  filter(organism == "Homo sapiens") %>% 
  mutate(lib_name = paste0(disease, " (", diseaseid, ") (", tissue, ", ", library_strategy, ") (", dsaid, ")"))


# Read the DEGs
DiSignAtlas_DEG <- read.table(file = "Databases/DiSignAtlas/Disease_information_DEGs.gmt", 
                              sep = "\t", 
                              fill = TRUE, 
                              header = FALSE, 
                              quote = "")

DiSignAtlas_DEG <- DiSignAtlas_DEG[DiSignAtlas_DEG$V1 %in% DiSignAtlas_diseases$dsaid, ]


result <- list()
for(i in DiSignAtlas_DEG$V1){
  tmp1 <- as.character(DiSignAtlas_DEG[DiSignAtlas_DEG$V1 == i, -c(1:2), drop = TRUE])
  tmp2 <- tmp1[tmp1 != "NA"]
  tmp2 <- tmp2[tmp2 != ""]
  if(length(tmp2) <= 5){ 
    next 
  }else{
    tmp3 <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, 
                                                   keys = unique(tmp2),
                                                   columns = "ENSEMBL", 
                                                   keytype = "ENTREZID", na.action = na.pass))
  }
  
  tmp4 <- DiSignAtlas_diseases[DiSignAtlas_diseases$dsaid %in% DiSignAtlas_DEG[DiSignAtlas_DEG$V1 == i,1], "lib_name", drop = TRUE]
  result[[tmp4]] <- unique(na.exclude(tmp3$ENSEMBL))
  rm(list = c("tmp1", "tmp2", "tmp3", "tmp4"))
}

DiSignAtlas_Disease2Gene_lib <- lapply(result, unique)
DiSignAtlas_Disease2Gene_lib <- DiSignAtlas_Disease2Gene_lib[lengths(DiSignAtlas_Disease2Gene_lib) >= 5]
rm(result)

if(!dir.exists("IOFiles/Geneset_libraries/")){dir.create("IOFiles/Geneset_libraries/", recursive = TRUE)}
saveRDS(DiSignAtlas_Disease2Gene_lib, "IOFiles/Geneset_libraries/DiSignAtlas_Disease2Gene_lib.rds")


#####


# DisGeNet (via OmniPathR)

# Direct download from DisGeNet is not possible.
# Using OmniPathR

DisGeNet_data <- import_omnipath_annotations(resources = "DisGeNet",
                                             wide = TRUE, 
                                             organism = "Homo sapiens")

DisGeNet_data <- DisGeNet_data %>% filter(entity_type == "protein")

# Add ensembl IDs
geneSymbol_2_ensemblId <- AnnotationDbi::select(org.Hs.eg.db, 
                                                keys = unique(DisGeNet_data$genesymbol), 
                                                columns = "ENSEMBL", 
                                                keytype = "SYMBOL")

DisGeNet_data <- DisGeNet_data %>% 
  left_join(geneSymbol_2_ensemblId, 
            by = c("genesymbol" = "SYMBOL"), 
            relationship = "many-to-many") %>% 
  filter(!is.na(ENSEMBL))


# Prepare for diseases
DisGeNet_Disease2Gene_lib <- DisGeNet_data[DisGeNet_data$type == "disease", ]
DisGeNet_Disease2Gene_lib <- split(DisGeNet_Disease2Gene_lib$ENSEMBL, DisGeNet_Disease2Gene_lib$disease)
DisGeNet_Disease2Gene_lib <- lapply(DisGeNet_Disease2Gene_lib, unique)
DisGeNet_Disease2Gene_lib <- DisGeNet_Disease2Gene_lib[lengths(DisGeNet_Disease2Gene_lib) >= 5]

if(!dir.exists("IOFiles/Geneset_libraries/")){dir.create("IOFiles/Geneset_libraries/", recursive = TRUE)}
saveRDS(DisGeNet_Disease2Gene_lib, "IOFiles/Geneset_libraries/DisGeNet_Disease2Gene_lib.rds")


# Prepare for groups
DisGeNet_Group2Gene_lib <- DisGeNet_data[DisGeNet_data$type == "group", ]
DisGeNet_Group2Gene_lib <- split(DisGeNet_Group2Gene_lib$ENSEMBL, DisGeNet_Group2Gene_lib$disease)
DisGeNet_Group2Gene_lib <- lapply(DisGeNet_Group2Gene_lib, unique)
DisGeNet_Group2Gene_lib <- DisGeNet_Group2Gene_lib[lengths(DisGeNet_Group2Gene_lib) >= 5]

if(!dir.exists("IOFiles/Geneset_libraries/")){dir.create("IOFiles/Geneset_libraries/", recursive = TRUE)}
saveRDS(DisGeNet_Group2Gene_lib, "IOFiles/Geneset_libraries/DisGeNet_Group2Gene_lib.rds")


# Prepare for phenotype
DisGeNet_Phenotype2Gene_lib <- DisGeNet_data[DisGeNet_data$type == "phenotype", ]
DisGeNet_Phenotype2Gene_lib <- split(DisGeNet_Phenotype2Gene_lib$ENSEMBL, DisGeNet_Phenotype2Gene_lib$disease)
DisGeNet_Phenotype2Gene_lib <- lapply(DisGeNet_Phenotype2Gene_lib, unique)
DisGeNet_Phenotype2Gene_lib <- DisGeNet_Phenotype2Gene_lib[lengths(DisGeNet_Phenotype2Gene_lib) >= 5]

if(!dir.exists("IOFiles/Geneset_libraries/")){dir.create("IOFiles/Geneset_libraries/", recursive = TRUE)}
saveRDS(DisGeNet_Phenotype2Gene_lib, "IOFiles/Geneset_libraries/DisGeNet_Phenotype2Gene_lib.rds")


#####


# DISEASES

print("Processing DISEASES")

if(!dir.exists("Databases/DISEASES/")){dir.create("Databases/DISEASES/", recursive = TRUE)}
if(!file.exists("Databases/DISEASES/human_disease_textmining_filtered.tsv")){
  download.file(url = "https://download.jensenlab.org/human_disease_textmining_filtered.tsv",
                destfile = "Databases/DISEASES/human_disease_textmining_filtered.tsv", method = "wget")
}
if(!file.exists("Databases/DISEASES/human_disease_knowledge_filtered.tsv")){
  download.file(url = "https://download.jensenlab.org/human_disease_knowledge_filtered.tsv",
                destfile = "Databases/DISEASES/human_disease_knowledge_filtered.tsv", method = "wget")
}
if(!file.exists("Databases/DISEASES/human_disease_experiments_filtered.tsv")){
  download.file(url = "https://download.jensenlab.org/human_disease_experiments_filtered.tsv",
                destfile = "Databases/DISEASES/human_disease_experiments_filtered.tsv", method = "wget")
}


## Disease gene enrichment library based on textmining
DISEASES_Gene_Disease <- read.table("Databases/DISEASES/human_disease_textmining_filtered.tsv", sep = "\t", quote = "", fill = TRUE)
colnames(DISEASES_Gene_Disease) <- c("gene_identifier", "gene_name", "disease_identifier",  "disease_name", "z_score", "confidence_score", "url")
DISEASES_Gene_Disease <- DISEASES_Gene_Disease %>% filter(confidence_score >= 4) %>% select(-c("url"))

# Add ensembl gene IDs
geneSymbol_2_ensemblId <- AnnotationDbi::select(org.Hs.eg.db, 
                                                keys = unique(DISEASES_Gene_Disease$gene_name), 
                                                columns = "ENSEMBL", 
                                                keytype = "SYMBOL")

DISEASES_Gene_Disease <- DISEASES_Gene_Disease %>% 
  left_join(geneSymbol_2_ensemblId, 
            by = c("gene_name" = "SYMBOL"), 
            relationship = "many-to-many") %>% 
  filter(!is.na(ENSEMBL)) %>%
  mutate(disease_name = paste0( disease_name, " (", disease_identifier, ")" ))


DISEASES_Disease2Gene_TM_lib <- split(DISEASES_Gene_Disease$ENSEMBL, f = DISEASES_Gene_Disease$disease_name)
DISEASES_Disease2Gene_TM_lib <- lapply(DISEASES_Disease2Gene_TM_lib, unique)
DISEASES_Disease2Gene_TM_lib <- DISEASES_Disease2Gene_TM_lib[lengths(DISEASES_Disease2Gene_TM_lib) >= 5]

if(!dir.exists("IOFiles/Geneset_libraries/")){dir.create("IOFiles/Geneset_libraries/", recursive = TRUE)}
saveRDS(DISEASES_Disease2Gene_TM_lib, "IOFiles/Geneset_libraries/DISEASES_Disease2Gene_TM_lib.rds")


## Disease gene enrichment library based on text mining
DISEASES_Gene_Disease <- read.table("Databases/DISEASES/human_disease_knowledge_filtered.tsv", sep = "\t", quote = "", fill = TRUE)
colnames(DISEASES_Gene_Disease) <- c("gene_identifier", "gene_name", "disease_identifier",  "disease_name", "source database", "evidence type", "confidence_score")
DISEASES_Gene_Disease <- DISEASES_Gene_Disease %>% filter(confidence_score >= 4) 

# Add ensembl gene IDs
geneSymbol_2_ensemblId <- AnnotationDbi::select(org.Hs.eg.db, 
                                                keys = unique(DISEASES_Gene_Disease$gene_name), 
                                                columns = "ENSEMBL", 
                                                keytype = "SYMBOL")

DISEASES_Gene_Disease <- DISEASES_Gene_Disease %>% 
  left_join(geneSymbol_2_ensemblId, 
            by = c("gene_name" = "SYMBOL"), 
            relationship = "many-to-many") %>% 
  filter(!is.na(ENSEMBL)) %>%
  mutate(disease_name = paste0( disease_name, " (", disease_identifier, ")" ))


DISEASES_Disease2Gene_KNOW_lib <- split(DISEASES_Gene_Disease$ENSEMBL, f = DISEASES_Gene_Disease$disease_name)
DISEASES_Disease2Gene_KNOW_lib <- lapply(DISEASES_Disease2Gene_KNOW_lib, unique)
DISEASES_Disease2Gene_KNOW_lib <- DISEASES_Disease2Gene_KNOW_lib[lengths(DISEASES_Disease2Gene_KNOW_lib) >= 5]

if(!dir.exists("IOFiles/Geneset_libraries/")){dir.create("IOFiles/Geneset_libraries/", recursive = TRUE)}
saveRDS(DISEASES_Disease2Gene_KNOW_lib, "IOFiles/Geneset_libraries/DISEASES_Disease2Gene_KNOW_lib.rds")


# ## Disease gene enrichment library based on experiments
# DISEASES_Gene_Disease <- read.table("Databases/DISEASES/human_disease_experiments_filtered.tsv", sep = "\t", quote = "", fill = TRUE)
# colnames(DISEASES_Gene_Disease) <- c("gene_identifier", "gene_name", "disease_identifier",  "disease_name", "source_database", " source_score", "confidence_score")
# DISEASES_Gene_Disease <- DISEASES_Gene_Disease %>% filter(confidence_score >= 4) 
# 
# # Add ensembl gene IDs
# geneSymbol_2_ensemblId <- AnnotationDbi::select(org.Hs.eg.db, 
#                                                 keys = unique(DISEASES_Gene_Disease$gene_name), 
#                                                 columns = "ENSEMBL", 
#                                                 keytype = "SYMBOL")
# 
# DISEASES_Gene_Disease <- DISEASES_Gene_Disease %>% 
#   left_join(geneSymbol_2_ensemblId, 
#             by = c("gene_name" = "SYMBOL"), 
#             relationship = "many-to-many") %>% 
#   filter(!is.na(ENSEMBL)) %>%
#   mutate(disease_name = paste0( disease_name, " (", disease_identifier, ")" ))
# 
# DISEASES_Disease2Gene_EXP_lib <- split(DISEASES_Gene_Disease$ENSEMBL, f = DISEASES_Gene_Disease$disease_name)
# DISEASES_Disease2Gene_EXP_lib <- lapply(DISEASES_Disease2Gene_EXP_lib, unique)
# DISEASES_Disease2Gene_EXP_lib <- DISEASES_Disease2Gene_EXP_lib[lengths(DISEASES_Disease2Gene_EXP_lib) >= 5]
# 
# if(!dir.exists("IOFiles/Geneset_libraries/")){dir.create("IOFiles/Geneset_libraries/", recursive = TRUE)}
# saveRDS(DISEASES_Disease2Gene_EXP_lib, "IOFiles/Geneset_libraries/DISEASES_Disease2Gene_EXP_lib.rds") 


#####


print(warnings())