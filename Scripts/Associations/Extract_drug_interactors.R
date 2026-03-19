set.seed(5081)


# Script to extract known interactors of drugs from DrugBank


# Load libraries
library(unixtools)
library(org.Hs.eg.db)
library(tidyverse)


# Set temporary directory
if(!dir.exists("tmp_dir/")){dir.create("tmp_dir/", recursive = TRUE)}
set.tempdir("tmp_dir/")


#####


# Download the complete DrugBank database as XML and parse it using dbparser

if(!dir.exists("Databases/DrugBank")){dir.create("Databases/DrugBank/", recursive = TRUE)} 
if(!file.exists("Databases/DrugBank/drugbank_all_full_database.xml.zip")){ 
  warning(paste0("ERROR: DrugBank database file not found !!! \n", 
                 "Download file in terminal using:\n", 
                 "\t curl -Lfv -o filename.zip -u EMAIL:PASSWORD https://go.drugbank.com/releases/5-1-13/downloads/all-full-database"))
  stop()
}


if(!file.exists("Databases/DrugBank/parsed_DrugBank_data.rds")){
  require(dbparser)
  dvobj <- parseDrugBank(db_path = "Databases/DrugBank/drugbank_all_full_database.xml.zip",
                         drug_options = drug_node_options(),
                         parse_salts = TRUE,
                         parse_products = TRUE,
                         references_options = references_node_options(),
                         cett_options = cett_nodes_options())
  saveRDS(dvobj, "Databases/DrugBank/parsed_DrugBank_data.rds")
}


######


# Extract all possible interactors

DrugBank_data <- readRDS("Databases/DrugBank/parsed_DrugBank_data.rds")

DrugBank_drugInteractors_idMap <- list()
for(type in c("carriers", "enzymes", "targets", "transporters")){
  tmp1 <- DrugBank_data[["cett"]][[type]][["polypeptides"]][["external_identy"]] 
  colnames(tmp1)[3] <- "interactor_id"
  DrugBank_drugInteractors_idMap[[type]] <- tmp1
  rm(tmp1)
}

DrugBank_drugInteractors_idMap <- DrugBank_drugInteractors_idMap %>% 
  bind_rows() %>% 
  distinct() %>% 
  pivot_wider(names_from = "resource", 
              values_from = "identifier")

# Filter to keep only human protein targets
DrugBank_drugInteractors_idMap <- DrugBank_drugInteractors_idMap[grep("_HUMAN$", 
                                                                      x = DrugBank_drugInteractors_idMap$`UniProt Accession`),]


# Add Ensembl Gene ID
uniprotID_2_ensemblID <- AnnotationDbi::select(org.Hs.eg.db, 
                                               keys = unique(DrugBank_drugInteractors_idMap$UniProtKB), 
                                               columns = c("UNIPROT", "ENSEMBL"), 
                                               keytype = "UNIPROT")


DrugBank_drugInteractors_idMap <- DrugBank_drugInteractors_idMap %>% 
  left_join(uniprotID_2_ensemblID, 
            by = c("UniProtKB" = "UNIPROT"), 
            relationship = "many-to-many")


#####


# Extract drug-interactor relationship
DrugBank_drugInteractors <- list()
for(type in c("carriers", "enzymes", "targets", "transporters")){
  tmp1 <- DrugBank_data[["cett"]][[type]][["general_information"]] 
  colnames(tmp1)[1] <- "interactor_id"
  DrugBank_drugInteractors[[type]] <- tmp1 %>% select(!c("position"))
  rm(tmp1)
}
DrugBank_drugInteractors <- DrugBank_drugInteractors %>% 
  bind_rows(.id = "type") %>% 
  filter(organism %in% "Humans")

# Add Ensembl Gene Id
DrugBank_drugInteractors <- DrugBank_drugInteractors %>% 
  left_join(DrugBank_drugInteractors_idMap %>% select(c("interactor_id", "ENSEMBL", "UniProtKB")),
            by = "interactor_id", 
            relationship = "many-to-many") %>% 
  dplyr::rename("ensembl_gene_id" = "ENSEMBL", "UniProtKB_protein_accession" = "UniProtKB") %>%
  filter( !(is.na(ensembl_gene_id) | is.na(UniProtKB_protein_accession)) ) # Filters non-protein targets also


#####


# Save the networks

if(!dir.exists("IOFiles/Associations/")){dir.create("IOFiles/Associations/", recursive = TRUE)}

saveRDS(DrugBank_drugInteractors, file = "IOFiles/Associations/DrugBank_drugInteractors.rds")


#####


print(warnings())