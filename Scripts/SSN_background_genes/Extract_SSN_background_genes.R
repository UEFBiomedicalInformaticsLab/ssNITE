set.seed(5081)


# Script to short-list the genes to be used as background in creating SSN (Cancer)


# Load libraries
library(unixtools)
library(KEGGREST)
library(org.Hs.eg.db)
library(igraph)
library(UpSetR)
library(tidyverse)


# Set temporary directory
if(!dir.exists("tmp_dir/")){dir.create("tmp_dir/", recursive = TRUE)}
set.tempdir("tmp_dir/")

if(!dir.exists("IOFiles/SSN_background_genes/")){dir.create("IOFiles/SSN_background_genes/", recursive = TRUE)}


#####


# Retrieve DrugBank data
DrugBank_data <- readRDS("Databases/DrugBank/parsed_DrugBank_data.rds")
DrugBank_drugInteractors <- readRDS("IOFiles/Associations/DrugBank_drugInteractors.rds")


#####


# Add known interactors of licensed cancer drugs as background genes

# Read the drug type information
DrugBank_drug_type <- DrugBank_data$drugs$general_information
DrugBank_drug_type <- DrugBank_drug_type[DrugBank_drug_type$type == "small molecule", ] # retain only small molecular drugs


# Download list of licensed anti-cancer drugs
if(!dir.exists("Databases/Cancer_Drug_Database/")){dir.create("Databases/Cancer_Drug_Database/", recursive = TRUE)}
if(!file.exists("Databases/Cancer_Drug_Database/cancerdrugsdb.txt")){
  download.file(url = "https://sciencedata.anticancerfund.org/pages//cancerdrugsdb.txt",
                destfile = "Databases/Cancer_Drug_Database/cancerdrugsdb.txt", method = "wget")
}


CancerDrugDb_drugs <- read.table("Databases/Cancer_Drug_Database/cancerdrugsdb.txt", sep = "\t", 
                                 fill = TRUE, header = TRUE, strip.white = TRUE, check.names = FALSE)
CancerDrugDb_drugs <- CancerDrugDb_drugs[, c("Product", "DrugBank ID", "Indications")]
CancerDrugDb_drugs <- CancerDrugDb_drugs %>% mutate(DrugBank_drug_id = gsub("<.*>(DB\\d{5})</.*", "\\1", `DrugBank ID`))

# Filter to keep only small molecular drugs
CancerDrugDb_drugs <- CancerDrugDb_drugs %>% filter(DrugBank_drug_id %in% DrugBank_drug_type$drugbank_id)

# Get the targets of the drugs
CancerDrugDb_drugs <- DrugBank_drugInteractors %>% filter(drugbank_id %in% CancerDrugDb_drugs$DrugBank_drug_id)


#####


# Add known interactors of antineoplastic drugs based on ATC code (L01) as background genes

DrugBank_atcL01_drugs <- DrugBank_data$drugs$atc_codes
DrugBank_atcL01_drugs <- DrugBank_atcL01_drugs %>%
  filter(code_3 == "L01") %>% 
  filter(drugbank_id %in% DrugBank_drug_type$drugbank_id)
DrugBank_atcL01_drugs <- DrugBank_drugInteractors %>%
  filter(drugbank_id %in% DrugBank_atcL01_drugs$drugbank_id)


#####


# Add genes in SMPDB drug action/metabolism pathways for anti-neoplastic drugs and licensed anti-cancer drugs as background genes

# Get the list of drug action and drug metabolism pathways linked to anti-cancer and anti-neoplastic drugs
DrugBank_drugPathways <- DrugBank_data$drugs$pathway$general_information
DrugBank_drugPathways <- DrugBank_drugPathways %>% 
  filter( category %in% c("drug_action", "drug_metabolism") ) %>% 
  filter( drugbank_id %in% unique(c(CancerDrugDb_drugs$drugbank_id, DrugBank_atcL01_drugs$drugbank_id)) ) 

# Download the SMPDB pathway linked protein info
if(!dir.exists("Databases/SMPDB")){dir.create("Databases/SMPDB", recursive = TRUE)}
if(!file.exists("Databases/SMPDB/smpdb_proteins.csv.zip")){
  download.file(url = "https://smpdb.ca/downloads/smpdb_proteins.csv.zip",
                destfile = "Databases/SMPDB/smpdb_proteins.csv.zip", method = "wget")
  unzip(zipfile = "Databases/SMPDB/smpdb_proteins.csv.zip", exdir = "Databases/SMPDB/smpdb_proteins")
}


# Get the list of proteins linked to the identified pathways
SMPDB_pathways <- lapply(unique(DrugBank_drugPathways$smpdb_id), function(x){
  
  file <- paste0("Databases/SMPDB/smpdb_proteins/", x, "_proteins.csv")
  if(file.exists(file)){
    read.csv(file, header = TRUE, check.names = FALSE)
  }
  
})

SMPDB_pathways <- data.table::rbindlist(SMPDB_pathways)

# Add ensembl IDs
uniprotProteinID_2_ensemblGeneID <- AnnotationDbi::select(org.Hs.eg.db,
                                                          keys = unique(SMPDB_pathways$`Uniprot ID`),
                                                          columns = c("UNIPROT", "ENSEMBL"),
                                                          keytype = "UNIPROT")

SMPDB_pathways <- SMPDB_pathways %>% 
  left_join(uniprotProteinID_2_ensemblGeneID, 
            by = c("Uniprot ID" = "UNIPROT"), 
            relationship = "many-to-many") %>% 
  filter(!is.na(ENSEMBL)) %>% 
  rename("Ensembl Gene ID" = "ENSEMBL") 


#####


# Read KEGG pathway gene list
KEGG_Pathway2Gene <- readRDS("IOFiles/Geneset_libraries/KEGG_Pathway2Gene_human_lib.rds")

# Add genes in Cancer pathway (hsa05200) in KEGG as background genes
KEGG_geneList <- KEGG_Pathway2Gene[grep(pattern = "hsa05200", x = names(KEGG_Pathway2Gene))]
KEGG_geneList <- sort(unique(unlist(KEGG_geneList, use.names = FALSE)))


#####


# Add genes in KEGG pathways associated with hallmarkers of cancer


# Download the HOC associated pathways
if(!dir.exists("Databases/Database_of_Cancer_Hallmark_Genes/")){dir.create("Databases/Database_of_Cancer_Hallmark_Genes/", recursive = TRUE)}
if(!file.exists("Databases/Database_of_Cancer_Hallmark_Genes/Table_1.xls")){
  download.file(url = "https://www.frontiersin.org/api/v4/articles/482153/file/Table_1.XLS/482153_supplementary-materials_tables_1_xls/1",
                destfile = "Databases/Database_of_Cancer_Hallmark_Genes/Table_1.xls", method = "wget")
}


HOC_pathways <- readxl::read_excel("Databases/Database_of_Cancer_Hallmark_Genes/Table_1.xls", skip = 2)
HOC_pathways <- HOC_pathways$Pathway
HOC_pathways <- sort(unique(unlist(str_split(HOC_pathways, "\n"))))
HOC_pathways <- HOC_pathways[HOC_pathways != ""]

HOC_geneList <- KEGG_Pathway2Gene[grep(pattern = paste(HOC_pathways, collapse = "|"), x = names(KEGG_Pathway2Gene))]
HOC_geneList <- sort(unique(unlist(HOC_geneList, use.names = FALSE)))
rm(HOC_pathways)


#####


# Add genes reported as pan cancer oncogenic genes in OncoVar as background genes

# Download list of licensed anti-cancer drugs
if(!dir.exists("Databases/OncoVar/")){dir.create("Databases/OncoVar/", recursive = TRUE)}
if(!file.exists("Databases/OncoVar/TCGA.PanCancer.onco.genes.OncoVar.tsv.gz")){
  download.file(url = "https://oncovar.org/resource/download/Onco_genes_OncoVar_TCGA/TCGA.PanCancer.onco.genes.OncoVar.tsv.gz",
                destfile = "Databases/OncoVar/TCGA.PanCancer.onco.genes.OncoVar.tsv.gz", method = "wget")
}
if(!file.exists("Databases/OncoVar/ICGC.PanCancer.onco.genes.OncoVar.tsv.gz")){
  download.file(url = "https://oncovar.org/resource/download/Onco_genes_OncoVar_ICGC/ICGC.PanCancer.onco.genes.OncoVar.tsv.gz",
                destfile = "Databases/OncoVar/ICGC.PanCancer.onco.genes.OncoVar.tsv.gz", method = "wget")
}

OncoVar_oncoGenes_tcga <- read.table("Databases/OncoVar/TCGA.PanCancer.onco.genes.OncoVar.tsv.gz", header = TRUE, row.names = NULL, sep = "\t")
OncoVar_oncoGenes_tcga <- OncoVar_oncoGenes_tcga %>% filter( (Gene_biotype == "protein_coding") )

OncoVar_oncoGenes_icgc <- read.table("Databases/OncoVar/ICGC.PanCancer.onco.genes.OncoVar.tsv.gz", header = TRUE, row.names = NULL, sep = "\t")
OncoVar_oncoGenes_icgc <- OncoVar_oncoGenes_icgc %>% filter( (Gene_biotype == "protein_coding") )


#####


# Add genes listed as oncogenes or tumor suppresor genes in Network of Cancer Genes database as background genes

if(!dir.exists("Databases/Network_of_Cancer_Genes/")){dir.create("Databases/Network_of_Cancer_Genes/", recursive = TRUE)}

if(!file.exists("Databases/Network_of_Cancer_Genes/NCG_cancerdrivers_annotation_supporting_evidence.tsv")){
  stop("Manually download the file 'NCG_cancerdrivers_annotation_supporting_evidence.tsv' (list of all cancer drivers and their annotation and supporting evidence) from http://network-cancer-genes.org/download.php")
}


NCG_genes <- read.table("Databases/Network_of_Cancer_Genes/NCG_cancerdrivers_annotation_supporting_evidence.tsv", header = TRUE, sep = "\t")

NCG_genes <- NCG_genes %>% 
  filter(coding_status %in% "coding") %>% 
  filter((NCG_oncogene == 1) | (NCG_tsg == 1))

NCG_genes$entrez <- as.character(NCG_genes$entrez)


NCG_genes <- NCG_genes %>%
  group_by(entrez, symbol) %>%
  summarise(across(c("pubmed_id", "type", "organ_system", "primary_site", 
                     "cancer_type", "method", "coding_status", "cgc_annotation", 
                     "vogelstein_annotation", "saito_annotation", "NCG_oncogene", "NCG_tsg"), 
                   ~ paste(sort(unique(.)), collapse = "; ")), .groups = "drop")


entrezGeneID_2_ensemblGeneID <- AnnotationDbi::select(org.Hs.eg.db,
                                                      keys = unique(NCG_genes$entrez),
                                                      columns = c("ENTREZID", "ENSEMBL"),
                                                      keytype = "ENTREZID")

NCG_genes <- NCG_genes %>% 
  left_join(entrezGeneID_2_ensemblGeneID, 
            by = c("entrez" = "ENTREZID")) 


#####


# Add key resistance molecules from DRESIS as background genes

if(!dir.exists("Databases/DRESIS/")){dir.create("Databases/DRESIS/", recursive = TRUE)}
if(!file.exists("Databases/DRESIS/3-1. The general information of disease related with resistance.txt")){
  download.file(url = "https://dresis.idrblab.net/sites/files/full_data/3-1.%20The%20general%20information%20of%20disease%20related%20with%20resistance.txt",
                destfile = "Databases/DRESIS/3-1. The general information of disease related with resistance.txt", method = "wget")
}
if(!file.exists("Databases/DRESIS/3-4. The panorama of resistant drugs for particular disease.txt")){
  download.file(url = "https://dresis.idrblab.net/sites/files/full_data/3-4.%20The%20panorama%20of%20resistant%20drugs%20for%20particular%20disease.txt",
                destfile = "Databases/DRESIS/3-4. The panorama of resistant drugs for particular disease.txt", method = "wget")
}
if(!file.exists("Databases/DRESIS/4-1. The general information of molecular associated with resistance.txt")){
  download.file(url = "https://dresis.idrblab.net/sites/files/full_data/4-1.%20The%20general%20information%20of%20molecular%20associated%20with%20resistance.txt",
                destfile = "Databases/DRESIS/4-1. The general information of molecular associated with resistance.txt", method = "wget")
}


# Identify all the cancers 
DRESIS_diseases <- read.table("Databases/DRESIS/3-1. The general information of disease related with resistance.txt", header = TRUE, sep = "\t")
DRESIS_diseases <- DRESIS_diseases %>% 
  filter(str_detect(Disease_ICD, "^ICD-11: 2")) # ICD-11:02 is for neoplasms


# Get molecules listes as key for resistance in the selected cancers
DRESIS_links <- read.table("Databases/DRESIS/3-4. The panorama of resistant drugs for particular disease.txt", header = TRUE, sep = "\t")
DRESIS_links <- DRESIS_links %>% filter(Disease_ID %in% DRESIS_diseases$Disease_ID)


# Get details about the molecules
DRESIS_molecules <- read.table("Databases/DRESIS/4-1. The general information of molecular associated with resistance.txt", 
                               header = TRUE, sep = "\t", row.names = NULL, fill = TRUE, quote = "")
DRESIS_molecules <- DRESIS_molecules %>% 
  filter(Molecule_ID %in% DRESIS_links$Molecule_ID) %>% 
  filter((Molecule_type %in% "Protein") & (Molecule_species %in% "Homo sapiens")) %>% 
  filter(str_detect(Ensembl.ID, "^ENSG"))


#####


# Combine the gene list from disease and drugs
network_background_genes <- list("CancerDrugDb" = CancerDrugDb_drugs$ensembl_gene_id, 
                                 "DrugBank_atcL01" = DrugBank_atcL01_drugs$ensembl_gene_id, 
                                 "SMPDB" = SMPDB_pathways$`Ensembl Gene ID`,
                                 "KEGG" = KEGG_geneList, 
                                 "HOC_paths" = HOC_geneList, 
                                 "OncoVar_tcga" = OncoVar_oncoGenes_tcga$Ensembl_ID, 
                                 "OncoVar_icgc" = OncoVar_oncoGenes_icgc$Ensembl_ID, 
                                 "NCG" = NCG_genes$ENSEMBL,
                                 "DRESIS" = DRESIS_molecules$Ensembl.ID)  


network_background_genes <- lapply(network_background_genes, function(x){
  
  tmp1 <- AnnotationDbi::select(org.Hs.eg.db, 
                                keys = unique(x), 
                                columns = c("ENSEMBL", "SYMBOL", "ENTREZID", "GENETYPE"), 
                                keytype = "ENSEMBL")
  colnames(tmp1) <- c("ensembl_gene_id", "gene_symbol", "entrez_gene_id", "gene_type")
  tmp1 %>% filter(gene_type == "protein-coding")
  
})


# Plot the overlap between gene sets

plot_data <- lapply(network_background_genes, function(x) unique(x$ensembl_gene_id))

tiff(paste0("IOFiles/SSN_background_genes/SSN_backgroundGenes_overlap.tiff"),
     width = 40, 
     height = 18,
     units = "cm", compression = "lzw", res = 1200)

upset(data = fromList(plot_data), 
      nsets = length(plot_data), 
      nintersects = NA, 
      keep.order = FALSE, 
      order.by = "freq", 
      text.scale = c(1, 1, 1, 1, 1, 0.75), 
      number.angles = 0)

dev.off()


# Aggregate and export
network_background_genes <- network_background_genes %>% 
  bind_rows(.id = "source") %>% 
  group_by(ensembl_gene_id, gene_symbol, entrez_gene_id) %>% 
  summarise("source" = paste(sort(unique(source)), collapse = ";"))

write.csv(network_background_genes, "IOFiles/SSN_background_genes/SSN_background_genes.csv", row.names = FALSE)


#####


# Read the STRING PPI
STRING_Net <- readRDS("IOFiles/Networks/STRING_PPI_Net.rds")


# Extract the subnetworks of the background genes from the STRING PPI
subnet <- list()

for(type in c("neighborhood", "fusion", "cooccurence", "coexpression", "experimental", "database", "textmining", "combined_score")){
  
  tmp1 <- switch( type, 
                  "neighborhood" = {delete_edges(graph = STRING_Net, edges = E(STRING_Net)[E(STRING_Net)$neighborhood <= 0])},
                  "fusion" = {delete_edges(graph = STRING_Net, edges = E(STRING_Net)[E(STRING_Net)$fusion <= 0])}, 
                  "cooccurence" = {delete_edges(graph = STRING_Net, edges = E(STRING_Net)[E(STRING_Net)$cooccurence <= 0])}, 
                  "coexpression" = {delete_edges(graph = STRING_Net, edges = E(STRING_Net)[E(STRING_Net)$coexpression <= 0])}, 
                  "experimental" = {delete_edges(graph = STRING_Net, edges = E(STRING_Net)[E(STRING_Net)$experimental <= 0])}, 
                  "database" = {delete_edges(graph = STRING_Net, edges = E(STRING_Net)[E(STRING_Net)$database <= 0])}, 
                  "textmining" = {delete_edges(graph = STRING_Net, edges = E(STRING_Net)[E(STRING_Net)$textmining <= 0])}, 
                  "combined_score" = {delete_edges(graph = STRING_Net, edges = E(STRING_Net)[E(STRING_Net)$combined_score <= 0])} )
  
  cat(paste0("\nEdge type : ", type))
  cat(paste0("\n --Nodes (ref) : ", vcount(tmp1)))
  cat(paste0("\n --Edges (ref) : ", ecount(tmp1)))
  cat(paste0("\n --Components (ref) : ", count_components(tmp1)))
  
  subnet[[type]] <- induced_subgraph(graph = tmp1, vids = V(tmp1)[V(tmp1)$name %in% unique(network_background_genes$ensembl_gene_id)])
  
  cat(paste0("\n --Nodes (subnet) : ", vcount(subnet[[type]])))
  cat(paste0("\n --Edges (subnet) : ", ecount(subnet[[type]])))
  cat(paste0("\n --Components (subnet) : ", count_components(subnet[[type]]), "\n"))
  
  # RCy3::createNetworkFromIgraph(subnet[[type]], title = type)
  
}
rm(tmp1)

saveRDS(subnet, "IOFiles/SSN_background_genes/SSN_backgroundGenes_PPI_subnetworks.rds")


#####


print(warnings())