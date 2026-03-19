# Functions for ID conversions


#Load libraries
require(tidyverse)


#####


# Function to convert list object to dataframe
func_list_2_df <- function(query_list){
  if(is.list(query_list)){
    df <- as.data.frame(matrix(nrow = length(names(query_list)), ncol = 2))
    for(i in 1:length(names(query_list))){
      df[i,1] <- names(query_list[i])
      df[i,2] <- paste(query_list[[i]], collapse = ";")
    }
    df <- as.data.frame(separate_rows(df, V2, sep = ";"))
    return(df)
  } else {
    print("Error: Not a list")
  }
}


#####


# Function cbind.fill
cbind.fill <- function(...){
  nm <- list(...) 
  nm<-lapply(nm, as.matrix)
  n <- max(sapply(nm, nrow)) 
  do.call(cbind, lapply(nm, function (x) 
    rbind(x, matrix(, n-nrow(x), ncol(x))))) 
}


#####


# Function to convert PharmGKB IDs
# Download required files for mapping

if(!dir.exists("Databases/PharmGKB/")){dir.create("Databases/PharmGKB/", recursive = TRUE)}

if(!file.exists("Databases/PharmGKB/genes.tsv")){
  download.file(url = "https://api.pharmgkb.org/v1/download/file/data/genes.zip",
                destfile = "Databases/PharmGKB/genes.zip", method = "wget")
  unzip("Databases/PharmGKB/genes.zip", files = "genes.tsv", exdir = "Databases/PharmGKB/")
}

if(!file.exists("Databases/PharmGKB/phenotypes.tsv")){
  download.file(url = "https://api.pharmgkb.org/v1/download/file/data/phenotypes.zip",
                destfile = "Databases/PharmGKB/phenotypes.zip", method = "wget")
  unzip("Databases/PharmGKB/phenotypes.zip", files = "phenotypes.tsv", exdir = "Databases/PharmGKB/")
}

if(!file.exists("Databases/PharmGKB/chemicals.tsv")){
  download.file(url = "https://api.pharmgkb.org/v1/download/file/data/chemicals.zip",
                destfile = "Databases/PharmGKB/chemicals.zip", method = "wget")
  unzip("Databases/PharmGKB/chemicals.zip", files = "chemicals.tsv", exdir = "Databases/PharmGKB/")
}


#####


# Funtion to map PharmGKB IDs to other type of IDs
func_pharmgkb_mapping <- function(query_list, map_to = c("entrez_gene_id", "cas_number")){
  switch (map_to,
          entrez_gene_id = {
            PharmGKB_genes <- read.table("Databases/PharmGKB/genes.tsv", sep = "\t", fill = TRUE, quote = "", header = TRUE)
            pharmgkbId_2_entrezGeneId <- PharmGKB_genes[PharmGKB_genes$PharmGKB.Accession.Id %in% query_list, c("PharmGKB.Accession.Id", "NCBI.Gene.ID")]
            colnames(pharmgkbId_2_entrezGeneId) <- c("Entity_id", "entrez_gene_id")
            return(pharmgkbId_2_entrezGeneId)
          },
          ensembl_gene_id = {
            PharmGKB_genes <- read.table("Databases/PharmGKB/genes.tsv", sep = "\t", fill = TRUE, quote = "", header = TRUE)
            pharmgkbId_2_ensemblGeneId <- PharmGKB_genes[PharmGKB_genes$PharmGKB.Accession.Id %in% query_list, c("PharmGKB.Accession.Id", "Ensembl.Id")]
            colnames(pharmgkbId_2_ensemblGeneId) <- c("Entity_id", "ensembl_gene_id")
            pharmgkbId_2_ensemblGeneId <- as.data.frame(separate_rows(pharmgkbId_2_ensemblGeneId, "ensembl_gene_id", sep = ", ", convert = TRUE))
            pharmgkbId_2_ensemblGeneId$ensembl_gene_id <- gsub(pattern = '^\\"|\\"$', replacement = "", x = pharmgkbId_2_ensemblGeneId$ensembl_gene_id)
            return(pharmgkbId_2_ensemblGeneId)
          },
          cas_number = {
            PharmGKB_chemicals <- read.table("Databases/PharmGKB/chemicals.tsv", sep = "\t", fill = TRUE, quote = "", header = TRUE)
            PharmGKB_chemicals <- separate_rows(PharmGKB_chemicals, c("Cross.references"), sep = ",")
            PharmGKB_chemicals$Cross.references <- gsub(pattern = '\"', replacement = "", x = PharmGKB_chemicals$Cross.references)
            PharmGKB_chemicals <- as.data.frame(separate(PharmGKB_chemicals, c("Cross.references"), into = c("Source", "ID"), sep = ":", extra = "merge", fill = "right"))
            
            require(metaboliteIDmapping)
            PharmGKB_chemicals$cas <- NA
            PharmGKB_chemicals[PharmGKB_chemicals$Source == "PubChem Compound",]$cas <- metabolitesMapping$CAS[match(PharmGKB_chemicals[PharmGKB_chemicals$Source == "PubChem Compound",]$ID, metabolitesMapping$CID)]
            PharmGKB_chemicals[PharmGKB_chemicals$Source == "HMDB",]$cas <- metabolitesMapping$CAS[match(PharmGKB_chemicals[PharmGKB_chemicals$Source == "HMDB",]$ID, metabolitesMapping$HMDB)]
            PharmGKB_chemicals[PharmGKB_chemicals$Source == "KEGG Compound",]$cas <- metabolitesMapping$CAS[match(PharmGKB_chemicals[PharmGKB_chemicals$Source == "KEGG Compound",]$ID, metabolitesMapping$KEGG)]
            PharmGKB_chemicals[PharmGKB_chemicals$Source == "DrugBank",]$cas <- metabolitesMapping$CAS[match(PharmGKB_chemicals[PharmGKB_chemicals$Source == "DrugBank",]$ID, metabolitesMapping$Drugbank)]
            PharmGKB_chemicals[PharmGKB_chemicals$Source == "ChEBI",]$ID <- gsub(pattern = "^CHEBI:", replacement = "", x = PharmGKB_chemicals[PharmGKB_chemicals$Source == "ChEBI",]$ID)
            PharmGKB_chemicals[PharmGKB_chemicals$Source == "ChEBI",]$cas <- metabolitesMapping$CAS[match(PharmGKB_chemicals[PharmGKB_chemicals$Source == "ChEBI",]$ID, metabolitesMapping$ChEBI)]
            pharmgkbId_2_casNumber <- PharmGKB_chemicals[PharmGKB_chemicals$PharmGKB.Accession.Id %in% query_list, c("PharmGKB.Accession.Id", "cas")]
            colnames(pharmgkbId_2_casNumber) <- c("Entity_id", "cas")
            return(pharmgkbId_2_casNumber)
          }
  )
  
}
