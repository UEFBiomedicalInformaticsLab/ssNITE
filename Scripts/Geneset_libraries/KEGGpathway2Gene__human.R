set.seed(5081)


# Script to retrieve human KEGG pathway genes


# Load libraries
library(unixtools)
library(KEGGREST)
library(org.Hs.eg.db)
library(tidyverse)

# Set temporary directory
if(!dir.exists("tmp_dir/")){dir.create("tmp_dir/", recursive = TRUE)}
set.tempdir("tmp_dir/")

if(!dir.exists("IOFiles/Geneset_libraries/")){ dir.create("IOFiles/Geneset_libraries/", recursive = TRUE) }


#####


# Download the list of KEGG pathways
KEGG_pathwayList <- keggList(database = "pathway", organism = "hsa")
KEGG_pathwayList <- tibble("pathway_id" = names(KEGG_pathwayList), "pathway_name" = KEGG_pathwayList)


# Download all genes liked to KEGG pathways
KEGG_pathway_gene <- keggLink(target = "hsa", source = "pathway") 
KEGG_pathway_gene <- tibble("pathway_id" = names(KEGG_pathway_gene), 
                            "KEGG_gene_id" = KEGG_pathway_gene)
KEGG_pathway_gene$pathway_id <- gsub("path:", "", KEGG_pathway_gene$pathway_id)


# Convert KEGG gene ID to Ensembl gene ID
KEGG_geneList <- keggConv(target = "ncbi-geneid", 
                          source = unique(KEGG_pathway_gene$KEGG_gene_id)) %>% 
  as.data.frame() %>%
  rownames_to_column("KEGG_gene_id") %>% 
  dplyr::rename("entrez_gene_ID" = ".") 

KEGG_geneList$entrez_gene_ID <- gsub("ncbi-geneid:", "", KEGG_geneList$entrez_gene_ID)

entrezGeneID_2_ensemblGeneID <- AnnotationDbi::select(org.Hs.eg.db,
                                                      keys = KEGG_geneList$entrez_gene_ID,
                                                      columns = c("ENTREZID", "ENSEMBL", "SYMBOL"),
                                                      keytype = "ENTREZID")

KEGG_geneList <- left_join(KEGG_geneList, 
                           entrezGeneID_2_ensemblGeneID,
                           by = c("entrez_gene_ID" = "ENTREZID"))


# Merge all data
KEGG_pathway_gene <- KEGG_pathway_gene %>% 
  left_join(KEGG_pathwayList, by = "pathway_id") %>% 
  left_join(KEGG_geneList, by = "KEGG_gene_id", relationship = "many-to-many")
KEGG_pathway_gene$pathway_name <- paste0(KEGG_pathway_gene$pathway_name, " (", KEGG_pathway_gene$pathway_id, ")")
KEGG_pathway_gene <- KEGG_pathway_gene %>% filter(!is.na(ENSEMBL))

# Create lists
KEGG_pathway_gene <- split(KEGG_pathway_gene$ENSEMBL, KEGG_pathway_gene$pathway_name)
KEGG_pathway_gene <- lapply(KEGG_pathway_gene, unique)

saveRDS(KEGG_pathway_gene, "IOFiles/Geneset_libraries/KEGG_Pathway2Gene_human_lib.rds")


#####


print(warnings())