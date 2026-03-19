set.seed(5081)


# Script to build protein-protein interaction network from STRING database


# Load libraries
library(unixtools)
library(biomaRt)
library(igraph)
library(tidyverse)


# Set temporary directory
if(!dir.exists("tmp_dir/")){dir.create("tmp_dir/", recursive = TRUE)}
set.tempdir("tmp_dir/")


#####


# Download the protein-protein interactions from stringdb
if(!dir.exists("Databases/StringDB/")){dir.create("Databases/StringDB", recursive = TRUE)}

if(!file.exists("Databases/StringDB/9606.protein.links.detailed.v12.0.onlyAB.txt.gz")){
  download.file(url = "https://stringdb-downloads.org/download/stream/protein.links.detailed.v12.0/9606.protein.links.detailed.v12.0.onlyAB.txt.gz",
                destfile = "Databases/StringDB/9606.protein.links.detailed.v12.0.onlyAB.txt.gz", method = "wget")
}


#####


# Read the STRING protein-protein interactions
String_ppi <- data.table::fread(file = "Databases/StringDB/9606.protein.links.detailed.v12.0.onlyAB.txt.gz", header = TRUE, sep = " ")


# Filter all interactions with scores greater than 400 combined score
String_ppi <- String_ppi %>% filter(combined_score > 400)


# Get list of proteins in the PPI
String_proteins <- sort(unique(c(String_ppi$protein1, String_ppi$protein2)))


#####


# Create mapping from Ensembl Peptide ID ID to Ensembl Gene ID
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
ensemblPeptideID_2_ensemblGeneID <- getBM(attributes = c("ensembl_gene_id", "ensembl_peptide_id", "external_gene_name"),
                                          mart = ensembl,
                                          filters = "ensembl_peptide_id",
                                          values = gsub("^9606.", "", String_proteins))
ensemblPeptideID_2_ensemblGeneID$ensembl_peptide_id <- paste0("9606.", ensemblPeptideID_2_ensemblGeneID$ensembl_peptide_id)


# Merge the IDs
String_ppi <- String_ppi %>% 
  left_join(ensemblPeptideID_2_ensemblGeneID, 
            by = c("protein1" = "ensembl_peptide_id")) %>%
  dplyr::rename("Node1_string_protein_id" = "protein1", "Node1_ensembl_gene_id" = "ensembl_gene_id", "Node1_gene_name" = "external_gene_name") %>% 
  left_join(ensemblPeptideID_2_ensemblGeneID, 
            by = c("protein2" = "ensembl_peptide_id")) %>%
  dplyr::rename("Node2_string_protein_id" = "protein2", "Node2_ensembl_gene_id" = "ensembl_gene_id", "Node2_gene_name" = "external_gene_name") %>%
  filter( !(is.na(Node1_ensembl_gene_id) | is.na(Node2_ensembl_gene_id)) ) 


#####


# Add additional attributes
String_ppi <- String_ppi %>% 
  mutate("Node1_type" = "gene", 
         "Node2_type" = "gene",
         "Edge_type" = "undirected")


# Create node attributes
tmp1 <- String_ppi %>% 
  dplyr::select(c("Node1_ensembl_gene_id", "Node1_type", "Node1_string_protein_id", "Node1_gene_name"))
colnames(tmp1) <- c("Node_ensembl_gene_id", "Node_type", "Node_string_protein_id", "Node_gene_name")

tmp2 <- String_ppi %>% 
  dplyr::select(c("Node2_ensembl_gene_id", "Node2_type", "Node2_string_protein_id", "Node2_gene_name"))
colnames(tmp2) <- c("Node_ensembl_gene_id", "Node_type", "Node_string_protein_id", "Node_gene_name")


node_attributes <- unique(rbind(tmp1, tmp2)) %>% 
  group_by(Node_ensembl_gene_id) %>% 
  summarise(Node_type = paste(unique(Node_type), collapse = "|"), 
            Node_string_protein_id = paste(unique(Node_string_protein_id), collapse = "|"),
            Node_gene_name = paste(unique(Node_gene_name), collapse = "|"))
rm(list = c("tmp1", "tmp2"))


#####


# Convert to graph object

String_ppi_Net <- String_ppi %>% 
  dplyr::select(!c("Node1_string_protein_id", "Node2_string_protein_id", "Node1_type", "Node2_type", "Node1_gene_name", "Node2_gene_name")) %>%
  distinct() %>%
  dplyr::select(c("Node1_ensembl_gene_id", "Node2_ensembl_gene_id", everything()))

String_ppi_Net <- graph_from_data_frame(String_ppi_Net,
                                        directed = FALSE, 
                                        vertices = node_attributes)

String_ppi_Net <- igraph::simplify(  String_ppi_Net, 
                                     remove.loops = TRUE, 
                                     remove.multiple	= TRUE,
                                     edge.attr.comb = list("mean", 
                                                           "Edge_type" = function(x){unique(x)})  ) # remove loops and multi-edges


# Find the connected components
ppi_comps <- components(String_ppi_Net)

# Identify the largest component
largest_ppi_comp <- which.max(ppi_comps$csize)

# Extract the largest component
String_ppi_Net <- induced_subgraph(String_ppi_Net, which(ppi_comps$membership == largest_ppi_comp))

print(paste("Network size (vertices, edges):", vcount(String_ppi_Net), ecount(String_ppi_Net)))

if(!dir.exists("IOFiles/Networks/")){dir.create("IOFiles/Networks/", recursive = TRUE)}
saveRDS(String_ppi_Net, "IOFiles/Networks/STRING_PPI_Net.rds")


#####


# Calculate the network parameters
network_parameters <- data.frame(t(data.frame(nodeNumber = vcount(String_ppi_Net),
                                              edgeNumber = ecount(String_ppi_Net),
                                              diameter = diameter(graph = String_ppi_Net, directed = FALSE),
                                              radius = radius(graph = String_ppi_Net),
                                              averagePathLength = mean_distance(graph = String_ppi_Net),
                                              clusteringCoefficient = transitivity(graph = String_ppi_Net, type = "global", isolates = "zero"),
                                              density = edge_density(graph = String_ppi_Net),
                                              connectedComponents = count_components(graph = String_ppi_Net),
                                              largest_component_size = components(String_ppi_Net)$csize[1],
                                              powerFit = fit_power_law(x = degree(String_ppi_Net)+1)$alpha,
                                              degree = centr_degree(graph = String_ppi_Net)$centralization,
                                              closeness = centr_clo(graph = String_ppi_Net)$centralization,
                                              betweenness = centr_betw(graph = String_ppi_Net)$centralization)))

network_parameters <- tibble::rownames_to_column(network_parameters, "Parameters")
colnames(network_parameters) <- c("Parameters", "Value")
network_parameters$Value <- round(network_parameters$Value, 3)


#####


if(!dir.exists("IOFiles/Networks/")){
  dir.create("IOFiles/Networks/", recursive = TRUE)
}
write.csv(network_parameters, "IOFiles/Networks/STRING_PPI_Net_params.csv", row.names = FALSE)


#####


print(warnings())