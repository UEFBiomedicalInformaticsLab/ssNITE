set.seed(5081)


# Script to summarise the similarity between SSNs  (for publication)


# Load libraries
library(unixtools)
library(data.table)
library(tidyverse)
library(parallel)

# Set temporary directory
if(!dir.exists("tmp_dir/")){dir.create("tmp_dir/", recursive = TRUE)}
set.tempdir("tmp_dir/")

if(!dir.exists(paste0("IOFiles/Compare_X_methods/TCGABRCA_prunedSSN_similarity/"))){
  dir.create({paste0("IOFiles/Compare_X_methods/TCGABRCA_prunedSSN_similarity/")})
}


#####


# Read sample information
print(paste0(Sys.time(), " | Reading the sample information"))

TCGA_data <- readRDS("IOFiles/GDC_data/vstCounts__DA_TCGA__DT_BreastCancer.rds")
sample_info <- TCGA_data$sample_info %>% filter(condition != "SolidTissueNormal_NA") %>% droplevels()
sample_info <- sample_info %>%
  select(c("barcode", "patient", "race", "paper_pathologic_stage", "paper_BRCA_Subtype_PAM50"))  %>% 
  rename(c("sample_id" = "barcode", 
           "Stage" = "paper_pathologic_stage", 
           "Subtype" = "paper_BRCA_Subtype_PAM50", 
           "Race" = "race"))
row.names(sample_info) <- NULL

sample_info[which(sample_info$Race == "not reported"), "Race"] <- NA
remove_race <- names(which(table(sample_info$Race) < 10))
sample_info[which(sample_info$Race %in% remove_race), "Race"] <- NA
rm(list = c("TCGA_data", "remove_race"))


#####


# Add the information about histological marker status

print(paste0(Sys.time(), " | Adding the histological marker status for the samples"))

clinicalInfo_patient <- readRDS("IOFiles/GDC_data/TCGABRCA_data.rds")
clinicalInfo_patient <- clinicalInfo_patient$clinicalInfo_patient
clinicalInfo_patient <- clinicalInfo_patient %>% 
  select(c("bcr_patient_barcode", 
           "breast_carcinoma_estrogen_receptor_status",
           "breast_carcinoma_progesterone_receptor_status", 
           "lab_proc_her2_neu_immunohistochemistry_receptor_status")) %>% 
  rename(c("ER status" = "breast_carcinoma_estrogen_receptor_status", 
           "PR status" = "breast_carcinoma_progesterone_receptor_status", 
           "Her2 status" = "lab_proc_her2_neu_immunohistochemistry_receptor_status")) %>%
  mutate(across(c("ER status", "PR status", "Her2 status"), 
                ~ ifelse(. == "Positive", 1, 0)))

sample_info <- sample_info %>%
  left_join(clinicalInfo_patient, 
            by = c("patient" = "bcr_patient_barcode"))
rm(clinicalInfo_patient)


#####

sample_info <- as.data.frame(lapply(sample_info, function(x){ factor(x, exclude = NULL) }), check.names = FALSE)

#####


# Function to make symmetric matrix 
make_symmetric_matrix <- function(df){
  
  sample_id_list <- sort(unique(c(df$Network1, df$Network2)))
  
  df <-  df %>%
    bind_rows(df %>% 
                rename(Network2 = Network1,
                       Network1 = Network2))
  
  df <- df %>% 
    bind_rows(data.frame("Network1" = sample_id_list, 
                         "Network2" = sample_id_list, 
                         "Similarity" = 1))
  
  df <- df %>% 
    pivot_wider(id_cols = "Network1", 
                names_from = "Network2", 
                values_from = "Similarity")  %>% 
    column_to_rownames("Network1") %>% 
    as.matrix()
  
  df <- df[sample_id_list, sample_id_list]
  
  return(df)
  
}


#####


# Create the matrix to store results
print(paste0(Sys.time(), " | Creating the test grid"))
test_grid <- expand_grid("network_type" = c("ssNITEpy", "BONOBO", "CSN", "LIONESS", "SWEET", "CSNsparse", "SWEETsparse"), 
                         "similarity_type" = c("Jaccard similarity (edge)", "Jaccard similarity (top degree)", "Adj. Rand Index (communities)"), 
                         "group_type" = c("Subtype", "ER status", "PR status", "Her2 status"))


#####


# Summarise the similarity

print(paste0(Sys.time(), " | Checking similarity vs clinical info relation"))

result <- mclapply(X = 1:nrow(test_grid), 
                   mc.cores = 15, 
                   FUN = function(i){
                     
                     network_type_select <- test_grid[i, ]$network_type
                     similarity_type_select <- test_grid[i, ]$similarity_type
                     group_type_select <- test_grid[i, ]$group_type
                     
                     
                     # Get the similarity data
                     sim_data <- readRDS(paste0("IOFiles/SSN_similarity/TCGABRCA/undir_", network_type_select, "/", network_type_select, "_prunedSSNsimilarity__DA_TCGA__DT_BreastCancer.rds"))
                     sim_data <- sim_data %>% 
                       select(c("Network1", "Network2", all_of(similarity_type_select))) %>% 
                       rename("Similarity" = similarity_type_select)
                     
                     if(all(is.na(sim_data$Similarity))){
                       
                       # Fill the results
                       tmp1 <- data.frame("Network" = network_type_select, 
                                          "Similarity" = similarity_type_select, 
                                          "Grouping" = group_type_select, 
                                          "Modularity" = NA, 
                                          "Silhoutte" = NA, 
                                          "ANOSIM (stat)" = NA, 
                                          "ANOSIM (signif)" = NA, 
                                          check.names = FALSE)
                       
                     }else{
                       
                       sim_data <- make_symmetric_matrix(sim_data)
                       
                       if(similarity_type_select == "Adj. Rand Index (communities)"){
                         sim_data[sim_data < 0] <- 0
                       }
                       
                       if(!all(row.names(sim_data) == colnames(sim_data))){ stop("---") }
                       
                       
                       # Get the groups
                       groups <- sample_info[[group_type_select]]
                       names(groups) <- sample_info$sample_id
                       groups <- groups[row.names(sim_data)]
                       
                       
                       # Calculate modularity
                       net <- igraph::graph_from_adjacency_matrix(adjmatrix = sim_data, mode = "upper", weighted = TRUE, diag = FALSE)
                       modularity_res <- igraph::modularity(x = net, membership = as.numeric(groups), weights = igraph::E(net)$weight)
                       
                       
                       # Calculate silhouette score
                       silhouette_res <- cluster::silhouette(as.numeric(groups), dist = as.dist(1 - sim_data))
                       silhouette_res <- summary(silhouette_res)
                       silhouette_res <- silhouette_res$si.summary[["Mean"]]
                       
                       
                       # ANOSIM
                       anosim_res <- vegan::anosim(x = (1 - sim_data), grouping = groups, permutations = 1000)  
                       
                       
                       # Fill the results
                       tmp1 <- data.frame("Network" = network_type_select, 
                                          "Similarity" = similarity_type_select, 
                                          "Grouping" = group_type_select, 
                                          "Modularity" = modularity_res, 
                                          "Silhoutte" = silhouette_res, 
                                          "ANOSIM (stat)" = anosim_res$statistic, 
                                          "ANOSIM (signif)" = anosim_res$signif, 
                                          check.names = FALSE)
                     }
                     
                     return(tmp1)
                     
                   })


result <- result %>% bind_rows()
write.csv(result, "IOFiles/Compare_X_methods/TCGABRCA_prunedSSN_similarity/TCGABRCA_prunedSSN_similarityVSgroups.csv", row.names = FALSE)


#####


# result <- read.csv("IOFiles/Compare_X_methods/TCGABRCA_prunedSSN_similarity/TCGABRCA_prunedSSN_similarityVSgroups.csv", check.names = FALSE)

plot_data <- result %>% 
  select(c("Network", "Similarity", "Grouping", "ANOSIM (stat)", "ANOSIM (signif)")) %>% 
  mutate( sig_label = ifelse(`ANOSIM (signif)` <= 0.05, "*", "") )

plot_data$Network <- factor(x = plot_data$Network, 
                            levels = c("ssNITEpy", "BONOBO", "CSN", "LIONESS", "SWEET", "CSNsparse", "SWEETsparse")) 

plot_data$Grouping <- factor(x = plot_data$Grouping, 
                             levels = c("Subtype", "ER status", "PR status", "Her2 status")) 

plot_data$Similarity <- factor(x = plot_data$Similarity, levels = unique(plot_data$Similarity))

plot <- ggplot(plot_data, aes(x = Grouping, y = Network, fill = `ANOSIM (stat)`)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sig_label), color = "black", size = 2, vjust = 0.75, na.rm = TRUE) +
  facet_grid( ~ Similarity, scales = "free") +
  scale_y_discrete(labels = c("ssNITEpy" = "ssNITE")) +
  scale_fill_gradient2(low = "#2166AC", 
                       mid = "white", 
                       high = "#B2182B",
                       midpoint = 0,
                       # limits = c(-1, 1), 
                       na.value = "#F5F5F5") +  
  labs( x = "Clinical groups",
        y = "SSN type",
        fill = "R" ) +
  theme( panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.25, linetype = NULL),
         panel.grid = element_blank(),
         strip.background = element_blank(),
         strip.text = element_text(size = 4, margin = margin(t = 1, b = 1, r = 1, l= 1)),
         text = element_text(size = 5),
         axis.text.x = element_text(size = 4, angle = 45, vjust = 1, hjust = 1, margin = margin(t = 1)),
         axis.text.y = element_text(size = 4, margin = margin(r = 1)),
         axis.ticks.length = unit(0, "cm"),
         legend.position = "bottom",
         legend.box.spacing = unit(0.05, "cm"),
         legend.key = element_rect(fill = NA),
         legend.key.size = unit(0.1, "cm"),
         legend.key.width = unit(0.2, "cm"),
         legend.key.spacing.y = unit(0.5, "cm"),
         legend.title = element_text(size = 2, face = "bold", margin = margin(t = 0, r = 1, b = 0, l = 0)),
         legend.text = element_text(size = 2, margin = margin(t = 1, r = 0, b = 0, l = 0)),
         legend.margin = margin(t = 1, r = 1, b = 1, l = 1),
         legend.spacing = unit(0.2, "cm") )

tiff(paste0("IOFiles/Compare_X_methods/TCGABRCA_prunedSSN_similarity/TCGABRCA_prunedSSN_similarityVSgroups.tiff"),
     width = 8, height = 4,
     units = "cm", compression = "lzw", res = 1200)

plot

dev.off()


#####


# Summarise the similarity between SSN within each method

print(paste0(Sys.time(), " | Summarising similarity"))

network_type_list <- c("ssNITEpy", "BONOBO", "CSN", "LIONESS", "SWEET", "CSNsparse", "SWEETsparse")
sim_data <- lapply(network_type_list, function(x){
  tmp1 <- readRDS(paste0("IOFiles/SSN_similarity/TCGABRCA/undir_", x, "/", x, "_prunedSSNsimilarity__DA_TCGA__DT_BreastCancer.rds"))
  tmp1 <- tmp1 %>% rename(c("Jaccard Similarity (node)" = "Jaccard similarity (node)",
                            "Jaccard Similarity (edge)" = "Jaccard similarity (edge)", 
                            "Jaccard Similarity (top degree)" = "Jaccard similarity (top degree)"))
  tmp1
})
names(sim_data) <- network_type_list
sim_data <- rbindlist(sim_data, idcol = "ssn_type")

sel_cols <- setdiff(colnames(sim_data), c("ssn_type", "Network1", "Network2"))
summary_results <- sim_data %>% 
  group_by(ssn_type) %>% 
  summarise(across( all_of(sel_cols), 
                    list( "minimum" = ~if(all(is.na(.x))){NA}else{min(.x, na.rm = TRUE)}, 
                          "maximum" = ~if(all(is.na(.x))){NA}else{max(.x, na.rm = TRUE)},
                          "median" = ~if(all(is.na(.x))){NA}else{median(.x, na.rm = TRUE)},
                          "mean" = ~if(all(is.na(.x))){NA}else{mean(.x, na.rm = TRUE)},
                          "sd" = ~if(all(is.na(.x))){NA}else{sd(.x, na.rm = TRUE)} ),
                    .names = "{.col}___{.fn}" ))


summary_results <- summary_results %>%
  pivot_longer(-ssn_type, 
               names_to = "metric", 
               values_to = "score") %>%
  separate(col = "metric", 
           into = c("similarity_type", "stat"), 
           sep = "___") %>% 
  pivot_wider(id_cols = c("ssn_type", "similarity_type"), 
              names_from = "stat", 
              values_from = "score")

write.csv(summary_results, "IOFiles/Compare_X_methods/TCGABRCA_prunedSSN_similarity/TCGABRCA_prunedSSN_summaryOFsimilarity.csv", row.names = FALSE)


#####


# Plot the similarity as box plots
plot_data <- sim_data %>% 
  select(!c("Jaccard Similarity (node)", )) %>% 
  pivot_longer(cols = c("Jaccard Similarity (edge)", "Jaccard Similarity (top degree)", "Adj. Rand Index (communities)"), 
               names_to = "metric", values_to = "score") 

plot_data$ssn_type <- factor(plot_data$ssn_type, levels = c("ssNITEpy", "BONOBO", "CSN", "LIONESS", "SWEET", "CSNsparse", "SWEETsparse"))
plot_data$metric <- factor(plot_data$metric, levels = c("Jaccard Similarity (edge)", "Jaccard Similarity (top degree)", "Adj. Rand Index (communities)"))


plot <- ggplot(plot_data, aes(x = metric, y = score, fill = ssn_type)) +
  geom_boxplot(position = position_dodge(preserve = "single"),
               width = 0.75,
               lwd = 0.1,
               outliers = FALSE,
               na.rm = TRUE)  +
  labs(x = "Metric", y = "Score", fill = "SSN type") +
  scale_fill_manual(label = c("ssNITEpy" = "ssNITE"), 
                    values = c("ssNITEpy" = "#33A02C", "ssNITE" = "#B2DF8A", 
                               "BONOBO" = "#1F78B4", "BONOBOsparse" = "#A6CEE3", 
                               "CSN" = "#E72F31", "CSNsparse" = "#FB9A99", 
                               "SWEET" = "#CC6600", "SWEETsparse" = "#FDBF6F", 
                               "LIONESS" = "#6A3D9A"),  
                    drop = FALSE) +
  theme(panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.25, linetype = NULL),
        panel.grid = element_blank(),
        panel.spacing = unit(0.1, "cm"),
        text = element_text(size = 6),
        axis.text.x = element_text(size = 4, angle = 0, vjust = 0, hjust = 0.5),
        axis.ticks = element_line(colour = "black", linewidth = 0.2),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.key = element_blank(),
        legend.key.spacing.x = unit(0.01, "cm"),
        legend.key.spacing.y = unit(0.01, "cm"),
        legend.key.size = unit(0.1, "cm"),
        legend.title = element_text(size = 4, margin = margin(t = 1, r = 1, b = 1, l = 1)),
        legend.text = element_text(size = 3, margin = margin(t = 1, r = 2, b = 1, l = 1)),
        legend.margin = margin(1,1,1,1),
        legend.box.spacing = unit(0.1, "cm"),
        legend.box.background = element_rect(colour = "black", linewidth = 0.25))


tiff(paste0("IOFiles/Compare_X_methods/TCGABRCA_prunedSSN_similarity/TCGABRCA_prunedSSN_similarity_boxplot.tiff"),
     width = 8, height = 4,
     units = "cm", compression = "lzw", res = 1200)

plot

dev.off()


#####


print(warnings())