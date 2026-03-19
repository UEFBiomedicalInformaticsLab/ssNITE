set.seed(5081)


# Script to plot summary table generated from prediction of TCGABRCA sample information from network topological data and gene expression (for publication)


# Load libraries
library(unixtools)
library(data.table)
library(tidyverse)
library(parallel)

# Set temporary directory
if(!dir.exists("tmp_dir/")){dir.create("tmp_dir/", recursive = TRUE)}
set.tempdir("tmp_dir/")

if(!dir.exists(paste0("IOFiles/Compare_X_methods/TCGABRCA_summariseMLres_SSNtopol_vs_sampleInfo/"))){
  dir.create({paste0("IOFiles/Compare_X_methods/TCGABRCA_summariseMLres_SSNtopol_vs_sampleInfo/")})
}


#####


# Set the network list
network_type_list <- c("ssNITEpy", "BONOBO", "CSN", "LIONESS", "SWEET", "CSNsparse", "SWEETsparse")


#####


# Read the accuracy scores from predictions by using network level topology

netTopol_ml <- lapply(network_type_list, function(network_type){
  tmp1 <- readRDS(paste0("IOFiles/ML_SSNtopol_vs_sampleInfo/TCGABRCA/undir_", network_type, "/MLres_ssnTopol_vs_sampleInfo__DA_TCGA__DT_BreastCancer.rds"))
  tmp1 <- tmp1[, c("ml_target", "n_repeat", "n_fold", "accuracy_score", "accuracy_metric")]
  
  tmp1 <- tmp1[, unlist(lapply(.SD, function(x){ list(median = median(x),
                                                      mean = mean(x), 
                                                      sd = sd(x), 
                                                      iqr = IQR(x)) }), 
                        recursive = FALSE), 
               by = .(ml_target, accuracy_metric), 
               .SDcols = !c("n_repeat", "n_fold") ]
  tmp1
})
names(netTopol_ml) <- network_type_list
netTopol_ml <- rbindlist(netTopol_ml, idcol = "network_type")

plot_data <- netTopol_ml %>% filter(ml_target %in% c("Subtype", "ER status", "PR status", "Her2 status", "ESTIMATE", "LUMP"))
plot_data$network_type <- factor(plot_data$network_type , 
                                 levels = c("ssNITEpy", "BONOBO", "CSN", "LIONESS", "SWEET", "CSNsparse", "SWEETsparse"))
plot_data$ml_target <- factor(plot_data$ml_target , levels = c("Subtype", "ER status", "PR status", "Her2 status", "ESTIMATE", "LUMP"))


plot <- ggplot(plot_data, aes(x = ml_target, y = network_type, fill = `accuracy_score.median`)) +
  geom_tile(color = "white", linewidth = 0.05, na.rm = TRUE) +
  geom_text(aes(label = round(`accuracy_score.median`, 2)), size = 1, vjust = -0.2, na.rm = TRUE) +
  geom_text(aes(label = paste0("(IQR = ", round(accuracy_score.iqr, 2), ")") ), size = 0.6, vjust = 2, na.rm = TRUE) +  
  # facet_wrap(~accuracy_metric, scales = "free_x") + 
  labs(x = "Prediction target", 
       y = "Data source", 
       fill = "Score") +
  scale_y_discrete(labels = c("ssNITEpy" = "ssNITE")) +
  scale_fill_distiller(palette = "Blues", direction = 1, na.value = "grey95") +  
  theme( panel.background = element_blank(),
         panel.grid = element_blank(),
         strip.background = element_rect(color = "black", linewidth = 0.25,),
         strip.text = element_text(size = 5, face = "bold", margin = margin(t = 1, b = 1, r = 1, l= 1)),
         text = element_text(size = 7),
         # plot.title = element_text(face = "bold", hjust = 0.5, size = 8),
         axis.text.x = element_text(size = 4, angle = 45, vjust = 1, hjust = 1, margin = margin(t = 0, b = 0, r = 0, l= 0)),
         axis.text.y = element_text(size = 4, margin = margin(t = 0, b = 0, r = 0, l= 0)),
         axis.ticks.length = unit(0, "cm"),
         legend.position = "right",
         legend.key = element_rect(fill = NA),
         legend.key.size = unit(0.3, "cm"),
         legend.key.spacing.y = unit(0.1, "cm"),
         legend.title = element_text(size = 4, face = "bold", margin = margin(t = 1, r = 0, b = 1.5, l = 0)),
         legend.text = element_text(size = 4, margin = margin(t = 0.5, r = 0, b = 0.1, l = 1)),
         legend.margin = margin(1,1,1,1),
         legend.spacing = unit(0.2, "cm") )

write.csv(netTopol_ml, paste0("IOFiles/Compare_X_methods/TCGABRCA_summariseMLres_SSNtopol_vs_sampleInfo/TCGABRCA_allMLres_ssnTopol_vs_sampleInfo.csv"), row.names = FALSE)
tiff(paste0("IOFiles/Compare_X_methods/TCGABRCA_summariseMLres_SSNtopol_vs_sampleInfo/TCGABRCA_allMLres_ssnTopol_vs_sampleInfo.tiff"),
     width = 8, height = 6,
     units = "cm", compression = "lzw", res = 1200)

plot

dev.off()


#####


# Read the accuracy scores from predictions by using node centralities

geneExpr_ml <- readRDS("IOFiles/ML_SSNtopol_vs_sampleInfo/TCGABRCA/MLres_geneExpr_vs_sampleInfo__DA_TCGA__DT_BreastCancer.rds")
geneExpr_ml <- geneExpr_ml[, c("ml_target", "n_repeat", "n_fold", "accuracy_score", "accuracy_metric")]

geneExpr_ml <- geneExpr_ml[, unlist(lapply(.SD, function(x){ list(median = median(x),
                                                                  mean = mean(x), 
                                                                  sd = sd(x), 
                                                                  iqr = IQR(x)) }), 
                                    recursive = FALSE), 
                           by = .(ml_target, accuracy_metric), 
                           .SDcols = !c("n_repeat", "n_fold") ]

geneExpr_ml[, "centrality_type" := "Gene expr."]
geneExpr_ml[, "network_type" := "Gene expr."]



# Read the accuracy scores from predictions by using node centralities

nodeCentral_ml <- lapply(network_type_list, function(network_type){
  tmp1 <- readRDS(paste0("IOFiles/ML_SSNtopol_vs_sampleInfo/TCGABRCA/undir_", network_type, "/MLres_ssnCentPC_vs_sampleInfo__DA_TCGA__DT_BreastCancer.rds"))
  tmp1 <- tmp1[, c("centrality_type", "ml_target", "n_repeat", "n_fold", "accuracy_score", "accuracy_metric")]
  
  tmp1 <- tmp1[, unlist(lapply(.SD, function(x){ list(median = median(x), 
                                                      mean = mean(x), 
                                                      sd = sd(x), 
                                                      iqr = IQR(x)) }), 
                        recursive = FALSE), 
               by = .(centrality_type, ml_target, accuracy_metric), 
               .SDcols = !c("n_repeat", "n_fold") ]
  tmp1
})
names(nodeCentral_ml) <- network_type_list
nodeCentral_ml <- rbindlist(nodeCentral_ml, idcol = "network_type")


plot_data <- rbindlist(list(nodeCentral_ml, geneExpr_ml), use.names = TRUE)
plot_data <- plot_data %>% filter(ml_target %in% c("Subtype", "ER status", "PR status", "Her2 status", "ESTIMATE", "LUMP"))
plot_data$network_type <- factor(plot_data$network_type , levels = c(network_type_list, "Gene expr."))
plot_data$centrality_type <- factor(plot_data$centrality_type , 
                                    levels = c("Degree", "Betweenness", "Closeness", "Clustering coef.", "Eigen centrality", 
                                               "Page rank", "Edge betweenness", "NNEC", "Gene expr."))
plot_data$ml_target <- factor(plot_data$ml_target , levels = c("Subtype", "ER status", "PR status", "Her2 status", "ESTIMATE", "LUMP"))


plot <- ggplot(plot_data, aes(x = ml_target, y = network_type, fill = `accuracy_score.median`)) +
  geom_tile(color = "white", linewidth = 0.05, na.rm = TRUE) +
  geom_text(aes(label = round(`accuracy_score.median`, 2)), size = 1, vjust = -0.2, na.rm = TRUE) +
  geom_text(aes(label = paste0("(IQR = ", round(accuracy_score.iqr, 2), ")") ), size = 0.6, vjust = 2, na.rm = TRUE) +  
  facet_wrap(~ centrality_type, scales = "free_y", space = "fixed") +
  labs(x = "Prediction target", 
       y = "Data source", 
       fill = "Score") +
  scale_y_discrete(labels = c("ssNITEpy" = "ssNITE")) +
  scale_fill_distiller(palette = "Blues", direction = 1, na.value = "grey95") +  
  theme( panel.background = element_blank(),
         panel.grid = element_blank(),
         strip.background = element_blank(),
         strip.text = element_text(size = 5, face = "bold", margin = margin(t = 1, b = 1, r = 1, l= 1)),
         text = element_text(size = 7),
         # plot.title = element_text(face = "bold", hjust = 0.5, size = 8),
         axis.text.x = element_text(size = 4, angle = 45, vjust = 1, hjust = 1, margin = margin(t = 0, b = 0, r = 0, l= 0)),
         axis.text.y = element_text(size = 4, margin = margin(t = 0, b = 0, r = 0, l= 0)),
         axis.ticks.length = unit(0, "cm"),
         legend.position = "right",
         legend.key = element_rect(fill = NA),
         legend.key.size = unit(0.4, "cm"),
         legend.key.spacing.y = unit(0.1, "cm"),
         legend.title = element_text(size = 4, face = "bold", margin = margin(t = 1, r = 0, b = 1.5, l = 0)),
         legend.text = element_text(size = 4, margin = margin(t = 0.5, r = 0, b = 0.1, l = 1)),
         legend.margin = margin(1,1,1,1),
         legend.spacing = unit(0.2, "cm") )

write.csv(nodeCentral_ml, paste0("IOFiles/Compare_X_methods/TCGABRCA_summariseMLres_SSNtopol_vs_sampleInfo/TCGABRCA_allMLres_ssnCentPC_vs_sampleInfo.csv"), row.names = FALSE)
tiff(paste0("IOFiles/Compare_X_methods/TCGABRCA_summariseMLres_SSNtopol_vs_sampleInfo/TCGABRCA_allMLres_ssnCentPC_vs_sampleInfo.tiff"),
     width = 18, height = 10,
     units = "cm", compression = "lzw", res = 1200)

plot

dev.off()


#####


# Generate winner plot

netTopol_ml[, "centrality_type" := "Net. topol."]
plot_data <- rbindlist(list(nodeCentral_ml, netTopol_ml, geneExpr_ml), use.names = TRUE, fill = TRUE)


plot_data <- plot_data %>%
  filter(centrality_type != "Gene expr.") %>%
  filter(ml_target %in% c("Subtype", "ER status", "PR status", "Her2 status", "ESTIMATE", "LUMP")) %>%
  group_by(ml_target, centrality_type) %>%
  slice_max(order_by = accuracy_score.median, n = 1) %>% 
  ungroup()
plot_data$label_score <- paste0(round(plot_data$accuracy_score.median, 2), " (", round(plot_data$accuracy_score.iqr, 2), ")")
plot_data$label_net <- recode(plot_data$network_type, "ssNITEpy" = "ssNITE")

plot_data$centrality_type <- factor(plot_data$centrality_type,
                                    levels = c("Betweenness", "Closeness", "Clustering coef.", "Degree", 
                                               "Eigen centrality", "Page rank", "Edge betweenness",  "NNEC", "Net. topol."))

plot_data$ml_target <- factor(plot_data$ml_target , 
                              levels = c("Subtype", "ER status", "PR status", "Her2 status", "ESTIMATE", "LUMP"))


# Generate an annotation bar 
top_annot <- plot_data %>%
  group_by(ml_target) %>%
  summarise(accuracy_metric = dplyr::first(accuracy_metric), .groups = "drop") %>% 
  mutate("annot" = "Target type")
top_annot$accuracy_metric <- recode(top_annot$accuracy_metric, "F1" = "Classification (F1)", "R2" = "Regression (R²)")


top_annot <- ggplot(top_annot, aes(x = ml_target, y = annot, fill = accuracy_metric)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = NULL, 
       y = NULL, fill = "Target type (metric): ") +
  theme( plot.margin = margin(t = 0, b = 0, r = 0, l = 0),
         panel.background = element_blank(),
         panel.grid = element_blank(),
         strip.background = element_rect(color = "black", linewidth = 0.25),
         strip.text = element_text(size = 5, face = "bold", margin = margin(t = 1, b = 1, r = 1, l= 1)),
         text = element_text(size = 7),
         axis.text.x = element_blank(),
         axis.text.y = element_text(size = 3, margin = margin(t = 0, b = 0, r = 0, l= 0)),
         axis.ticks.length = unit(0, "cm"),
         legend.position = "bottom" )


# generate the main annotation plot
main_plot <- ggplot(plot_data, aes(x = ml_target, y = centrality_type, fill = network_type)) +
  geom_tile(color = "white", linewidth = 0.5, show.legend = FALSE) +
  geom_text(aes(label = label_net), size = 1, color = "black", fontface = "bold") +
  geom_text(aes(label = label_score), size = 1, vjust = 2.5, na.rm = TRUE) +
  # facet_wrap(~accuracy_metric, scale = "free_x") +
  labs(y = "Centrality type", 
       x = "Prediction target") +
  scale_fill_manual(values = c("ssNITEpy" = "#33A02C", "ssNITE" = "#B2DF8A", 
                               "BONOBO" = "#1F78B4", "BONOBOsparse" = "#A6CEE3", 
                               "CSN" = "#E72F31", "CSNsparse" = "#FB9A99", 
                               "SWEET" = "#CC6600", "SWEETsparse" = "#FDBF6F", 
                               "LIONESS" = "#6A3D9A"),  
                    drop = FALSE) +  
  theme( plot.margin = margin(t = 0, b = 0),
         panel.background = element_blank(),
         panel.grid = element_blank(),
         strip.background = element_rect(color = "black", linewidth = 0.25,),
         strip.text = element_text(size = 5, face = "bold", margin = margin(t = 1, b = 1, r = 1, l= 1)),
         text = element_text(size = 7),
         axis.text.x = element_text(size = 4, angle = 0, vjust = 1, hjust = 0.5, margin = margin(t = 0, b = 0, r = 0, l= 0)),
         axis.text.y = element_text(size = 4, margin = margin(t = 0, b = 0, r = 0, l= 0)),
         axis.ticks.length = unit(0, "cm"),
         legend.position = "none")


require(patchwork)
plot <- top_annot / main_plot + 
  plot_layout(heights = c(0.02, 1), guides = "collect") & 
  theme(legend.position = "bottom",
        legend.box.spacing = unit(0.1, "cm"),
        legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0),
        legend.key.size = unit(0.3, "cm"),
        legend.key.spacing.x = unit(0.1, "cm"),
        legend.title = element_text(size = 3, margin = margin(t = 0, r = 1, b = 0, l = 0)),
        legend.text = element_text(size = 3, margin = margin(t = 0, r = 0, b = 0, l = 0.5)),
        legend.margin = margin(t = 0, r = 1, b = 0, l = 1),
        legend.spacing = unit(0.2, "cm") )

tiff(paste0("IOFiles/Compare_X_methods/TCGABRCA_summariseMLres_SSNtopol_vs_sampleInfo/TCGABRCA_winnerPlot_sampleInfo_pred.tiff"),
     width = 8, height = 8,
     units = "cm", compression = "lzw", res = 1200)

plot

dev.off()


#####


print(warnings())