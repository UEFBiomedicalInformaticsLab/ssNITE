set.seed(5081)


# Script to summarise the corelation between SSN topological properties and sample information (for publication)


# Load libraries
library(unixtools)
library(data.table)
library(tidyverse)
library(parallel)

# Set temporary directory
if(!dir.exists("tmp_dir/")){dir.create("tmp_dir/", recursive = TRUE)}
set.tempdir("tmp_dir/")

if(!dir.exists(paste0("IOFiles/Compare_X_methods/TCGABRCA_summariseCorel_SSNtopol_vs_sampleInfo/"))){
  dir.create({paste0("IOFiles/Compare_X_methods/TCGABRCA_summariseCorel_SSNtopol_vs_sampleInfo/")})
}


#####


# Set the network list
network_type_list <- c("ssNITEpy", "BONOBO", "CSN", "LIONESS", "SWEET", "CSNsparse", "SWEETsparse")


#####


# Read the corelation results 

netTopol_corel <- lapply(network_type_list, function(network_type){
  tmp1 <- read.csv(paste0("IOFiles/Corel_SSNtopol_vs_sampleInfo/TCGABRCA/undir_", network_type, "/", network_type, "_TopoSampleCor__DA_TCGA__DT_BreastCancer.csv"), check.names = FALSE, row.names = c(1))
  tmp1 <- tmp1 %>% rownames_to_column("network_topol")
  tmp1
})
names(netTopol_corel) <- network_type_list
netTopol_corel <- rbindlist(netTopol_corel, idcol = "network_type")


netTopol_corel <- netTopol_corel %>% 
  filter( network_topol %in% c("Density", "Connected component count", "Largest component size", 
                               "Number of isolates", "Transitivity", "Assortativity", 
                               "Mean distance", "Diameter", "Modularity", 
                               "Number of modules", "Gini coefficient", "Hill index") ) %>% 
  select( c("network_type", "network_topol", 
            "Subtype::Basal", "Subtype::Her2", "Subtype::LumA", "Subtype::LumB", 
            "ER status", "PR status", "Her2 status", "ESTIMATE", "LUMP") ) 


#####


tmp1 <- str_wrap(colnames(netTopol_corel[,-c(1:2)]), 30)
plot_data <- netTopol_corel %>% 
  pivot_longer(-c("network_type", "network_topol"), 
               names_to = "sample_info", 
               values_to = "correlation")

plot_data$network_type <- factor(plot_data$network_type, levels = c("ssNITEpy", "BONOBO", "CSN", "LIONESS", "SWEET", "CSNsparse", "SWEETsparse"))
plot_data$sample_info <- factor(str_wrap(plot_data$sample_info, 30), levels = tmp1)
plot_data$network_topol <- factor(plot_data$network_topol, level = unique(plot_data$network_topol))
rm(tmp1)

plot <- ggplot(plot_data, aes(x = sample_info, y = network_topol, fill = correlation, label = round(correlation, 2))) +
  geom_tile() +
  # geom_text(size = 1,
  #           na.rm = TRUE) +
  facet_wrap(~network_type, nrow = 2, labeller = labeller(network_type = c("ssNITEpy" = "ssNITE"))) +
  labs(fill = "Spearman ρ", 
       x = "Sample info.", 
       y = "Topology") +
  scale_fill_gradient2(low = "#2166AC", 
                       mid = "white", 
                       high = "#B2182B",
                       midpoint = 0,
                       # limits = c(-1, 1),
                       oob = scales::squish, 
                       na.value = "#F5F5F5") +
  theme(panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.25, linetype = NULL),
        panel.grid = element_blank(),
        panel.spacing = unit(0.1, "cm"),
        strip.background = element_blank(),
        strip.text = element_text(size = 5, face = "bold", margin = margin(t = 1, b = 1, r = 1, l= 1)),
        text = element_text(size = 6), 
        plot.title = element_text(size = 5, hjust = 0.5, face = "bold"),
        axis.text.x = element_text(size = 4, angle = 90, vjust = 0.5, hjust = 1), 
        axis.text.y = element_text(size = 4),
        axis.ticks = element_line(colour = "black", linewidth = 0.2),
        legend.position = c(0.9, 0.2),
        legend.title = element_text(size = 4),
        legend.direction = "horizontal",
        legend.key = element_blank(),
        legend.key.size = unit(0.25, "cm"),
        legend.text = element_text(size = 4, margin = margin(t = 2, b = 0, r = 0, l = 0)),
        legend.margin = margin(t = 3, b = 3, r = 3, l = 3),
        legend.box.spacing = unit(0.3, "cm"),
        legend.box.background = element_rect(colour = "black", linewidth = 0.25))


tiff(paste0("IOFiles/Compare_X_methods/TCGABRCA_summariseCorel_SSNtopol_vs_sampleInfo/TCGABRCA_allCorel_SSNtopol_vs_sampleInfo.tiff"),
     width = 20, height = 8,
     units = "cm", compression = "lzw", res = 1200)

plot

dev.off()


#####


# Generate winner plot
plot_data <- netTopol_corel %>%
  pivot_longer(-c("network_type", "network_topol"),
               names_to = "sample_info",
               values_to = "correlation") %>%
  group_by(network_topol, sample_info) %>%
  slice_max(order_by = abs(correlation), n = 1, with_ties = TRUE, na_rm = TRUE) %>%
  ungroup()
plot_data$label_net <- recode(plot_data$network_type, "ssNITEpy" = "ssNITE")


plot_data$network_topol <- factor(plot_data$network_topol,
                                  levels = c("Density", "Connected component count", "Largest component size",
                                             "Number of isolates", "Transitivity", "Assortativity",
                                             "Mean distance", "Diameter", "Modularity",
                                             "Number of modules", "Gini coefficient", "Hill index"))

plot_data$sample_info <- factor(plot_data$sample_info,
                                levels = c("Subtype::Basal", "Subtype::Her2", "Subtype::LumA", "Subtype::LumB",
                                           "ER status", "PR status", "Her2 status", "ESTIMATE", "LUMP"))


plot <- ggplot(plot_data, aes(x = sample_info, y = network_topol, fill = network_type)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = label_net), size = 1, color = "black", fontface = "bold") +
  geom_text(aes(label = round(correlation, 2)), size = 1, vjust = 2.5, na.rm = TRUE) +
  labs(x = "Clinical info.",
       y = "Network topol.") +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 10)) +
  scale_fill_manual(values = c("ssNITEpy" = "#33A02C", "ssNITE" = "#B2DF8A",
                               "BONOBO" = "#1F78B4", "BONOBOsparse" = "#A6CEE3",
                               "CSN" = "#E72F31", "CSNsparse" = "#FB9A99",
                               "SWEET" = "#CC6600", "SWEETsparse" = "#FDBF6F",
                               "LIONESS" = "#6A3D9A"),
                    drop = FALSE) +
  theme( panel.background = element_blank(),
         panel.grid = element_blank(),
         strip.background = element_rect(color = "black", linewidth = 0.25,),
         strip.text = element_text(size = 5, face = "bold", margin = margin(t = 1, b = 1, r = 1, l= 1)),
         text = element_text(size = 7),
         axis.text.x = element_text(size = 4, angle = 0, vjust = 1, hjust = 0.5, margin = margin(t = 0, b = 0, r = 0, l= 0)),
         axis.text.y = element_text(size = 4, margin = margin(t = 0, b = 0, r = 0, l= 0)),
         axis.ticks.length = unit(0, "cm"),
         legend.position = "none")

tiff(paste0("IOFiles/Compare_X_methods/TCGABRCA_summariseCorel_SSNtopol_vs_sampleInfo/TCGABRCA_winnerPlot_sampleInfo_corel.tiff"),
     width = 10, height = 8,
     units = "cm", compression = "lzw", res = 1200)

plot

dev.off()


#####


print(warnings())