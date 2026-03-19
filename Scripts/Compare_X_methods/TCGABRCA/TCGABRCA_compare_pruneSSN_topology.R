set.seed(5081)


# Script to compile the pruned SSN topology from different methods into a single figure (for publication)


# Load libraries
library(unixtools)
library(data.table)
library(tidyverse)

# Set temporary directory
if(!dir.exists("tmp_dir/")){dir.create("tmp_dir/", recursive = TRUE)}
set.tempdir("tmp_dir/")


if(!dir.exists("IOFiles/Compare_X_methods/TCGABRCA_compare_pruneSSN_topology/")){
  dir.create("IOFiles/Compare_X_methods/TCGABRCA_compare_pruneSSN_topology/", recursive = TRUE)
}


#####


# Read data
SSN_type_list <- c("ssNITEpy", "BONOBO", "CSN", "LIONESS", "SWEET", "CSNsparse", "SWEETsparse")

pruned_net_stats <- list()
for(ssn_type in SSN_type_list){
  
  pruned_net_stats[[ssn_type]] <- read.csv(paste0("IOFiles/Pruned_SSN_analysis/TCGABRCA/undir_", ssn_type, "/", ssn_type, "_prunedSSNstats__DA_TCGA__DT_BreastCancer.csv"), 
                                           header = TRUE, 
                                           check.names = FALSE)
  
}

pruned_net_stats <- rbindlist(pruned_net_stats, idcol = "SSN type")


#####


# Generate plot network topologies
plot_metrics <- c("Density", "Connected component count", "Largest component size", 
                  "Number of isolates", "Transitivity", "Assortativity", 
                  "Mean distance", "Diameter", "Modularity", 
                  "Number of modules", "Gini coefficient", "Hill index")



plot <- melt(pruned_net_stats, 
             id.vars = c("SSN type", "Sample ID"), 
             variable.name = "Metric", 
             value.name = "Score") %>% 
  filter(Metric %in% plot_metrics)

plot$`SSN type` <- factor(plot$`SSN type`, levels = c("ssNITEpy", "BONOBO", "CSN", "LIONESS", "SWEET", "CSNsparse", "SWEETsparse"))

plot <- ggplot(plot, aes(x = `SSN type`, y = Score, fill = `SSN type`)) +
  geom_boxplot(width = 0.75, 
               lwd = 0.2,
               outliers = FALSE, 
               na.rm = TRUE) +
  facet_wrap(~ Metric, scale = "free_y", nrow = 2) +
  labs(x = "SSN type",
       y = "Value") +
  scale_x_discrete(labels = c("ssNITEpy" = "ssNITE")) +
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
        strip.background = element_rect(color = "black", linewidth = 0.25,),
        strip.text = element_text(size = 5.5, margin = margin(t = 1, b = 1, r = 1, l= 1)),
        text = element_text(size = 6),
        axis.text.x = element_text(size = 4, angle = 45, vjust = 1, hjust = 1),
        axis.ticks = element_line(colour = "black", linewidth = 0.2),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.key = element_blank(),
        legend.key.spacing.x = unit(0.01, "cm"),
        legend.key.spacing.y = unit(0.01, "cm"),
        legend.key.size = unit(0.2, "cm"),
        legend.title = element_text(size = 4, margin = margin(t = 1, r = 1, b = 1, l = 1)),
        legend.text = element_text(size = 3, margin = margin(t = 1, r = 2, b = 1, l = 1)),
        legend.margin = margin(1,1,1,1),
        legend.box.spacing = unit(0.1, "cm"),
        legend.box.background = element_rect(colour = "black", linewidth = 0.25))


tiff(paste0("IOFiles/Compare_X_methods/TCGABRCA_compare_pruneSSN_topology/TCGABRCA_compare_pruneSSN_topology_plot1.tiff"),
     width = 20, height = 8,
     units = "cm", compression = "lzw", res = 1200)

plot

dev.off()


#####


# Generate plot SSN prunning stats

plot_metrics <- c("Number of nodes", "Number of edges", "Largest component density", "Selected threshold", 
                  "Used threshold",  "Mean degree", "Girth", "Reciprocity", 
                  "Power law exp.", "Power law xmin", "Rsquared", "Slope") 

plot <- melt(pruned_net_stats, 
             id.vars = c("SSN type", "Sample ID"), 
             variable.name = "Metric", 
             value.name = "Score") %>% 
  filter(Metric %in% plot_metrics)

plot$`SSN type` <- factor(plot$`SSN type`, levels = c("ssNITEpy", "BONOBO", "CSN", "LIONESS", "SWEET", "CSNsparse", "SWEETsparse"))

plot <- ggplot(plot, aes(x = `SSN type`, y = Score, fill = `SSN type`)) +
  geom_boxplot(width = 0.75, 
               lwd = 0.2, 
               outliers = FALSE, 
               na.rm = TRUE) +
  facet_wrap(~ Metric, scale = "free_y", nrow = 2) +
  labs(x = "SSN type",
       y = "Value") +
  scale_x_discrete(labels = c("ssNITEpy" = "ssNITE")) +
  scale_fill_manual(values = c("ssNITEpy" = "#33A02C", "ssNITE" = "#B2DF8A", 
                               "BONOBO" = "#1F78B4", "BONOBOsparse" = "#A6CEE3", 
                               "CSN" = "#E72F31", "CSNsparse" = "#FB9A99", 
                               "SWEET" = "#CC6600", "SWEETsparse" = "#FDBF6F", 
                               "LIONESS" = "#6A3D9A"), 
                    drop = FALSE) +  
  theme(panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.25, linetype = NULL),
        panel.grid = element_blank(),
        panel.spacing = unit(0.1, "cm"),
        strip.background = element_rect(color = "black", linewidth = 0.25,),
        strip.text = element_text(size = 5.5, margin = margin(t = 1, b = 1, r = 1, l= 1)),
        text = element_text(size = 6),
        axis.text.x = element_text(size = 4, angle = 45, vjust = 1, hjust = 1),
        axis.ticks = element_line(colour = "black", linewidth = 0.2),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.key = element_blank(),
        legend.key.spacing.x = unit(0.01, "cm"),
        legend.key.spacing.y = unit(0.01, "cm"),
        legend.key.size = unit(0.2, "cm"),
        legend.title = element_text(size = 4, margin = margin(t = 1, r = 1, b = 1, l = 1)),
        legend.text = element_text(size = 3, margin = margin(t = 1, r = 2, b = 1, l = 1)),
        legend.margin = margin(1,1,1,1),
        legend.box.spacing = unit(0.1, "cm"),
        legend.box.background = element_rect(colour = "black", linewidth = 0.25))

tiff(paste0("IOFiles/Compare_X_methods/TCGABRCA_compare_pruneSSN_topology/TCGABRCA_compare_pruneSSN_topology_plot2.tiff"),
     width = 20, height = 8,
     units = "cm", compression = "lzw", res = 1200)

plot

dev.off()


#####


print(warnings())