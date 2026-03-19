set.seed(5081)


# Script to check the correlation between SSN topological features and sample information for TCGA-BRCA samples

# Load libraries
library(unixtools)
library(optparse)
library(tidyverse)


# Set temporary directory
if(!dir.exists("tmp_dir/")){dir.create("tmp_dir/", recursive = TRUE)}
set.tempdir("tmp_dir/")


#####


# Get arguments
option_list = list(
  make_option(c("--network_type"), type = "character", default = NULL,
              help = "Network type to be used. Possible values: ssNITE, ssNITEpy, BONOBO, CSN, LIONESS, SWEET, BONOBOsparse, CSNsparse, SWEETsparse", metavar = "character")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)


if(is.null(opt$network_type)){
  print_help(opt_parser)
  stop("--network_type argument needed", call. = FALSE)
}

if(!opt$network_type %in% c("ssNITE", "ssNITEpy", "BONOBO", "CSN", "LIONESS", "SWEET", "BONOBOsparse", "CSNsparse", "SWEETsparse")){
  print_help(opt_parser)
  stop("--network_type must be one of the following: ssNITE, ssNITEpy, BONOBO, CSN, LIONESS, SWEET, BONOBOsparse, CSNsparse, SWEETsparse", call. = FALSE)
}


# Define global options for this script
network_type <- opt$network_type

print(paste0(Sys.time(), " | Using the following parameters: "))
print(paste0(Sys.time(), " | Network type: ", network_type))


#####


print(paste0(Sys.time(), " | Reading the sample information"))

TCGA_data <- readRDS("IOFiles/GDC_data/vstCounts__DA_TCGA__DT_BreastCancer.rds")
sample_info <- TCGA_data$sample_info %>% filter(condition != "SolidTissueNormal_NA") %>% droplevels()
sample_info <- sample_info %>%
  select(c("barcode", "patient", "race", "paper_pathologic_stage", "paper_BRCA_Subtype_PAM50", "paper_age_at_initial_pathologic_diagnosis"))  %>% 
  rename(c("sample_id" = "barcode", 
           "Stage" = "paper_pathologic_stage", 
           "Subtype" = "paper_BRCA_Subtype_PAM50", 
           "Age at diagnosis" = "paper_age_at_initial_pathologic_diagnosis",
           "Race" = "race"))
row.names(sample_info) <- NULL

sample_info[which(sample_info$Race == "not reported"), "Race"] <- NA
remove_race <- names(which(table(sample_info$Race) < 10))
sample_info[which(sample_info$Race %in% remove_race), "Race"] <- NA
rm(list = c("TCGA_data", "remove_race"))


#####


# Add the survival information for the patients

print(paste0(Sys.time(), " | Adding the survival data for the samples"))

if(!dir.exists("Databases/Genomic_Data_Commons/")){dir.create("Databases/Genomic_Data_Commons/", recursive = TRUE)}
if(!file.exists("Databases/Genomic_Data_Commons/TCGA-CDR-SupplementalTableS1.xlsx")){
  download.file(url = "https://api.gdc.cancer.gov/data/1b5f413e-a8d1-4d10-92eb-7c4ae739ed81",
                destfile = "Databases/Genomic_Data_Commons/TCGA-CDR-SupplementalTableS1.xlsx", 
                method = "wget")
}

TCGA_CDR <- readxl::read_xlsx("Databases/Genomic_Data_Commons/TCGA-CDR-SupplementalTableS1.xlsx", sheet = "TCGA-CDR", col_names = TRUE, guess_max = 10000)

sample_info <- sample_info %>% 
  left_join(TCGA_CDR %>% 
              select(c("bcr_patient_barcode", 
                       "OS", "OS.time", "DSS", "DSS.time", 
                       "DFI", "DFI.time", "PFI", "PFI.time")), 
            by = c("patient" = "bcr_patient_barcode")) %>% 
  rename(c("OS time" = "OS.time", 
           "DSS time" = "DSS.time", 
           "DFI time" = "DFI.time", 
           "PFI time" = "PFI.time"))
rm(TCGA_CDR)


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


# Add the tumour sample purity information

print(paste0(Sys.time(), " | Adding the tumour purity information for the samples"))

purity_data <- TCGAbiolinks::Tumor.purity %>%
  mutate(across(c(ESTIMATE, ABSOLUTE, LUMP, IHC, CPE), function(x) {
    x <- as.character(x)
    x <- str_replace(x, ",", ".")
    as.numeric(x)
  }))

sample_info <- sample_info %>% 
  mutate("tmp_id" = substr(sample_id, 1, 16)) %>% 
  left_join(purity_data %>% select(!c("Cancer.type")), 
            by = c("tmp_id" = "Sample.ID")) %>% 
  select(!c("tmp_id")) 

rm(purity_data)

sample_info <- sample_info %>% select(!c("patient"))
sample_info[sample_info == "NA"] <- NA
sample_info[sample_info == "NaN"] <- NA


#####


# Read the network-level topology info
network_topol <- read.csv(paste0("IOFiles/Pruned_SSN_analysis/TCGABRCA/undir_", network_type, "/", network_type, "_prunedSSNstats__DA_TCGA__DT_BreastCancer.csv"), check.names = FALSE)
network_topol <- network_topol %>% 
  select(!c("Number of nodes", "Number of edges", 
            "Selected threshold", "Used threshold")) %>% 
  rename("sample_id" = "Sample ID")


#####


print(paste0(Sys.time(), " | Computing the correlation between network-level topology and sample information"))

source("Scripts/Functions/correlational_analysis/Functions_SSNtopology_vs_sampleInfo.R")

corel_res <- func_SSNtopol_vs_sampleInfo_cor(SSN_topol = network_topol, 
                                             sample_info = sample_info)

IO_path <- paste0("IOFiles/Corel_SSNtopol_vs_sampleInfo/TCGABRCA/undir_", network_type, "/")
if(!dir.exists(IO_path)){dir.create(IO_path, recursive = TRUE)}

write.csv(corel_res, paste0(IO_path, network_type, "_TopoSampleCor__DA_TCGA__DT_BreastCancer.csv"), row.names = TRUE)


#####


# Plot the correlation as heatmap

tmp1 <- str_wrap(colnames(corel_res), 30)
plot_data <- corel_res %>% 
  as.data.frame() %>% 
  rownames_to_column("SSN_topol") %>%
  pivot_longer(-SSN_topol, 
               names_to = "sample_info", 
               values_to = "correlation")
plot_data$sample_info <- factor(str_wrap(plot_data$sample_info, 30), levels = tmp1)
rm(tmp1)

plot <- ggplot(plot_data, aes(x = sample_info, y = SSN_topol, fill = correlation, label = round(correlation, 2))) +
  geom_tile() +
  geom_text(size = 1, 
            na.rm = TRUE) +
  labs(fill = "Correlation", 
       x = "Sample info.", 
       y = "Topology") +
  scale_fill_gradient2(low = "blue", 
                       mid = "white", 
                       high = "red",
                       midpoint = 0,
                       limits = c(-1, 1), 
                       na.value = "#F5F5F5") +
  theme(panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.25, linetype = NULL),
        panel.grid = element_blank(),
        panel.spacing = unit(0.1, "cm"),
        strip.background = element_rect(color = "black", linewidth = 0.25,),
        strip.text = element_text(size = 5, margin = margin(t = 1, b = 1, r = 1, l= 1)),
        text = element_text(size = 4), 
        plot.title = element_text(size = 5, hjust = 0.5, face = "bold"),
        axis.text.x = element_text(size = 4, angle = 45, vjust = 1, hjust = 1), 
        axis.ticks = element_line(colour = "black", linewidth = 0.2),
        legend.position = "right",
        legend.title = element_text(size = 3),
        legend.key = element_blank(),
        legend.key.size = unit(0.2, "cm"),
        legend.text = element_text(size = 2.5, margin = margin(t = 0, b = 0, r = 0, l = 0.75)),
        legend.margin = margin(1,1,1,1),
        legend.box.spacing = unit(0.1, "cm"),
        legend.box.background = element_rect(colour = "black", linewidth = 0.1))

tiff( paste0(IO_path, network_type, "_TopoSampleCor__DA_TCGA__DT_BreastCancer.tiff"),
      width = 15, height = 7,
      units = "cm", compression = "lzw", res = 1200 )

print(plot)

dev.off() 


#####


print(warnings())